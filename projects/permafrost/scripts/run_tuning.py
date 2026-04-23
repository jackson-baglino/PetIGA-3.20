#!/usr/bin/env python3
"""
run_tuning.py — Master orchestration script for penalty parameter tuning.

Runs both vapor and sediment penalty sweeps across three geometries:
  1. test_1D_IceSlab           (1-D, fast)
  2. test_2D_IceSlab           (2-D thin slab, medium)
  3. test3_EnclosedGrainPair   (2-D enclosed grains, slow)

Results are written to:
  ./test/tune_vapor/<geometry>/   — per-geometry vapor sweep CSV + PNG
  ./test/tune_sed/<geometry>/     — per-geometry sediment sweep CSV + PNG
  ./test/tune_vapor/optimal_params.json
  ./test/tune_sed/optimal_params.json

Usage (run from the project root):
  python scripts/run_tuning.py [--binary ./permafrost] [--skip-vapor] [--skip-sed]
  python scripts/run_tuning.py --geom 1D_IceSlab        # single geometry
"""

import argparse
import json
import os
import sys

# ─── locate modules ────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

import tune_vapor_penalty as tvp
import tune_sed_penalty   as tsp

# ─── constants ─────────────────────────────────────────────────────────────────
EPS = 9.3295e-7   # interface width [m]

# Step counts per geometry (number of time steps to run).
# t_final = nsteps * delt_t;  all opts files use delt_t = 1e-4 s.
NSTEPS_1D   = 50    # fast; 1-D mesh 152 nodes
NSTEPS_2D   = 25    # medium; 2-D mesh 114×23
NSTEPS_T3   = 15    # slow; 2-D mesh 143×285


def _t_final(nsteps, delt_t=1e-4):
    return nsteps * delt_t


def _p(name):
    return os.path.join(PROJECT_ROOT, name)


# ─── geometry configurations ──────────────────────────────────────────────────
#
# opts_vpa       — opts file for VP-A (flat saturated, hum≈1.0)
# opts_vpb       — opts file for VP-B/C (sublimation; run_sweep adds -humidity 0.5)
# opts_sed       — opts file for sediment sweep
# extra_vpa      — extra binary flags prepended to every VP-A invocation
# extra_vpb      — extra binary flags prepended to every VP-B invocation
#                  (note: run_sweep will ALSO append -humidity 0.5 for VP-B)
# extra_sed      — extra binary flags prepended to every sediment-sweep invocation
# timeout_vapor  — per-run timeout for vapor sweep [s]
# timeout_sed    — per-run timeout for sediment sweep [s]
#
GEOMETRIES = [
    {
        "name": "1D_IceSlab",
        # VP-A uses the dedicated flat-stable test (hum=1.0, t_final=1e-2 built-in)
        "opts_vpa": _p("inputs/tests/test_T06_flat_stable.opts"),
        "opts_vpb": _p("inputs/tests/test_1D_IceSlab.opts"),
        "opts_sed": _p("inputs/tests/test_1D_IceSlab.opts"),
        "extra_vpa": [],                                          # flat_stable self-contained
        "extra_vpb": ["-t_final", str(_t_final(NSTEPS_1D))],
        "extra_sed": ["-t_final", str(_t_final(NSTEPS_1D))],
        "timeout_vapor": 180,
        "timeout_sed":   180,
    },
    {
        "name": "2D_IceSlab",
        "opts_vpa": _p("inputs/tests/test_2D_IceSlab.opts"),
        "opts_vpb": _p("inputs/tests/test_2D_IceSlab.opts"),
        "opts_sed": _p("inputs/tests/test_2D_IceSlab.opts"),
        # VP-A: override humidity to 1.0 (flat-saturated check for 2D geometry)
        # VP-B: run_sweep appends -humidity 0.5, which overrides the base hum
        "extra_vpa": ["-t_final", str(_t_final(NSTEPS_2D)), "-humidity", "1.0"],
        "extra_vpb": ["-t_final", str(_t_final(NSTEPS_2D))],
        "extra_sed": ["-t_final", str(_t_final(NSTEPS_2D))],
        "timeout_vapor": 420,
        "timeout_sed":   420,
    },
    {
        "name": "test3_EnclosedGrainPair",
        "opts_vpa": _p("inputs/tests/test3_EnclosedGrainPair.opts"),
        "opts_vpb": _p("inputs/tests/test3_EnclosedGrainPair.opts"),
        "opts_sed": _p("inputs/tests/test3_EnclosedGrainPair.opts"),
        # Override initial_cond with "" so the code uses the generated ic_type=enclosed IC
        # (the IC file listed in the opts file does not exist on this machine).
        # test3 opts already has humidity=1.0, so VP-A works without extra hum flag.
        "extra_vpa": ["-t_final", str(_t_final(NSTEPS_T3)), "-initial_cond", ""],
        "extra_vpb": ["-t_final", str(_t_final(NSTEPS_T3)), "-initial_cond", ""],
        "extra_sed": ["-t_final", str(_t_final(NSTEPS_T3)), "-initial_cond", ""],
        "timeout_vapor": 720,
        "timeout_sed":   720,
    },
]

# ─── sweep value sets ─────────────────────────────────────────────────────────
SWEEP_DIFVAP   = [1e-7, 1e-5, 1e-4, 1e-3, 1e-2]     # m²/s
SWEEP_KPEN     = [1e3, 1e6, 1e9, 1e11]                # m⁻²
SWEEP_SED_PF   = [0.0, 1e-9, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2]  # k_sed = pf/eps²


# ─── helpers ──────────────────────────────────────────────────────────────────

def _passing_range(results_list, key="value"):
    """Return (min_val, max_val) among all dicts with PASS==True, or (None, None)."""
    passed = [r for r in results_list if r.get("PASS")]
    if not passed:
        return None, None
    vals = [r[key] for r in passed]
    return min(vals), max(vals)


def _intersect_ranges(ranges):
    """Given a list of (lo, hi) tuples (None means absent), return the intersection."""
    valid = [(lo, hi) for lo, hi in ranges if lo is not None and hi is not None]
    if not valid:
        return None, None
    lo = max(r[0] for r in valid)
    hi = min(r[1] for r in valid)
    if lo > hi:
        return None, None
    return lo, hi


def _recommended(lo, hi):
    """Pick a round recommended value near the low end of the passing range."""
    if lo is None:
        return None
    import math
    # geometric midpoint, rounded to one significant figure
    mid = math.sqrt(lo * hi) if hi else lo
    exp = math.floor(math.log10(mid))
    return round(mid / 10**exp) * 10**exp


def _run_vapor_sweep(binary, geom, out_root, sweep_difvap, sweep_kpen):
    geo_dir = os.path.join(out_root, geom["name"])
    os.makedirs(geo_dir, exist_ok=True)

    all_results = {}

    for param_name, values in [("difvap_pen", sweep_difvap), ("k_pen", sweep_kpen)]:
        print(f"\n{'#'*70}")
        print(f"# VAPOR SWEEP: {param_name}  |  geometry: {geom['name']}")
        print(f"{'#'*70}")
        sub_dir = os.path.join(geo_dir, param_name)
        results = tvp.run_sweep(
            binary        = binary,
            opts_vpa      = geom["opts_vpa"],
            opts_vpb      = geom["opts_vpb"],
            out_dir       = sub_dir,
            values        = values,
            param_name    = param_name,
            eps           = EPS,
            timeout       = geom["timeout_vapor"],
            extra_vpa_flags = geom["extra_vpa"],
            extra_vpb_flags = geom["extra_vpb"],
        )
        tvp._print_summary(results, param_name, param_name)
        tvp._save_csv(results, sub_dir, param_name)
        tvp._plot_sweep(results, sub_dir, param_name, param_name)
        all_results[param_name] = results

    return all_results


def _run_sed_sweep(binary, geom, out_root, sweep_pf):
    geo_dir = os.path.join(out_root, geom["name"])
    os.makedirs(geo_dir, exist_ok=True)

    print(f"\n{'#'*70}")
    print(f"# SED SWEEP: k_sed_pen  |  geometry: {geom['name']}")
    print(f"{'#'*70}")
    results = tsp.run_sweep(
        binary           = binary,
        opts_file        = geom["opts_sed"],
        out_dir          = geo_dir,
        prefactors       = sweep_pf,
        eps              = EPS,
        timeout          = geom["timeout_sed"],
        extra_base_flags = geom["extra_sed"],
    )
    tsp._print_summary(results)
    tsp._save_csv(results, geo_dir)
    tsp._plot_sweep(results, geo_dir)
    return results


# ─── aggregation ──────────────────────────────────────────────────────────────

def _aggregate_vapor(all_geo_results):
    """
    all_geo_results: {geo_name: {"difvap_pen": [...], "k_pen": [...]}}
    Returns optimal dict for vapor params.
    """
    optimal = {}
    for param_name in ("difvap_pen", "k_pen"):
        ranges = []
        for geo_name, pmap in all_geo_results.items():
            if param_name not in pmap:
                continue
            lo, hi = _passing_range(pmap[param_name])
            ranges.append((lo, hi))
            status = f"{lo:.2e}–{hi:.2e}" if lo is not None else "NONE"
            print(f"  {geo_name:30s}  {param_name}: {status}")
        lo, hi = _intersect_ranges(ranges)
        rec = _recommended(lo, hi)
        optimal[param_name] = {
            "min_passing":   lo,
            "max_passing":   hi,
            "recommended":   rec,
            "n_geometries":  len(ranges),
        }
        if lo is not None:
            print(f"  → consensus passing range: {lo:.2e} – {hi:.2e}  (recommended: {rec:.2e})")
        else:
            print(f"  → no consensus passing range for {param_name}")
    return optimal


def _aggregate_sed(all_geo_results):
    """
    all_geo_results: {geo_name: [...results...]}
    Returns optimal dict for sediment param.
    """
    ranges = []
    for geo_name, results in all_geo_results.items():
        lo, hi = _passing_range(results, key="prefactor")
        ranges.append((lo, hi))
        status = f"{lo:.2e}–{hi:.2e}" if lo is not None else "NONE"
        print(f"  {geo_name:30s}  k_sed_prefactor: {status}")
    lo, hi = _intersect_ranges(ranges)
    rec = _recommended(lo, hi)
    optimal = {
        "min_passing_prefactor": lo,
        "max_passing_prefactor": hi,
        "recommended_prefactor": rec,
        "recommended_k_sed_pen": rec / (EPS * EPS) if rec else None,
        "eps":                   EPS,
        "n_geometries":          len(ranges),
    }
    if lo is not None:
        print(f"  → consensus: prefactor {lo:.2e} – {hi:.2e}  (recommended: {rec:.2e})")
        print(f"               k_sed_pen = {optimal['recommended_k_sed_pen']:.3e} m⁻²")
    else:
        print("  → no consensus passing range for k_sed_pen")
    return optimal


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Master orchestration for penalty parameter tuning.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",     default=os.path.join(PROJECT_ROOT, "permafrost"),
                   help="Path to the permafrost binary")
    p.add_argument("--vapor-out",  default=os.path.join(PROJECT_ROOT, "test", "tune_vapor"),
                   help="Root output directory for vapor sweep results")
    p.add_argument("--sed-out",    default=os.path.join(PROJECT_ROOT, "test", "tune_sed"),
                   help="Root output directory for sediment sweep results")
    p.add_argument("--skip-vapor", action="store_true", help="Skip vapor penalty sweep")
    p.add_argument("--skip-sed",   action="store_true", help="Skip sediment penalty sweep")
    p.add_argument("--geom",       nargs="+",
                   help="Restrict to named geometries (e.g. --geom 1D_IceSlab)")
    return p.parse_args()


def main():
    args = parse_args()

    # Filter geometry list if --geom was supplied
    geoms = GEOMETRIES
    if args.geom:
        geoms = [g for g in GEOMETRIES if g["name"] in args.geom]
        if not geoms:
            print(f"ERROR: none of {args.geom} matched known geometries: "
                  f"{[g['name'] for g in GEOMETRIES]}")
            sys.exit(1)

    binary = args.binary
    if not os.path.isfile(binary):
        print(f"ERROR: binary not found: {binary}")
        sys.exit(1)

    os.makedirs(args.vapor_out, exist_ok=True)
    os.makedirs(args.sed_out,   exist_ok=True)

    # ── vapor sweeps ──────────────────────────────────────────────────────────
    vapor_geo_results = {}
    if not args.skip_vapor:
        print("\n" + "="*70)
        print("  PHASE 1 — VAPOR PENALTY SWEEPS")
        print("="*70)
        for geom in geoms:
            vapor_geo_results[geom["name"]] = _run_vapor_sweep(
                binary, geom, args.vapor_out, SWEEP_DIFVAP, SWEEP_KPEN
            )

        print("\n" + "="*70)
        print("  VAPOR AGGREGATE ACROSS GEOMETRIES")
        print("="*70)
        vapor_optimal = _aggregate_vapor(vapor_geo_results)
        vapor_optimal["eps"]      = EPS
        vapor_optimal["geometries"] = [g["name"] for g in geoms]
        opt_path = os.path.join(args.vapor_out, "optimal_params.json")
        with open(opt_path, "w") as f:
            json.dump(vapor_optimal, f, indent=2, default=str)
        print(f"\nVapor optimal params saved to: {opt_path}")
    else:
        print("Skipping vapor sweep (--skip-vapor).")

    # ── sediment sweeps ───────────────────────────────────────────────────────
    sed_geo_results = {}
    if not args.skip_sed:
        print("\n" + "="*70)
        print("  PHASE 2 — SEDIMENT PENALTY SWEEPS")
        print("="*70)
        for geom in geoms:
            sed_geo_results[geom["name"]] = _run_sed_sweep(
                binary, geom, args.sed_out, SWEEP_SED_PF
            )

        print("\n" + "="*70)
        print("  SEDIMENT AGGREGATE ACROSS GEOMETRIES")
        print("="*70)
        sed_optimal = _aggregate_sed(sed_geo_results)
        sed_optimal["geometries"] = [g["name"] for g in geoms]
        opt_path = os.path.join(args.sed_out, "optimal_params.json")
        with open(opt_path, "w") as f:
            json.dump(sed_optimal, f, indent=2, default=str)
        print(f"\nSediment optimal params saved to: {opt_path}")
    else:
        print("Skipping sediment sweep (--skip-sed).")

    print("\n" + "="*70)
    print("  TUNING COMPLETE")
    print("="*70)
    if not args.skip_vapor:
        print(f"  Vapor results : {args.vapor_out}/")
    if not args.skip_sed:
        print(f"  Sediment results: {args.sed_out}/")


if __name__ == "__main__":
    main()
