#!/usr/bin/env python3
"""
tune_sed_penalty.py — Sweep k_sed_pen to find the optimal sediment shape-restoring
penalty for the permafrost phase-field solver.

For each value of k_sed_pen the same 1D IceSlab run is used; metrics are
extracted from both the Monitor stdout and the sol_*.dat snapshot files.

Pass criteria:
  |Δtot_sed| / tot_sed(0)   < 0.05 %   (SP-A: sediment drift)
  max|φ_s(x,t) − φ_s0(x)|  < 0.02     (SP-B: spatial shape preservation)
  max SNES Newton iters/step ≤ 7        (SP-C: solver health)
  |Δtot_ice| / tot_ice(0)   < 5 %      (coexistence: penalty must not suppress ice)

Usage
-----
  python scripts/tune_sed_penalty.py \\
      --binary ./permafrost \\
      --opts inputs/tests/test_1D_IceSlab.opts \\
      --eps 9.3295e-7 \\
      --out-dir /tmp/tune_sed

  # Custom prefactor sweep (k_sed = prefactor / eps²):
  python scripts/tune_sed_penalty.py \\
      --sweep-prefactor 0 1e-9 1e-7 1e-5 1e-4 1e-3 1e-2
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from igakit.io import PetIGA
    _IGAKIT = True
except ImportError:
    _IGAKIT = False

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
EPS_DEFAULT = 9.3295e-7    # m

# k_sed_pen sweep values expressed as prefactors p where k_sed = p / eps²
SWEEP_PREFACTORS_DEFAULT = [0.0, 1e-9, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2]

PASS_SED_DRIFT_PCT  = 0.05   # % — max |Δtot_sed| / tot_sed(0)
PASS_SHAPE_ERR      = 0.02   # — max |φ_s(x,t) − φ_s0(x)|
PASS_SNES_ITERS     = 7      # — max Newton iterations per step
PASS_ICE_DRIFT_PCT  = 5.0    # % — max |Δtot_ice| / tot_ice(0)


# ---------------------------------------------------------------------------
# Monitor parser (same pattern as tune_vapor_penalty.py)
# ---------------------------------------------------------------------------
_NUM = r"[-+]?[eE\d.]+(?:[eE][+-]?\d+)?"
_MONITOR_RE = re.compile(
    r"^\s+(\d+)\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*$",
    re.MULTILINE,
)


def _parse_monitor(stdout):
    rows = []
    for m in _MONITOR_RE.finditer(stdout):
        rows.append({
            "step":     int(m.group(1)),
            "time":     float(m.group(2)),
            "tot_ice":  float(m.group(4)),
            "tot_sed":  float(m.group(6)),
        })
    return rows


def _parse_max_snes_iters(stdout):
    iters = []
    for line in stdout.splitlines():
        m = re.search(r"converged.*?iterations\s+(\d+)", line, re.I)
        if m:
            iters.append(int(m.group(1)))
    return max(iters) if iters else -1


# ---------------------------------------------------------------------------
# Sol file utilities
# ---------------------------------------------------------------------------

def _step_number(path):
    base   = os.path.splitext(os.path.basename(path))[0]
    digits = base.lstrip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_")
    return int(digits) if digits else 0


def _load_sed_profiles(run_dir):
    """
    Load φ_s(x) for every sol_*.dat file in run_dir.
    Returns (x_array, list_of_sed_arrays, list_of_steps) or (None, [], []).
    """
    if not _IGAKIT:
        return None, [], []
    iga_path = os.path.join(run_dir, "igasol.dat")
    if not os.path.isfile(iga_path):
        return None, [], []
    sol_files = sorted(glob.glob(os.path.join(run_dir, "sol_*.dat")),
                       key=_step_number)
    if not sol_files:
        return None, [], []

    nrb = PetIGA().read(iga_path)
    x = nrb.control[:, 0]

    sed_list  = []
    step_list = []
    for sf in sol_files:
        try:
            sol = PetIGA().read_vec(sf, nrb)
        except Exception:
            continue
        if sol.ndim < 2 or sol.shape[1] < 4:
            continue
        sed_list.append(np.clip(sol[:, 3], 0.0, 1.0))
        step_list.append(_step_number(sf))

    return x, sed_list, step_list


def _max_shape_error(sed_list):
    """max over t of max over x of |φ_s(x,t) − φ_s0(x)|."""
    if len(sed_list) < 2:
        return None
    phi_s0 = sed_list[0]
    return float(max(np.max(np.abs(s - phi_s0)) for s in sed_list[1:]))


# ---------------------------------------------------------------------------
# Run helper
# ---------------------------------------------------------------------------

def _run(binary, opts_file, extra_flags, run_dir, timeout=400):
    cmd = ["mpirun", "-n", "1", binary,
           "-options_file", os.path.abspath(opts_file)] + extra_flags
    env = os.environ.copy()
    env["folder"] = run_dir   # monitoring.c reads getenv("folder") for output dir
    result = subprocess.run(
        cmd, capture_output=True, text=True, env=env,
        timeout=timeout, cwd=os.path.dirname(os.path.abspath(binary)) or "."
    )
    return result.stdout + result.stderr, result.returncode


# ---------------------------------------------------------------------------
# Sweep engine
# ---------------------------------------------------------------------------

def run_sweep(binary, opts_file, out_dir, prefactors, eps, timeout=400,
              extra_base_flags=None):
    os.makedirs(out_dir, exist_ok=True)
    results = []

    for pf in prefactors:
        k_val   = pf / (eps * eps) if pf > 0 else 0.0
        label   = f"pf={pf:.2e}"
        print(f"\n{'='*60}")
        print(f"  k_sed_prefactor = {pf:.2e}  →  k_sed_pen = {k_val:.3e}")
        print(f"{'='*60}")

        run_dir_sp = os.path.join(out_dir, f"SP_{label}")
        os.makedirs(run_dir_sp, exist_ok=True)

        extra = list(extra_base_flags or []) + ["-k_sed_pen", str(k_val)]
        stdout, rc = _run(binary, opts_file, extra, run_dir_sp, timeout)

        rows    = _parse_monitor(stdout)
        iters   = _parse_max_snes_iters(stdout)
        _, sed_list, _ = _load_sed_profiles(run_dir_sp)

        # SP-A: sediment volume drift
        sed_drift_pct = None
        if rows:
            s0 = rows[0]["tot_sed"]
            if s0 != 0:
                drifts = [abs(r["tot_sed"] - s0) / abs(s0) * 100 for r in rows]
                sed_drift_pct = max(drifts)

        # SP-B: spatial shape preservation
        shape_err = _max_shape_error(sed_list)

        # Coexistence: ice must still evolve
        ice_drift_pct = None
        if rows:
            i0 = rows[0]["tot_ice"]
            if i0 != 0:
                drifts_i = [abs(r["tot_ice"] - i0) / abs(i0) * 100 for r in rows]
                ice_drift_pct = max(drifts_i)

        pass_sed   = (sed_drift_pct  is not None) and (sed_drift_pct  < PASS_SED_DRIFT_PCT)
        pass_shape = (shape_err      is None)     or  (shape_err      < PASS_SHAPE_ERR)
        pass_iters = (iters          >= 0)        and (iters          <= PASS_SNES_ITERS)
        pass_ice   = (ice_drift_pct  is None)     or  (ice_drift_pct  < PASS_ICE_DRIFT_PCT)
        passed     = pass_sed and pass_shape and pass_iters and pass_ice

        row = {
            "prefactor":      pf,
            "k_sed_pen":      k_val,
            "rc":             rc,
            "sed_drift_pct":  sed_drift_pct,
            "shape_err":      shape_err,
            "max_snes_iters": iters,
            "ice_drift_pct":  ice_drift_pct,
            "pass_sed":       pass_sed,
            "pass_shape":     pass_shape,
            "pass_iters":     pass_iters,
            "pass_ice":       pass_ice,
            "PASS":           passed,
        }
        results.append(row)

        verdict = "PASS" if passed else "FAIL"
        print(f"  sed drift:   {sed_drift_pct:.4f}%" if sed_drift_pct is not None else "  sed drift:   N/A")
        print(f"  shape err:   {shape_err:.4f}"      if shape_err     is not None else "  shape err:   N/A")
        print(f"  max SNES it: {iters}")
        print(f"  ice drift:   {ice_drift_pct:.3f}%" if ice_drift_pct is not None else "  ice drift:   N/A")
        print(f"  → {verdict}")

    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_sweep(results, out_dir):
    # Exclude prefactor=0 for log-scale axes (plot separately or use linear)
    nonzero = [r for r in results if r["prefactor"] > 0]
    zero    = [r for r in results if r["prefactor"] == 0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Sediment penalty sweep: k_sed_pen prefactor (k_sed = pf / eps²)", fontsize=13)

    def _panel(ax, results_list, ykey, ylabel, threshold, use_log=True):
        if not results_list:
            return
        x = np.array([r["prefactor"] for r in results_list])
        y = np.array([r[ykey] if r[ykey] is not None else np.nan for r in results_list])
        colors = ["green" if r["PASS"] else "red" for r in results_list]
        if use_log:
            ax.semilogx(x, y, "-", color="steelblue", lw=2, zorder=2)
        else:
            ax.plot(x, y, "-", color="steelblue", lw=2, zorder=2)
        for xi, yi, c in zip(x, y, colors):
            ax.plot(xi, yi, "o", color=c, ms=9, zorder=3)
        if threshold is not None:
            ax.axhline(threshold, color="gray", ls="--", lw=1.2,
                       label=f"threshold = {threshold}")
            ax.legend(fontsize=9)
        ax.set_xlabel("k_sed prefactor (k_sed = pf / eps²)", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=9)

    _panel(axes[0], nonzero, "sed_drift_pct",
           r"$|\Delta\phi_s|/\phi_{s,0}$  [%]", PASS_SED_DRIFT_PCT)
    _panel(axes[1], nonzero, "shape_err",
           r"$\max|\phi_s(x,t) - \phi_{s,0}(x)|$", PASS_SHAPE_ERR)
    _panel(axes[2], nonzero, "max_snes_iters",
           "Max SNES Newton iters / step", PASS_SNES_ITERS)

    plt.tight_layout()
    path = os.path.join(out_dir, "sweep_k_sed_pen.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved to: {path}")


def _save_csv(results, out_dir):
    if not results:
        return
    path = os.path.join(out_dir, "sweep_k_sed_pen_results.csv")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"CSV saved to:  {path}")


def _print_summary(results):
    print(f"\n{'─'*75}")
    print("  SUMMARY — k_sed_pen sweep")
    print(f"  {'Prefactor':>12}  {'sed drift%':>10}  {'shape err':>10}  "
          f"{'SNES its':>8}  {'ice drift%':>10}  {'PASS':>6}")
    print(f"  {'─'*65}")
    for r in results:
        sstr = f"{r['sed_drift_pct']:.4f}" if r["sed_drift_pct"] is not None else "N/A"
        estr = f"{r['shape_err']:.4f}"     if r["shape_err"]     is not None else "N/A"
        istr = f"{r['ice_drift_pct']:.3f}" if r["ice_drift_pct"] is not None else "N/A"
        print(f"  {r['prefactor']:>12.2e}  {sstr:>10}  {estr:>10}  "
              f"{r['max_snes_iters']:>8}  {istr:>10}  "
              f"{'✓' if r['PASS'] else '✗':>6}")

    passed = [r for r in results if r["PASS"]]
    if passed:
        lo = min(r["prefactor"] for r in passed)
        hi = max(r["prefactor"] for r in passed)
        print(f"\n  Recommended prefactor range:  {lo:.2e} – {hi:.2e}  "
              f"({len(passed)}/{len(results)} pass)")
        print(f"  (k_sed_pen = prefactor / eps²,  eps = {EPS_DEFAULT:.4e} m)")
    else:
        print("\n  No values passed all criteria — adjust thresholds or investigate.")
    print(f"{'─'*75}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Sweep k_sed_pen (sediment shape-restoring penalty).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",          default="./permafrost")
    p.add_argument("--opts",            default="inputs/tests/test_1D_IceSlab.opts")
    p.add_argument("--eps",             type=float, default=EPS_DEFAULT)
    p.add_argument("--out-dir",         default="/tmp/tune_sed")
    p.add_argument("--sweep-prefactor", nargs="+", type=float,
                   default=SWEEP_PREFACTORS_DEFAULT,
                   help="Prefactor values p where k_sed_pen = p / eps²")
    p.add_argument("--timeout",         type=int, default=400)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 70)
    print("  SWEEPING k_sed_pen")
    print(f"  eps = {args.eps:.4e} m   →   eps² = {args.eps**2:.4e} m²")
    print("=" * 70)

    results = run_sweep(
        args.binary, args.opts, args.out_dir,
        args.sweep_prefactor, args.eps, args.timeout
    )
    _print_summary(results)
    _save_csv(results, args.out_dir)
    _plot_sweep(results, args.out_dir)


if __name__ == "__main__":
    main()
