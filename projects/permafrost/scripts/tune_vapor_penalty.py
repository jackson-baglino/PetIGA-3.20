#!/usr/bin/env python3
"""
tune_vapor_penalty.py — Sweep difvap_pen and k_pen to find the optimal vapor
penalty parameters for the permafrost phase-field solver.

Two independent sweeps are run:
  A) Vary difvap_pen, hold k_pen at its default (difvap_pen/eps²).
  B) Vary k_pen directly, hold difvap_pen at its default (1e-3 m²/s).

For each parameter value three cases are evaluated:
  VP-A  Flat saturated interface (hum=1.0, flag_tIC=2)  →  rhov drift, SNES iters
  VP-B  Undersaturated slab (hum=0.5)                   →  sublimation rate, SNES iters
  VP-C  (same run as VP-B)                              →  in-ice rhov error from sol files

Usage
-----
  python scripts/tune_vapor_penalty.py \\
      --binary ./permafrost \\
      --opts-vpa inputs/tests/test_T06_flat_stable.opts \\
      --opts-vpb inputs/tests/test_1D_IceSlab.opts \\
      --eps 9.3295e-7 \\
      --out-dir /tmp/tune_vapor

  # Custom sweeps:
  python scripts/tune_vapor_penalty.py --sweep-difvap 1e-7 1e-5 1e-4 1e-3 1e-2 \\
      --sweep-kpen 1e3 1e6 1e9 1e11 ...

Pass criteria (printed per run):
  |Δtot_rhov| / tot_rhov(0)        < 1 %      (VP-A flat saturated)
  d(tot_ice)/dt < 0                             (VP-B sublimation sign)
  max SNES Newton iterations        ≤ 7
  max|rhov − ρ_vs| / ρ_vs (in ice) < 5 %      (VP-C in-ice error)
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile

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
# Physical constants
# ---------------------------------------------------------------------------
EPS_DEFAULT = 9.3295e-7      # m  (interface width)
DIFVAP_PHYS = 2.178e-5       # m²/s  (physical vapour diffusivity in air)

SWEEP_DIFVAP_DEFAULT = [1e-7, 1e-5, 1e-4, 1e-3, 1e-2]   # m²/s
SWEEP_KPEN_DEFAULT   = [1e3,  1e6,  1e9,  1e11]           # dimensionless prefactor (× 1/eps²)

PASS_RHOV_DRIFT  = 1.0    # % — VP-A: max |Δtot_rhov| / tot_rhov(0)
PASS_SNES_ITERS  = 7      # — max Newton iterations per step
PASS_ICE_ERR_PCT = 5.0    # % — VP-C: max|rhov − ρ_vs| / ρ_vs inside ice


# ---------------------------------------------------------------------------
# Saturation vapour density (matches material_properties.c: RhoVS_I)
# ---------------------------------------------------------------------------
def rho_vs(T_C):
    """Saturated vapour density [kg/m³] — matches RhoVS_I() in material_properties.c."""
    T = T_C + 273.15
    K0, K1, K2 = -0.5865e4, 0.2224e2, 0.1375e-1
    K3, K4, K5 = -0.3403e-4, 0.2697e-7, 0.6918
    Patm, bb, rho_air = 101325.0, 0.62, 1.341
    Pvs = np.exp(K0/T + K1 + K2*T + K3*T**2 + K4*T**3 + K5*np.log(T))
    return rho_air * bb * Pvs / (Patm - Pvs)


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------
_NUM = r"[-+]?[eE\d.]+(?:[eE][+-]?\d+)?"
_MONITOR_RE = re.compile(
    r"^\s+(\d+)\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*$",
    re.MULTILINE,
)


def _parse_monitor(stdout: str):
    """Return list of dicts with monitor row data."""
    rows = []
    for m in _MONITOR_RE.finditer(stdout):
        rows.append({
            "step":      int(m.group(1)),
            "time":      float(m.group(2)),
            "dt":        float(m.group(3)),
            "tot_ice":   float(m.group(4)),
            "tot_air":   float(m.group(5)),
            "tot_sed":   float(m.group(6)),
            "temp":      float(m.group(7)),
            "tot_rhov":  float(m.group(8)),
            "ia_interf": float(m.group(9)),
            "trip_junc": float(m.group(10)),
        })
    return rows


def _parse_max_snes_iters(stdout: str) -> int:
    """Extract maximum Newton iteration count from snes_converged_reason output."""
    iters = []
    for line in stdout.splitlines():
        # PETSc prints: "Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations N"
        m = re.search(r"converged.*?iterations\s+(\d+)", line, re.I)
        if m:
            iters.append(int(m.group(1)))
    return max(iters) if iters else -1


def _run(binary, opts_file, extra_flags, run_dir, timeout=300):
    """Run the permafrost binary, return (stdout, returncode). Returns ("", -2) on timeout."""
    cmd = ["mpirun", "-n", "1", binary,
           "-options_file", os.path.abspath(opts_file)] + extra_flags
    env = os.environ.copy()
    env["folder"] = run_dir   # monitoring.c reads getenv("folder") for output dir
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, env=env,
            timeout=timeout, cwd=os.path.dirname(os.path.abspath(binary)) or "."
        )
        return result.stdout + result.stderr, result.returncode
    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT after {timeout}s — marking run as N/A]")
        return "", -2


def _read_sol_rhov_ice(run_dir):
    """
    Read the final sol_*.dat file and return arrays (phi_ice, rhov).
    Requires igakit.
    """
    if not _IGAKIT:
        return None, None
    import glob
    iga_path = os.path.join(run_dir, "igasol.dat")
    if not os.path.isfile(iga_path):
        return None, None
    sol_files = sorted(glob.glob(os.path.join(run_dir, "sol_*.dat")))
    if not sol_files:
        return None, None
    nrb = PetIGA().read(iga_path)
    sol = PetIGA().read_vec(sol_files[-1], nrb)
    # igakit returns (nx, n_dof) for 1D or (nx, ny, n_dof) for 2D.
    # Flatten spatial dims so indexing is always (n_nodes, n_dof).
    if sol.ndim == 3:
        sol = sol.reshape(-1, sol.shape[-1])
    phi_ice = sol[:, 0]
    rhov    = sol[:, 2]
    return phi_ice, rhov


def _in_ice_rhov_error(phi_ice, rhov, T_C=-20.0, threshold=0.5):
    """
    Compute max|rhov − ρ_vs(T)| / ρ_vs(T) inside the ice region.
    Returns percentage error, or None if no ice found.
    """
    rvs = rho_vs(T_C)
    mask = phi_ice > threshold
    if not np.any(mask):
        return None
    err = np.max(np.abs(rhov[mask] - rvs) / rvs) * 100.0
    return err


# ---------------------------------------------------------------------------
# Single-parameter sweep engine
# ---------------------------------------------------------------------------

def run_sweep(binary, opts_vpa, opts_vpb, out_dir, values, param_name,
              eps, timeout=300, extra_vpa_flags=None, extra_vpb_flags=None):
    """
    Sweep a single penalty parameter over `values`.
    Returns a list of result dicts.
    """
    os.makedirs(out_dir, exist_ok=True)
    results = []

    for val in values:
        label = f"{param_name}={val:.2e}"
        print(f"\n{'='*60}")
        print(f"  {label}")
        print(f"{'='*60}")

        # Build extra flags depending on which parameter we're sweeping
        _base_a = list(extra_vpa_flags or [])
        _base_b = list(extra_vpb_flags or [])
        if param_name == "difvap_pen":
            extra_a = _base_a + ["-difvap_pen", str(val)]
            extra_b = _base_b + ["-difvap_pen", str(val)]
        else:  # k_pen
            extra_a = _base_a + ["-k_pen", str(val)]
            extra_b = _base_b + ["-k_pen", str(val)]

        # ---- VP-A: flat saturated (hum=1.0) --------------------------------
        run_dir_a = os.path.join(out_dir, f"VPA_{label}")
        os.makedirs(run_dir_a, exist_ok=True)
        stdout_a, rc_a = _run(binary, opts_vpa, extra_a, run_dir_a, timeout)
        rows_a  = _parse_monitor(stdout_a)
        iters_a = _parse_max_snes_iters(stdout_a)

        rhov_drift_pct = None
        if rows_a:
            rhov0 = rows_a[0]["tot_rhov"]
            if rhov0 != 0:
                drifts = [abs(r["tot_rhov"] - rhov0) / abs(rhov0) * 100 for r in rows_a]
                rhov_drift_pct = max(drifts)

        # ---- VP-B: sublimation (hum=0.5) ------------------------------------
        run_dir_b = os.path.join(out_dir, f"VPB_{label}")
        os.makedirs(run_dir_b, exist_ok=True)
        stdout_b, rc_b = _run(binary, opts_vpb, extra_b + ["-humidity", "0.5"],
                              run_dir_b, timeout)
        rows_b  = _parse_monitor(stdout_b)
        iters_b = _parse_max_snes_iters(stdout_b)

        # VP-B sublimation signal: rhov should INCREASE when starting undersaturated
        # (ice sublimes → releases vapour → tot_rhov rises).  Ice volume change is
        # unmeasurable at millisecond timescales, so we use the vapour accumulation
        # rate instead.
        rhov_increase = None
        if len(rows_b) >= 2:
            dt_run = rows_b[-1]["time"] - rows_b[0]["time"] + 1e-300
            rhov_increase = (rows_b[-1]["tot_rhov"] - rows_b[0]["tot_rhov"]) / dt_run

        # ---- VP-C: in-ice rhov error from sol files -------------------------
        T_ref = -20.0
        phi_ice, rhov = _read_sol_rhov_ice(run_dir_b)
        ice_err_pct = _in_ice_rhov_error(phi_ice, rhov, T_C=T_ref) if phi_ice is not None else None

        max_iters = max(i for i in [iters_a, iters_b] if i >= 0) if any(i >= 0 for i in [iters_a, iters_b]) else -1

        # ---- Pass/fail -------------------------------------------------------
        pass_rhov    = (rhov_drift_pct is not None) and (rhov_drift_pct < PASS_RHOV_DRIFT)
        pass_sublim  = (rhov_increase  is not None) and (rhov_increase  > 0)   # vapour rises → sublimation active
        pass_iters   = (max_iters      >= 0)         and (max_iters      <= PASS_SNES_ITERS)
        pass_ice     = (ice_err_pct    is None)       or  (ice_err_pct    < PASS_ICE_ERR_PCT)
        passed = pass_rhov and pass_sublim and pass_iters and pass_ice

        row = {
            "param":            param_name,
            "value":            val,
            "rc_a":             rc_a,
            "rc_b":             rc_b,
            "rhov_drift_pct":   rhov_drift_pct,
            "rhov_increase":    rhov_increase,
            "max_snes_iters":   max_iters,
            "ice_rhov_err_pct": ice_err_pct,
            "pass_rhov":        pass_rhov,
            "pass_sublim":      pass_sublim,
            "pass_iters":       pass_iters,
            "pass_ice":         pass_ice,
            "PASS":             passed,
        }
        results.append(row)

        verdict = "PASS" if passed else "FAIL"
        print(f"  rhov drift:  {rhov_drift_pct:.3f}%" if rhov_drift_pct is not None else "  rhov drift:  N/A")
        print(f"  rhov incr:   {rhov_increase:.3e} /s" if rhov_increase is not None else "  rhov incr:   N/A")
        print(f"  max SNES it: {max_iters}")
        print(f"  in-ice err:  {ice_err_pct:.2f}%" if ice_err_pct is not None else "  in-ice err:  N/A")
        print(f"  → {verdict}")

    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_sweep(results, out_dir, param_name, param_label):
    vals          = [r["value"] for r in results]
    rhov_drift    = [r["rhov_drift_pct"] for r in results]
    rhov_incr     = [r["rhov_increase"] if r["rhov_increase"] is not None else np.nan for r in results]
    snes          = [r["max_snes_iters"] if r["max_snes_iters"] >= 0 else np.nan for r in results]
    ice_err       = [r["ice_rhov_err_pct"] if r["ice_rhov_err_pct"] is not None else np.nan for r in results]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"Vapor penalty sweep: {param_label}", fontsize=13)

    def _plot_panel(ax, ydata, ylabel, threshold, pass_above=False):
        x = np.array(vals)
        y = np.array([v if v is not None else np.nan for v in ydata])
        colors = []
        for r in results:
            colors.append("green" if r["PASS"] else "red")
        ax.semilogx(x, y, "o-", color="steelblue", lw=2, zorder=3)
        for xi, yi, c in zip(x, y, colors):
            ax.plot(xi, yi, "o", color=c, ms=8, zorder=4)
        if threshold is not None:
            ax.axhline(threshold, color="gray", ls="--", lw=1.2,
                       label=f"threshold = {threshold}")
            ax.legend(fontsize=9)
        ax.set_xlabel(param_label, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=9)

    _plot_panel(axes[0, 0], rhov_drift,
                r"$|\Delta\rho_v|/\rho_{v,0}$  [%]",
                PASS_RHOV_DRIFT)
    _plot_panel(axes[0, 1], rhov_incr,
                r"$d(\rho_{v,tot})/dt$  [kg/m³·s⁻¹]  (sublimation signal)",
                None)
    _plot_panel(axes[1, 0], snes,
                "Max SNES Newton iters / step",
                PASS_SNES_ITERS)
    _plot_panel(axes[1, 1], ice_err,
                r"In-ice $|\rho_v - \rho_{vs}|/\rho_{vs}$  [%]",
                PASS_ICE_ERR_PCT)

    plt.tight_layout()
    path = os.path.join(out_dir, f"sweep_{param_name}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved to: {path}")


def _save_csv(results, out_dir, param_name):
    path = os.path.join(out_dir, f"sweep_{param_name}_results.csv")
    if not results:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(f"CSV saved to:  {path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Sweep vapor penalty parameters and evaluate pass/fail metrics.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",   default="./permafrost",
                   help="Path to the permafrost binary")
    p.add_argument("--opts-vpa", default="inputs/tests/test_T06_flat_stable.opts",
                   help="Opts file for VP-A (flat saturated, hum=1.0)")
    p.add_argument("--opts-vpb", default="inputs/tests/test_1D_IceSlab.opts",
                   help="Opts file for VP-B/C (sublimation run)")
    p.add_argument("--eps",      type=float, default=EPS_DEFAULT,
                   help=f"Interface width [m] (default: {EPS_DEFAULT})")
    p.add_argument("--out-dir",  default="/tmp/tune_vapor",
                   help="Output directory for results")
    p.add_argument("--sweep-difvap", nargs="+", type=float,
                   default=SWEEP_DIFVAP_DEFAULT,
                   help="Values of difvap_pen to sweep [m²/s]")
    p.add_argument("--sweep-kpen", nargs="+", type=float,
                   default=SWEEP_KPEN_DEFAULT,
                   help="Values of k_pen to sweep [m⁻²]")
    p.add_argument("--timeout",  type=int, default=300,
                   help="Per-run timeout in seconds (default: 300)")
    p.add_argument("--skip-difvap", action="store_true",
                   help="Skip the difvap_pen sweep")
    p.add_argument("--skip-kpen",   action="store_true",
                   help="Skip the k_pen sweep")
    return p.parse_args()


def _print_summary(results, param_name, param_label):
    print(f"\n{'─'*70}")
    print(f"  SUMMARY — {param_label} sweep")
    print(f"  {'Value':>12}  {'rhov drift%':>12}  {'rhov incr>0':>11}  "
          f"{'SNES its':>8}  {'ice err%':>9}  {'PASS':>6}")
    print(f"  {'─'*62}")
    for r in results:
        rstr = f"{r['rhov_drift_pct']:.3f}" if r['rhov_drift_pct'] is not None else "N/A"
        estr = f"{r['ice_rhov_err_pct']:.2f}" if r['ice_rhov_err_pct'] is not None else "N/A"
        print(f"  {r['value']:>12.2e}  {rstr:>12}  "
              f"{'YES' if r['pass_sublim'] else 'NO':>11}  "
              f"{r['max_snes_iters']:>8}  {estr:>9}  "
              f"{'✓' if r['PASS'] else '✗':>6}")
    passed = [r for r in results if r["PASS"]]
    if passed:
        lo = min(r["value"] for r in passed)
        hi = max(r["value"] for r in passed)
        print(f"\n  Recommended range:  {lo:.2e} – {hi:.2e}  ({len(passed)}/{len(results)} pass)")
    else:
        print("\n  No values passed all criteria — consider relaxing thresholds.")
    print(f"{'─'*70}")


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    if not args.skip_difvap:
        print("\n" + "="*70)
        print("  SWEEPING difvap_pen")
        print("="*70)
        res = run_sweep(
            args.binary, args.opts_vpa, args.opts_vpb,
            os.path.join(args.out_dir, "difvap_pen"),
            args.sweep_difvap, "difvap_pen", args.eps, args.timeout
        )
        _print_summary(res, "difvap_pen", "difvap_pen [m²/s]")
        _save_csv(res, args.out_dir, "difvap_pen")
        _plot_sweep(res, args.out_dir, "difvap_pen", "difvap_pen [m²/s]")

    if not args.skip_kpen:
        print("\n" + "="*70)
        print("  SWEEPING k_pen")
        print("="*70)
        res = run_sweep(
            args.binary, args.opts_vpa, args.opts_vpb,
            os.path.join(args.out_dir, "k_pen"),
            args.sweep_kpen, "k_pen", args.eps, args.timeout
        )
        _print_summary(res, "k_pen", "k_pen [m⁻²]")
        _save_csv(res, args.out_dir, "k_pen")
        _plot_sweep(res, args.out_dir, "k_pen", "k_pen [m⁻²]")


if __name__ == "__main__":
    main()
