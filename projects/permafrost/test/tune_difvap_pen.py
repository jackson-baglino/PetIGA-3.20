#!/usr/bin/env python3
"""
tune_difvap_pen.py — Sweep the vapor diffusion penalty parameter difvap_pen.

For each value the script runs two short simulations:

  Run A  (flat saturated, hum=1.0, T06):
    Measures rhov conservation — the total vapour mass must not drift when
    the domain is already at saturation.  Metric: |Δtot_rhov| / tot_rhov(0).

  Run B  (sublimation, hum=0.5, T05):
    Measures the sublimation signal and SNES solver health.
    Metrics: d(tot_rhov)/dt > 0 (ice sublimes → vapour released),
             max SNES Newton iters ≤ threshold.

Sweep values (default)
  difvap_pen = 0  (no penalty — baseline)
               then 1e-10 to 1e0 in decade steps

Usage
-----
  cd <project_root>
  python test/tune_difvap_pen.py

  # Custom sweep:
  python test/tune_difvap_pen.py --sweep 0 1e-8 1e-6 1e-4 1e-2 1

  # Skip the sublimation run:
  python test/tune_difvap_pen.py --skip-runB

  # Use a non-default binary or opts:
  python test/tune_difvap_pen.py \\
      --binary ./permafrost \\
      --opts-a inputs/tests/test_T06_flat_stable.opts \\
      --opts-b inputs/tests/test_T05_sublimation.opts \\
      --out-dir test/tune_difvap
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import threading
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── project root (one level up from test/) ───────────────────────────────────
ROOT = Path(__file__).parent.parent.resolve()

# ── defaults ─────────────────────────────────────────────────────────────────
BINARY_DEFAULT = str(ROOT / "permafrost")
OPTS_A_DEFAULT = str(ROOT / "inputs" / "tests" / "test_1D_IceSlab.opts")
OPTS_B_DEFAULT = str(ROOT / "inputs" / "tests" / "test_2D_IceSlab.opts")
OUTDIR_DEFAULT = str(ROOT / "test" / "tune_difvap")
TIMEOUT_DEFAULT = 300  # s wall-clock per run
T_FINAL_DEFAULT = 86400.0  # s simulated time

SWEEP_DEFAULT = [0.0] + [10.0 ** e for e in range(-10, 1)]

# Pass/fail thresholds
PASS_RHOV_DRIFT_PCT = 1.0    # %  — Run A: max |Δtot_rhov| / tot_rhov(0)
PASS_SNES_ITERS     = 450      # —  — both runs
EPS_DEFAULT         = 9.3295e-7   # m  (used only to report k_pen)


# ── monitor line parser ───────────────────────────────────────────────────────
_NUM = r"[-+]?[eE\d.]+(?:[eE][+-]?\d+)?"
_MROW = re.compile(
    r"^\s+(\d+)\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*$",
    re.MULTILINE,
)


def _parse_monitor(stdout: str) -> list[dict]:
    rows = []
    for m in _MROW.finditer(stdout):
        rows.append({
            "step":     int(m.group(1)),
            "time":     float(m.group(2)),
            "tot_ice":  float(m.group(4)),
            "tot_rhov": float(m.group(8)),
        })
    return rows


_BOUNDS_LINE = re.compile(
    r"BOUNDS:\s+phi_ice\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
    r"\s+phi_sed\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
    r"\s+phi_air\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
)

PHASE_LO = -0.25
PHASE_HI =  1.25


def _parse_bounds_violation(stdout: str) -> bool:
    """Return True if any phase field is reported outside [PHASE_LO, PHASE_HI]."""
    for m in _BOUNDS_LINE.finditer(stdout):
        vals = [float(m.group(i)) for i in range(1, 7)]
        if any(v < PHASE_LO or v > PHASE_HI for v in vals):
            return True
    return False


def _parse_snes_iters(stdout: str) -> int:
    iters = [int(m.group(1))
             for line in stdout.splitlines()
             for m in [re.search(r"converged.*?iterations\s+(\d+)", line, re.I)]
             if m]
    return max(iters) if iters else -1


# ── post-processing scripts ───────────────────────────────────────────────────
SCRIPT_VTK     = str(ROOT / "scripts"     / "plotpermafrost.py")
SCRIPT_SCALARS = str(ROOT / "postprocess" / "plot_scalars.py")
SCRIPT_PROFS1D = str(ROOT / "postprocess" / "plot1D_profiles.py")
SCRIPT_PROFS2D = str(ROOT / "postprocess" / "plot2D_snapshot.py")


def _get_dim(opts_file: str) -> int:
    """Parse -dim N from an opts file; returns 1 if not found."""
    try:
        with open(opts_file) as f:
            for line in f:
                m = re.match(r"^\s*-dim\s+(\d+)", line)
                if m:
                    return int(m.group(1))
    except OSError:
        pass
    return 1


def _postprocess(run_dir: str, dim: int = 1, label: str = "") -> None:
    """Run the standard post-processing suite on a completed run directory."""
    if not os.path.isfile(os.path.join(run_dir, "igasol.dat")):
        print(f"  [post] No igasol.dat in {run_dir} — skipping post-processing.")
        return

    tag = f" [{label}]" if label else ""
    print(f"  [post]{tag} Running post-processing in {run_dir} (dim={dim})")

    def _sp(cmd: list[str], desc: str) -> None:
        try:
            r = subprocess.run(
                cmd, capture_output=True, text=True,
                cwd=run_dir, timeout=300,
            )
            if r.returncode != 0:
                print(f"    [post] WARNING: {desc} exited {r.returncode}")
                if r.stderr.strip():
                    print(f"      {r.stderr.strip()[:300]}")
            else:
                print(f"    [post] {desc} done.")
        except subprocess.TimeoutExpired:
            print(f"    [post] WARNING: {desc} timed out after 300 s.")
        except FileNotFoundError as e:
            print(f"    [post] WARNING: {desc} could not start: {e}")

    # 1. Convert sol_*.dat → VTK
    os.makedirs(os.path.join(run_dir, "vtkOut"), exist_ok=True)
    _sp(["python", SCRIPT_VTK, "--vtk-dir", "vtkOut"],
        "VTK conversion (plotpermafrost.py)")

    # 2. Scalar time-series from SSA_evo.dat
    if os.path.isfile(os.path.join(run_dir, "SSA_evo.dat")):
        _sp(["python", SCRIPT_SCALARS,
             "--file", "SSA_evo.dat",
             "--time-unit", "h",
             "--save", "scalars.png"],
            "Scalar time-series (plot_scalars.py)")
    else:
        print("    [post] No SSA_evo.dat — skipping scalar plot.")

    # 3. Field profiles — branched on dimensionality
    if dim == 1:
        # plot1D_profiles.py uses early-exit flags, so each mode needs its own call.
        _sp(["python", SCRIPT_PROFS1D, "--dir", ".", "--thermal"],
            "1D phase+thermal step PNGs + thermal overlay")
        _sp(["python", SCRIPT_PROFS1D, "--dir", ".", "--first-last"],
            "1D first/last comparison")
        _sp(["python", SCRIPT_PROFS1D, "--dir", ".", "--derived",
             "--save", "derived.png"],
            "1D derived quantities")
    else:
        # For 2D: plot every sol_*.dat snapshot that exists
        import glob as _glob
        sol_files = sorted(_glob.glob(os.path.join(run_dir, "sol_*.dat")))
        steps = []
        for sf in sol_files:
            base   = os.path.splitext(os.path.basename(sf))[0]
            digits = base.lstrip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_")
            if digits:
                steps.append(int(digits))
        if steps:
            for step in steps:
                _sp(["python", SCRIPT_PROFS2D,
                     "--dir", ".", "--step", str(step),
                     "--save", f"snap_{step:05d}.png"],
                    f"2D snapshot step {step}")
        else:
            print("    [post] No sol_*.dat files found for 2D plotting.")


# ── build helper ─────────────────────────────────────────────────────────────

def _build(root: Path, skip_clean: bool = False) -> None:
    """Run make clean && make in *root*, streaming output. Aborts on failure."""
    targets = (["make", "clean"], ["make"]) if not skip_clean else (["make"],)
    for cmd in targets:
        label = " ".join(cmd)
        print(f"\n  [build] {label}")
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, cwd=str(root),
        )
        for line in proc.stdout:
            sys.stdout.write(f"    {line}")
            sys.stdout.flush()
        proc.wait()
        if proc.returncode != 0:
            sys.exit(f"\n  [build] ERROR: '{label}' failed (exit {proc.returncode}). Aborting sweep.")
    print("  [build] Build successful.")


# ── run helper ────────────────────────────────────────────────────────────────

def _run(binary: str, opts_file: str, extra_flags: list[str],
         run_dir: str, timeout: int) -> tuple[str, int]:
    """Run the permafrost binary, streaming output to terminal and outp.txt."""
    cmd = [
        "mpiexec", "-np", "4", binary,
        "-options_file", os.path.abspath(opts_file),
    ] + extra_flags
    env = os.environ.copy()
    env["folder"] = run_dir
    os.makedirs(run_dir, exist_ok=True)
    outp_path = os.path.join(run_dir, "outp.txt")
    collected: list[str] = []

    try:
        with open(outp_path, "w") as f:
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, env=env, cwd=str(ROOT),
            )

            def _reader():
                for line in proc.stdout:
                    sys.stdout.write(line)
                    sys.stdout.flush()
                    f.write(line)
                    f.flush()
                    collected.append(line)

            t = threading.Thread(target=_reader, daemon=True)
            t.start()
            try:
                proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                t.join(timeout=5)
                print(f"\n    [TIMEOUT after {timeout}s]")
                return "".join(collected), -2
            t.join()

        return "".join(collected), proc.returncode

    except Exception as e:
        print(f"    [ERROR launching process: {e}]")
        return "", -3


# ── sweep engine ──────────────────────────────────────────────────────────────

def run_sweep(binary: str, opts_a: str, opts_b: str,
              sweep: list[float], out_dir: str,
              timeout: int, skip_b: bool,
              t_final: float = T_FINAL_DEFAULT,
              skip_postprocess: bool = False) -> list[dict]:
    results = []

    for val in sweep:
        label    = f"difvap_{val:.2e}" if val > 0 else "difvap_0"
        k_pen    = val / (EPS_DEFAULT ** 2) if val > 0 else 0.0
        pen_flag = ["-difvap_pen", str(val), "-k_pen", str(k_pen),
                    "-t_final", str(t_final)]

        print(f"\n{'='*60}")
        print(f"  difvap_pen = {val:.2e}   (k_pen = {k_pen:.3e})")
        print(f"{'='*60}")

        # ── Run A: flat saturated (hum = 1.0) ────────────────────────────────
        dir_a  = os.path.join(out_dir, f"A_{label}")
        dim_a  = _get_dim(opts_a)
        stdout_a, rc_a = _run(binary, opts_a, pen_flag, dir_a, timeout)
        rows_a   = _parse_monitor(stdout_a)
        iters_a  = _parse_snes_iters(stdout_a)
        if not skip_postprocess:
            _postprocess(dir_a, dim=dim_a, label=f"A difvap={val:.2e}")

        rhov_drift_pct = None
        if rows_a:
            rhov0 = rows_a[0]["tot_rhov"]
            if rhov0 != 0:
                rhov_drift_pct = max(
                    abs(r["tot_rhov"] - rhov0) / abs(rhov0) * 100.0
                    for r in rows_a
                )

        # ── Run B: sublimation (hum = 0.5) ───────────────────────────────────
        iters_b       = -1
        rhov_increase = None
        rc_b          = None

        if not skip_b:
            dir_b  = os.path.join(out_dir, f"B_{label}")
            dim_b  = _get_dim(opts_b)
            stdout_b, rc_b = _run(binary, opts_b, pen_flag, dir_b, timeout)
            rows_b  = _parse_monitor(stdout_b)
            iters_b = _parse_snes_iters(stdout_b)
            if not skip_postprocess:
                _postprocess(dir_b, dim=dim_b, label=f"B difvap={val:.2e}")

            if len(rows_b) >= 2:
                dt_run = rows_b[-1]["time"] - rows_b[0]["time"] + 1e-300
                rhov_increase = (rows_b[-1]["tot_rhov"] - rows_b[0]["tot_rhov"]) / dt_run

        # ── Metrics ──────────────────────────────────────────────────────────
        iters_all = [i for i in [iters_a, iters_b] if i >= 0]
        max_snes  = max(iters_all) if iters_all else -1

        bounds_viol_a = _parse_bounds_violation(stdout_a)
        bounds_viol_b = False if skip_b else _parse_bounds_violation(stdout_b)

        pass_drift  = rhov_drift_pct is not None and rhov_drift_pct < PASS_RHOV_DRIFT_PCT
        pass_sublim = skip_b or (rhov_increase is not None and rhov_increase > 0)
        pass_iters  = max_snes >= 0 and max_snes <= PASS_SNES_ITERS
        pass_bounds = not bounds_viol_a and not bounds_viol_b
        passed      = pass_drift and pass_sublim and pass_iters and pass_bounds

        row = {
            "difvap_pen":      val,
            "k_pen":           k_pen,
            "rc_a":            rc_a,
            "rc_b":            rc_b,
            "rhov_drift_pct":  rhov_drift_pct,
            "rhov_increase":   rhov_increase,
            "max_snes_iters":  max_snes,
            "pass_drift":      pass_drift,
            "pass_sublim":     pass_sublim,
            "pass_iters":      pass_iters,
            "pass_bounds":     pass_bounds,
            "PASS":            passed,
        }
        results.append(row)

        drift_str = f"{rhov_drift_pct:.3f}%" if rhov_drift_pct is not None else "N/A"
        incr_str  = (f"{rhov_increase:.3e} kg/m³·s" if rhov_increase is not None
                     else ("skipped" if skip_b else "N/A"))
        bounds_str = "OK" if pass_bounds else f"VIOLATION (A={bounds_viol_a}, B={bounds_viol_b})"
        print(f"  rhov drift (A): {drift_str}   threshold < {PASS_RHOV_DRIFT_PCT}%")
        print(f"  rhov incr  (B): {incr_str}")
        print(f"  max SNES iters: {max_snes}   threshold ≤ {PASS_SNES_ITERS}")
        print(f"  phase bounds:   {bounds_str}  [{PHASE_LO}, {PHASE_HI}]")
        print(f"  → {'PASS' if passed else 'FAIL'}")

    return results


# ── plotting ──────────────────────────────────────────────────────────────────

def _plot(results: list[dict], out_dir: str):
    nonzero = [r for r in results if r["difvap_pen"] > 0]
    zero    = next((r for r in results if r["difvap_pen"] == 0), None)

    if not nonzero:
        print("  No nonzero difvap_pen values — skipping plot.")
        return

    vals   = np.array([r["difvap_pen"] for r in nonzero])
    drift  = np.array([r["rhov_drift_pct"] if r["rhov_drift_pct"] is not None else np.nan
                       for r in nonzero])
    incr   = np.array([r["rhov_increase"] if r["rhov_increase"] is not None else np.nan
                       for r in nonzero])
    snes   = np.array([r["max_snes_iters"] if r["max_snes_iters"] >= 0 else np.nan
                       for r in nonzero])
    colors = ["green" if r["PASS"] else "red" for r in nonzero]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("difvap_pen sweep", fontsize=13)

    def _panel(ax, ydata, ylabel, threshold, pass_above=False):
        ax.semilogx(vals, ydata, "-", color="steelblue", lw=2, zorder=2)
        for xi, yi, c in zip(vals, ydata, colors):
            ax.plot(xi, yi, "o", color=c, ms=9, zorder=3)
        if threshold is not None:
            ax.axhline(threshold, color="gray", ls="--", lw=1.2,
                       label=f"threshold = {threshold}")
            ax.legend(fontsize=9)
        # annotate baseline (difvap_pen = 0)
        if zero is not None:
            z_val = zero.get(ylabel.split("[")[0].strip(), None)
            # use a text box instead of a point since log-scale can't show 0
            baseline_txt = ""
            if "drift" in ylabel.lower() and zero["rhov_drift_pct"] is not None:
                baseline_txt = f"difvap=0: {zero['rhov_drift_pct']:.2f}%"
            elif "snes" in ylabel.lower() and zero["max_snes_iters"] >= 0:
                baseline_txt = f"difvap=0: {zero['max_snes_iters']}"
            elif "rhov" in ylabel.lower() and zero["rhov_increase"] is not None:
                baseline_txt = f"difvap=0: {zero['rhov_increase']:.2e}"
            if baseline_txt:
                ax.text(0.03, 0.97, baseline_txt, transform=ax.transAxes,
                        fontsize=8, va="top", ha="left",
                        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", alpha=0.8))
        ax.set_xlabel("difvap_pen [m²/s]", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=9)

    _panel(axes[0], drift, r"|$\Delta\rho_{v,tot}$| / $\rho_{v,tot,0}$  [%]",
           PASS_RHOV_DRIFT_PCT)
    _panel(axes[1], incr,  r"$d\rho_{v,tot}/dt$  [kg/m³·s⁻¹]  (B: sublimation signal)",
           None)
    _panel(axes[2], snes,  "Max SNES Newton iters / step",
           PASS_SNES_ITERS)

    plt.tight_layout()
    path = os.path.join(out_dir, "sweep_difvap_pen.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Plot saved to: {path}")


# ── CSV + summary ─────────────────────────────────────────────────────────────

def _save_csv(results: list[dict], out_dir: str):
    if not results:
        return
    path = os.path.join(out_dir, "sweep_difvap_pen.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    print(f"  CSV saved to:  {path}")


def _print_summary(results: list[dict]):
    print(f"\n{'─'*88}")
    print("  SUMMARY — difvap_pen sweep")
    print(f"  {'difvap_pen':>12}  {'k_pen':>12}  {'rhov_drift%':>12}  "
          f"{'rhov_incr':>12}  {'SNES its':>8}  {'bounds':>7}  {'PASS':>6}")
    print(f"  {'─'*80}")
    for r in results:
        dstr = f"{r['rhov_drift_pct']:.3f}" if r["rhov_drift_pct"] is not None else "N/A"
        istr = (f"{r['rhov_increase']:.2e}" if r["rhov_increase"] is not None
                else "N/A")
        bstr = "✓" if r["pass_bounds"] else "✗ OOB"
        print(f"  {r['difvap_pen']:>12.2e}  {r['k_pen']:>12.3e}  {dstr:>12}  "
              f"{istr:>12}  {r['max_snes_iters']:>8}  {bstr:>7}  "
              f"{'✓' if r['PASS'] else '✗':>6}")

    passed = [r for r in results if r["PASS"]]
    if passed:
        nz_pass = [r for r in passed if r["difvap_pen"] > 0]
        if nz_pass:
            lo = min(r["difvap_pen"] for r in nz_pass)
            hi = max(r["difvap_pen"] for r in nz_pass)
            print(f"\n  Recommended range: {lo:.2e} – {hi:.2e}  "
                  f"({len(passed)}/{len(results)} pass)")
        else:
            print(f"\n  Only difvap_pen = 0 passed all criteria.")
    else:
        print("\n  No values passed all criteria — consider relaxing thresholds.")
    print(f"{'─'*78}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Sweep difvap_pen (vapor diffusion penalty) for the permafrost solver.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",    default=BINARY_DEFAULT,
                   help="Path to the permafrost binary")
    p.add_argument("--opts-a",    default=OPTS_A_DEFAULT,
                   help="Opts file for run A (flat saturated, hum=1.0)")
    p.add_argument("--opts-b",    default=OPTS_B_DEFAULT,
                   help="Opts file for run B (sublimation, hum=0.5)")
    p.add_argument("--out-dir",   default=OUTDIR_DEFAULT,
                   help="Output directory for results")
    p.add_argument("--sweep",     nargs="+", type=float, default=SWEEP_DEFAULT,
                   help="difvap_pen values to test (0 = no penalty)")
    p.add_argument("--t-final",   type=float, default=T_FINAL_DEFAULT,
                   help="Simulated end time [s] (default: 86400)")
    p.add_argument("--timeout",   type=int, default=TIMEOUT_DEFAULT,
                   help="Wall-clock timeout per run (s)")
    p.add_argument("--skip-runB", action="store_true",
                   help="Skip the sublimation run (run A only)")
    p.add_argument("--skip-postprocess", action="store_true",
                   help="Skip post-processing (VTK, scalars, 1D profiles) after each run")
    p.add_argument("--skip-build", action="store_true",
                   help="Skip 'make clean && make' before the sweep (use existing binary)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 65)
    print("  difvap_pen SWEEP")
    print(f"  binary:  {args.binary}")
    print(f"  opts A:  {args.opts_a}")
    if not args.skip_runB:
        print(f"  opts B:  {args.opts_b}")
    print(f"  out_dir: {args.out_dir}")
    print(f"  sweep:   {args.sweep}")
    print(f"  t_final: {args.t_final:.4g} s")
    print(f"  timeout: {args.timeout} s per run")
    print(f"  post-processing: {'disabled' if args.skip_postprocess else 'enabled'}")
    print(f"  build:           {'skipped' if args.skip_build else 'make clean && make'}")
    print("=" * 65)

    if not args.skip_build:
        _build(ROOT)

    results = run_sweep(
        args.binary, args.opts_a, args.opts_b,
        args.sweep, args.out_dir, args.timeout, args.skip_runB,
        t_final=args.t_final,
        skip_postprocess=args.skip_postprocess,
    )

    _print_summary(results)
    _save_csv(results, args.out_dir)
    _plot(results, args.out_dir)

    print(f"\n  All output saved to: {args.out_dir}")


if __name__ == "__main__":
    main()
