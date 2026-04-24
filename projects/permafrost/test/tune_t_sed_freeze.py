#!/usr/bin/env python3
"""
tune_t_sed_freeze.py — Find the optimal t_sed_freeze value.

Strategy
--------
Phase 1 — Probe:
  Run with flag_sed_mode=-1 (never freeze) at a practical dt to observe the
  free 3-phase evolution.  Parse relax_monitor.dat to find:
    t_stable  : simulated time when the ice-sed interface relative change drops
                below sed_freeze_tol  (→ safe to freeze)
    t_sinter  : simulated time when tot_sed first declines by > SINTER_TOL %
                (→ Ostwald ripening / grain sintering has begun)

Phase 2 — Sweep:
  Run flag_sed_mode=1 at discrete t_sed_freeze values in the range
  [0.5·t_stable, 2·t_sinter].  For each point measure:
    · sed_drift_pct   — |Δtot_sed| / tot_sed(0)  after penalty activates
    · ice_drift_pct   — |Δtot_ice| / tot_ice(0)  (penalty must not suppress ice)
    · max_snes_iters  — solver health

Phase 3 — Plot + report saved to test/tune_tsf/

Usage
-----
  cd <project_root>
  python test/tune_t_sed_freeze.py [options]

  python test/tune_t_sed_freeze.py --binary ./permafrost \\
      --opts inputs/tests/test_1D_IceSlab.opts \\
      --probe-dt 1e-2 --probe-tfinal 120.0 \\
      --sweep 0.5 1.0 5.0 15.0 60.0 \\
      --out-dir test/tune_tsf \\
      --timeout 600
"""
from __future__ import annotations
import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── project root (one level up from test/) ───────────────────────────────────
ROOT = Path(__file__).parent.parent.resolve()

# ── defaults ─────────────────────────────────────────────────────────────────
BINARY_DEFAULT   = str(ROOT / "permafrost")
OPTS_DEFAULT     = str(ROOT / "inputs" / "tests" / "test_1D_IceSlab.opts")
OUTDIR_DEFAULT   = str(ROOT / "test" / "tune_tsf")
PROBE_DT_DEFAULT = 1.0e-2    # s  — larger than the default 1e-4 to keep run fast
PROBE_TFINAL     = 120.0     # s  — long enough to see sintering onset
SWEEP_DEFAULT    = [0.5, 1.0, 2.0, 5.0, 15.0, 60.0]   # s
TIMEOUT_DEFAULT  = 600       # s  wall-clock per run

# Criteria
SED_FREEZE_TOL  = 1.0e-3   # relative per-step change in ice_sed_interf → stable
SINTER_TOL_PCT  = 0.5      # % drop in tot_sed  → sintering has started
POST_FREEZE_DT  = 1.0e-2   # s  — timestep used after penalty activates
POST_FREEZE_DUR = 20.0     # s  — how long after freeze to run the validation

PASS_SED_DRIFT_PCT  = 0.5   # % max sediment drift
PASS_ICE_DRIFT_PCT  = 5.0   # % max ice drift
PASS_SNES_ITERS     = 8     # max Newton iters / step


# ── helpers ───────────────────────────────────────────────────────────────────

def _run(binary: str, opts_file: str, extra_flags: list[str],
         run_dir: str, timeout: int) -> tuple[str, int]:
    """Run the permafrost binary; return (stdout+stderr, returncode)."""
    cmd = [
        "mpiexec", "-np", "4", binary,
        "-options_file", os.path.abspath(opts_file),
    ] + extra_flags
    env = os.environ.copy()
    env["folder"] = run_dir
    os.makedirs(run_dir, exist_ok=True)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, env=env,
                           timeout=timeout, cwd=str(ROOT))
        return r.stdout + r.stderr, r.returncode
    except subprocess.TimeoutExpired:
        print(f"    [TIMEOUT after {timeout}s]")
        return "", -2


def _parse_relax(path: str) -> dict | None:
    """Parse relax_monitor.dat → dict of numpy arrays."""
    try:
        data = np.loadtxt(path, comments="#")
    except Exception:
        return None
    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.shape[0] == 0 or data.shape[1] < 8:
        return None
    return {
        "step":             data[:, 0].astype(int),
        "t":                data[:, 1],
        "tot_ice":          data[:, 2],
        "tot_sed":          data[:, 3],
        "ice_air_interf":   data[:, 4],
        "sed_air_interf":   data[:, 5],
        "ice_sed_interf":   data[:, 6],
        "tot_trip":         data[:, 7],
    }


_NUM  = r"[-+]?[eE\d.]+(?:[eE][+-]?\d+)?"
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
            "step":    int(m.group(1)),
            "time":    float(m.group(2)),
            "tot_ice": float(m.group(4)),
            "tot_sed": float(m.group(6)),
        })
    return rows


def _parse_snes_iters(stdout: str) -> int:
    iters = [int(m.group(1))
             for line in stdout.splitlines()
             for m in [re.search(r"converged.*?iterations\s+(\d+)", line, re.I)]
             if m]
    return max(iters) if iters else -1


# ── Phase 1: probe ────────────────────────────────────────────────────────────

def run_probe(binary, opts_file, probe_dt, probe_tfinal, out_dir, timeout) -> dict | None:
    """Free-evolution probe: flag_sed_mode=-1."""
    run_dir = os.path.join(out_dir, "probe_free_evolution")
    print(f"\n{'='*65}")
    print(f"  PHASE 1 — Free-evolution probe")
    print(f"  dt = {probe_dt:.2e} s,  t_final = {probe_tfinal:.1f} s")
    print(f"  flag_sed_mode = -1  (no freeze)")
    print(f"{'='*65}")

    extra = [
        "-flag_sed_mode", "-1",
        "-delt_t",        str(probe_dt),
        "-t_final",       str(probe_tfinal),
        "-Permafrost_output", "0",   # skip heavy file output
        "-sed_freeze_tol", str(SED_FREEZE_TOL),
    ]
    stdout, rc = _run(binary, opts_file, extra, run_dir, timeout)

    if rc not in (0,):
        print(f"  Probe run failed (rc={rc}). Check stdout in {run_dir}.")
        log = os.path.join(run_dir, "probe_stdout.txt")
        with open(log, "w") as f:
            f.write(stdout)
        print(f"  Stdout saved to {log}")
        return None

    relax_path = os.path.join(run_dir, "relax_monitor.dat")
    rd = _parse_relax(relax_path)
    if rd is None:
        print("  Could not parse relax_monitor.dat")
        return None

    t           = rd["t"]
    ice_sed     = rd["ice_sed_interf"]
    tot_sed     = rd["tot_sed"]

    # Relative per-step change in ice_sed_interf
    rel_change = np.zeros_like(t)
    for i in range(1, len(t)):
        if ice_sed[i - 1] > 0:
            rel_change[i] = abs(ice_sed[i] - ice_sed[i - 1]) / ice_sed[i - 1]
        else:
            rel_change[i] = 1.0

    # t_stable: first time rel_change drops below SED_FREEZE_TOL (skip step 0/1)
    stable_mask = (rel_change < SED_FREEZE_TOL) & (np.arange(len(t)) > 1)
    t_stable = float(t[stable_mask][0]) if np.any(stable_mask) else float(t[-1])

    # t_sinter: first time tot_sed drops by > SINTER_TOL_PCT from its initial value
    sed0 = tot_sed[0]
    if sed0 > 0:
        sinter_mask = (sed0 - tot_sed) / sed0 * 100.0 > SINTER_TOL_PCT
        t_sinter = float(t[sinter_mask][0]) if np.any(sinter_mask) else float(t[-1])
    else:
        t_sinter = float(t[-1])

    print(f"\n  Probe results:")
    print(f"    t_stable  = {t_stable:.3f} s  (ice-sed rel change < {SED_FREEZE_TOL:.0e})")
    print(f"    t_sinter  = {t_sinter:.3f} s  (tot_sed drop > {SINTER_TOL_PCT:.1f}%)")
    if t_sinter > t_stable:
        print(f"    Optimal window: [{t_stable:.2f}, {t_sinter:.2f}] s")
    else:
        print(f"    WARNING: t_sinter ≤ t_stable — sintering may precede stabilization")

    return {
        "rd":         rd,
        "t_stable":   t_stable,
        "t_sinter":   t_sinter,
        "rel_change": rel_change,
        "run_dir":    run_dir,
    }


# ── Phase 2: sweep ────────────────────────────────────────────────────────────

def run_sweep(binary, opts_file, sweep_values, out_dir, probe_dt, timeout) -> list[dict]:
    """Run flag_sed_mode=1 at each t_sed_freeze value."""
    results = []

    for tsf in sweep_values:
        t_final = tsf + POST_FREEZE_DUR
        label   = f"tsf_{tsf:.2e}s"
        run_dir = os.path.join(out_dir, label)
        print(f"\n  {'─'*55}")
        print(f"  t_sed_freeze = {tsf:.2f} s   t_final = {t_final:.1f} s")

        extra = [
            "-flag_sed_mode", "1",
            "-t_sed_freeze",  str(tsf),
            "-delt_t",        str(probe_dt),
            "-t_final",       str(t_final),
            "-Permafrost_output", "0",
            "-sed_freeze_tol", str(SED_FREEZE_TOL),
        ]
        stdout, rc = _run(binary, opts_file, extra, run_dir, timeout)

        relax_path = os.path.join(run_dir, "relax_monitor.dat")
        rd = _parse_relax(relax_path)
        rows = _parse_monitor(stdout)
        max_snes = _parse_snes_iters(stdout)

        freeze_line_found = "SEDIMENT PENALTY ACTIVATED" in stdout

        # Metrics from monitor rows
        sed_drift_pct = ice_drift_pct = None
        if rows:
            s0 = rows[0]["tot_sed"]
            i0 = rows[0]["tot_ice"]
            if s0 != 0:
                sed_drift_pct = max(abs(r["tot_sed"] - s0) / abs(s0) * 100 for r in rows)
            if i0 != 0:
                ice_drift_pct = max(abs(r["tot_ice"] - i0) / abs(i0) * 100 for r in rows)

        # Determine actual freeze time from stdout
        t_freeze_actual = None
        if rd is not None and freeze_line_found:
            # The freeze fires at the first step where the banner appears.
            # Approximate from relax monitor: last step with step <= where it fired.
            # Easier: look for the time just before freeze in relax_monitor.dat
            # The sed0 reference is snapped at freeze, so tot_sed flattens after.
            t_freeze_actual = tsf  # best estimate without parsing banner

        passed = (
            (sed_drift_pct is not None and sed_drift_pct < PASS_SED_DRIFT_PCT) and
            (ice_drift_pct is None     or ice_drift_pct < PASS_ICE_DRIFT_PCT)  and
            (max_snes >= 0             and max_snes     <= PASS_SNES_ITERS)    and
            (rc == 0)
        )

        row = {
            "t_sed_freeze":   tsf,
            "rc":             rc,
            "sed_drift_pct":  sed_drift_pct,
            "ice_drift_pct":  ice_drift_pct,
            "max_snes_iters": max_snes,
            "freeze_fired":   freeze_line_found,
            "PASS":           passed,
        }
        results.append(row)

        sstr = f"{sed_drift_pct:.4f}%" if sed_drift_pct is not None else "N/A"
        istr = f"{ice_drift_pct:.3f}%" if ice_drift_pct is not None else "N/A"
        print(f"    sed drift:   {sstr}")
        print(f"    ice drift:   {istr}")
        print(f"    max SNES it: {max_snes}")
        print(f"    freeze fired:{freeze_line_found}")
        print(f"    → {'PASS' if passed else 'FAIL'}")

    return results


# ── Plotting ─────────────────────────────────────────────────────────────────

def _plot_probe(probe: dict, out_dir: str, sweep_values: list[float]):
    rd  = probe["rd"]
    t   = rd["t"]
    rc  = probe["rel_change"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle("Free-evolution probe  (flag_sed_mode = -1)", fontsize=13)

    t_stable = probe["t_stable"]
    t_sinter = probe["t_sinter"]

    def _vlines(ax):
        ax.axvline(t_stable, color="steelblue",  ls="--", lw=1.5, label=f"t_stable = {t_stable:.2f} s")
        ax.axvline(t_sinter, color="firebrick",  ls="--", lw=1.5, label=f"t_sinter = {t_sinter:.2f} s")
        for sv in sweep_values:
            ax.axvline(sv, color="gray", ls=":", lw=0.8, alpha=0.6)

    # (0,0) ice-sed interface
    ax = axes[0, 0]
    ax.plot(t, rd["ice_sed_interf"], color="steelblue", lw=1.8)
    _vlines(ax)
    ax.set_xlabel("Simulated time [s]");  ax.set_ylabel("ice-sed interface integral")
    ax.set_title("ice-sed interface evolution");  ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    # (0,1) per-step relative change
    ax = axes[0, 1]
    ax.semilogy(t[1:], rc[1:], color="darkorange", lw=1.5)
    ax.axhline(SED_FREEZE_TOL, color="gray", ls="--", lw=1.2, label=f"tol = {SED_FREEZE_TOL:.0e}")
    _vlines(ax)
    ax.set_xlabel("Simulated time [s]");  ax.set_ylabel("Relative Δice_sed / ice_sed")
    ax.set_title("Interface stability criterion");  ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    # (1,0) tot_sed
    ax = axes[1, 0]
    sed0 = rd["tot_sed"][0] if rd["tot_sed"][0] != 0 else 1.0
    ax.plot(t, (rd["tot_sed"] - sed0) / sed0 * 100.0, color="sienna", lw=1.8)
    ax.axhline(-SINTER_TOL_PCT, color="firebrick", ls="--", lw=1.2,
               label=f"sinter threshold = −{SINTER_TOL_PCT}%")
    _vlines(ax)
    ax.set_xlabel("Simulated time [s]");  ax.set_ylabel("Δtot_sed / tot_sed(0)  [%]")
    ax.set_title("Sediment volume change (sintering)");  ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    # (1,1) tot_ice
    ax = axes[1, 1]
    ice0 = rd["tot_ice"][0] if rd["tot_ice"][0] != 0 else 1.0
    ax.plot(t, (rd["tot_ice"] - ice0) / ice0 * 100.0, color="dodgerblue", lw=1.8)
    _vlines(ax)
    ax.set_xlabel("Simulated time [s]");  ax.set_ylabel("Δtot_ice / tot_ice(0)  [%]")
    ax.set_title("Ice volume change");  ax.legend(fontsize=8);  ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(out_dir, "probe_free_evolution.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Probe plot saved to: {path}")


def _plot_sweep(results: list[dict], probe: dict, out_dir: str):
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle("t_sed_freeze sweep — sediment quality after penalty activates", fontsize=13)

    tsf_vals = np.array([r["t_sed_freeze"] for r in results])
    colors   = ["green" if r["PASS"] else "red" for r in results]

    t_stable = probe["t_stable"] if probe else None
    t_sinter = probe["t_sinter"] if probe else None

    def _vband(ax):
        if t_stable is not None and t_sinter is not None and t_sinter > t_stable:
            ax.axvspan(t_stable, t_sinter, color="lightgreen", alpha=0.25, label="optimal window")
        if t_stable is not None:
            ax.axvline(t_stable, color="steelblue", ls="--", lw=1.2, label=f"t_stable={t_stable:.2f}s")
        if t_sinter is not None:
            ax.axvline(t_sinter, color="firebrick", ls="--", lw=1.2, label=f"t_sinter={t_sinter:.2f}s")

    # Panel 1: sed drift
    ax = axes[0]
    y  = [r["sed_drift_pct"] if r["sed_drift_pct"] is not None else np.nan for r in results]
    ax.semilogx(tsf_vals, y, "-", color="steelblue", lw=1.8, zorder=2)
    for xi, yi, c in zip(tsf_vals, y, colors):
        ax.plot(xi, yi, "o", color=c, ms=9, zorder=3)
    ax.axhline(PASS_SED_DRIFT_PCT, color="gray", ls="--", lw=1.2,
               label=f"threshold {PASS_SED_DRIFT_PCT}%")
    _vband(ax)
    ax.set_xlabel("t_sed_freeze [s]");  ax.set_ylabel("|Δtot_sed| / tot_sed(0)  [%]")
    ax.set_title("Sediment drift");  ax.legend(fontsize=7);  ax.grid(True, alpha=0.3)

    # Panel 2: ice drift
    ax = axes[1]
    y  = [r["ice_drift_pct"] if r["ice_drift_pct"] is not None else np.nan for r in results]
    ax.semilogx(tsf_vals, y, "-", color="darkorange", lw=1.8, zorder=2)
    for xi, yi, c in zip(tsf_vals, y, colors):
        ax.plot(xi, yi, "o", color=c, ms=9, zorder=3)
    ax.axhline(PASS_ICE_DRIFT_PCT, color="gray", ls="--", lw=1.2,
               label=f"threshold {PASS_ICE_DRIFT_PCT}%")
    _vband(ax)
    ax.set_xlabel("t_sed_freeze [s]");  ax.set_ylabel("|Δtot_ice| / tot_ice(0)  [%]")
    ax.set_title("Ice volume change");  ax.legend(fontsize=7);  ax.grid(True, alpha=0.3)

    # Panel 3: SNES iters
    ax = axes[2]
    y  = [r["max_snes_iters"] for r in results]
    ax.semilogx(tsf_vals, y, "-", color="purple", lw=1.8, zorder=2)
    for xi, yi, c in zip(tsf_vals, y, colors):
        ax.plot(xi, yi, "o", color=c, ms=9, zorder=3)
    ax.axhline(PASS_SNES_ITERS, color="gray", ls="--", lw=1.2,
               label=f"threshold {PASS_SNES_ITERS}")
    _vband(ax)
    ax.set_xlabel("t_sed_freeze [s]");  ax.set_ylabel("Max SNES Newton iters / step")
    ax.set_title("Solver health");  ax.legend(fontsize=7);  ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(out_dir, "sweep_t_sed_freeze.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Sweep plot saved to: {path}")


def _save_csv(results: list[dict], out_dir: str):
    if not results:
        return
    path = os.path.join(out_dir, "sweep_t_sed_freeze.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    print(f"  CSV saved to:  {path}")


def _print_summary(results: list[dict], probe: dict | None):
    print(f"\n{'─'*72}")
    print("  SUMMARY — t_sed_freeze sweep")
    if probe:
        print(f"  Probe found:  t_stable = {probe['t_stable']:.3f} s,  "
              f"t_sinter = {probe['t_sinter']:.3f} s")
        if probe["t_sinter"] > probe["t_stable"]:
            print(f"  Optimal range: [{probe['t_stable']:.2f}, {probe['t_sinter']:.2f}] s")
    print(f"\n  {'t_sed_freeze':>14}  {'sed_drift%':>10}  {'ice_drift%':>10}  "
          f"{'SNES its':>8}  {'freeze?':>7}  {'PASS':>6}")
    print(f"  {'─'*65}")
    for r in results:
        sstr = f"{r['sed_drift_pct']:.4f}" if r["sed_drift_pct"] is not None else "N/A"
        istr = f"{r['ice_drift_pct']:.3f}" if r["ice_drift_pct"] is not None else "N/A"
        fstr = "yes" if r["freeze_fired"] else "no"
        print(f"  {r['t_sed_freeze']:>14.2f}  {sstr:>10}  {istr:>10}  "
              f"{r['max_snes_iters']:>8}  {fstr:>7}  "
              f"{'✓' if r['PASS'] else '✗':>6}")

    passed = [r for r in results if r["PASS"]]
    if passed:
        lo = min(r["t_sed_freeze"] for r in passed)
        hi = max(r["t_sed_freeze"] for r in passed)
        rec = lo   # conservative: freeze as soon as it's safe
        print(f"\n  Recommended t_sed_freeze range:  {lo:.2f} – {hi:.2f} s")
        print(f"  Conservative recommendation:     {rec:.2f} s  (freeze early, penalty activates)")
    else:
        print("\n  No values passed all criteria — review thresholds or run longer.")
    print(f"{'─'*72}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Find optimal t_sed_freeze for the permafrost phase-field model.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--binary",       default=BINARY_DEFAULT)
    p.add_argument("--opts",         default=OPTS_DEFAULT)
    p.add_argument("--out-dir",      default=OUTDIR_DEFAULT)
    p.add_argument("--probe-dt",     type=float, default=PROBE_DT_DEFAULT,
                   help="timestep for the free-evolution probe (s)")
    p.add_argument("--probe-tfinal", type=float, default=PROBE_TFINAL,
                   help="duration of the free-evolution probe (s)")
    p.add_argument("--sweep",        nargs="+", type=float, default=SWEEP_DEFAULT,
                   help="t_sed_freeze values to sweep (s)")
    p.add_argument("--timeout",      type=int,   default=TIMEOUT_DEFAULT,
                   help="wall-clock timeout per run (s)")
    p.add_argument("--skip-probe",   action="store_true",
                   help="skip probe run and go straight to sweep")
    p.add_argument("--skip-sweep",   action="store_true",
                   help="run probe only, skip the sweep")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 65)
    print("  t_sed_freeze TUNING SCRIPT")
    print(f"  binary:  {args.binary}")
    print(f"  opts:    {args.opts}")
    print(f"  out_dir: {args.out_dir}")
    print(f"  probe dt = {args.probe_dt:.2e} s,  probe t_final = {args.probe_tfinal:.1f} s")
    print(f"  sweep values: {args.sweep} s")
    print("=" * 65)

    probe = None

    # ── Phase 1 ──────────────────────────────────────────────────────────────
    if not args.skip_probe:
        probe = run_probe(args.binary, args.opts,
                          args.probe_dt, args.probe_tfinal,
                          args.out_dir, args.timeout)
        if probe:
            _plot_probe(probe, args.out_dir, args.sweep)

    # ── Phase 2 ──────────────────────────────────────────────────────────────
    if not args.skip_sweep:
        print(f"\n{'='*65}")
        print("  PHASE 2 — Discrete sweep")
        print(f"{'='*65}")
        results = run_sweep(args.binary, args.opts, args.sweep,
                            args.out_dir, args.probe_dt, args.timeout)
        _print_summary(results, probe)
        _save_csv(results, args.out_dir)
        if probe:
            _plot_sweep(results, probe, args.out_dir)

    print(f"\n  All output saved to: {args.out_dir}")


if __name__ == "__main__":
    main()
