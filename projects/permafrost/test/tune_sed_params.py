#!/usr/bin/env python3
"""
tune_sed_params.py — Joint sweep of t_sed_freeze, k_sed_pen, and difvap_pen.

Opts file: test_1D_IceSlab.opts (1D, ice slab + adjacent sediment slab).

The simulation runs in flag_sed_mode=1: sediment evolves freely (3-phase) for
t < t_sed_freeze, then switches to the penalty formulation at t_sed_freeze.

Root cause of spurious air: if the sediment interface has not equilibrated
before the penalty activates, the unresolved mismatch drives a parasitic air
phase at the ice-sed boundary. This shows up as:
  - Δtot_air    = tot_air(final) - tot_air(at t_sed_freeze)  → should ≈ 0
  - ice_sed_drop = relative drop in ice_sed_interf after freeze → should ≈ 0

Pass criteria:
  1. Clean exit (rc == 0); non-zero → blow-up or phase-bounds violation.
  2. |Δtot_air| < --pass-air-drift   (default 5e-6, ~5 % of 1D domain air)
  3. ice_sed_drop < --pass-icesd-drop (default 20 %)

Usage
-----
  cd <project_root>
  python test/tune_sed_params.py

  # Custom sweep:
  python test/tune_sed_params.py \\
      --t-freeze 10 50 100 300 900 \\
      --k-prefactor 1e-5 1e-4 1e-3 \\
      --difvap-pen 0 1e-5 1e-3 \\
      --t-final 3600 --timeout 600

  # Skip build and post-processing for speed:
  python test/tune_sed_params.py --skip-build
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import threading
from itertools import product
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent.parent.resolve()

# ── defaults ──────────────────────────────────────────────────────────────────
BINARY_DEFAULT  = str(ROOT / "permafrost")
OPTS_DEFAULT    = str(ROOT / "inputs" / "tests" / "test_1D_IceSlab.opts")
OUTDIR_DEFAULT  = str(ROOT / "test" / "tune_sed_params")
TIMEOUT_DEFAULT = 600    # s wall-clock per run
T_FINAL_DEFAULT = 3600.0 # s simulated time (shorter than opts default for sweep speed)

EPS_DEFAULT = 9.3295e-7  # m — interface width (for k_sed_pen = prefactor / eps²)

# Default sweep values
T_SED_FREEZE_DEFAULT    = [10.0, 50.0, 100.0, 300.0, 900.0]   # s
K_SED_PREFACTOR_DEFAULT = [1e-5, 1e-4, 1e-3]                   # k_sed_pen = val / eps²
DIFVAP_PEN_DEFAULT      = [0.0, 1e-5, 1e-3]                    # m²/s

# Pass/fail thresholds (configurable via CLI)
PASS_AIR_DRIFT_DEFAULT   = 5e-6   # absolute Δtot_air (m in 1D) after freeze
PASS_ICESD_DROP_DEFAULT  = 20.0   # % drop in ice_sed_interf after freeze


# ── monitor line parser ───────────────────────────────────────────────────────
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"

# Matches the 10-column data row printed by monitoring.c
_MROW = re.compile(
    r"^\s+(\d+)\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*$",
    re.MULTILINE,
)

# Matches the BOUNDS line printed by monitoring.c
_BOUNDS_LINE = re.compile(
    r"BOUNDS:\s+phi_ice\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
    r"\s+phi_sed\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
    r"\s+phi_air\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
)
PHASE_LO, PHASE_HI = -0.25, 1.25


def _parse_monitor_rows(stdout: str) -> list[dict]:
    """Parse every 10-column monitor table row from stdout."""
    rows = []
    for m in _MROW.finditer(stdout):
        rows.append({
            "step":         int(m.group(1)),
            "time":         float(m.group(2)),
            "tot_ice":      float(m.group(4)),
            "tot_air":      float(m.group(5)),
            "tot_sed":      float(m.group(6)),
            "ice_air_int":  float(m.group(9)),
        })
    return rows


def _parse_bounds_violation(stdout: str) -> bool:
    """Return True if any phase field is reported outside [PHASE_LO, PHASE_HI]."""
    for m in _BOUNDS_LINE.finditer(stdout):
        vals = [float(m.group(i)) for i in range(1, 7)]
        if any(v < PHASE_LO or v > PHASE_HI for v in vals):
            return True
    return False


def _parse_relax_monitor(run_dir: str) -> list[dict]:
    """Read relax_monitor.dat written by monitoring.c each step."""
    path = os.path.join(run_dir, "relax_monitor.dat")
    rows = []
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    rows.append({
                        "step":          int(parts[0]),
                        "t":             float(parts[1]),
                        "ice_sed_interf": float(parts[6]),
                    })
    except OSError:
        pass
    return rows


def _compute_metrics(rows: list[dict], relax_rows: list[dict],
                     t_sed_freeze: float) -> dict:
    """
    Compute Δtot_air and ice_sed_interf stability after t_sed_freeze.

    Both metrics are evaluated relative to the first step where t >= t_sed_freeze.
    """
    if not rows:
        return {"delta_tot_air": None, "ice_sed_drop_pct": None}

    # Index of first step at or after t_sed_freeze
    freeze_idx = next(
        (i for i, r in enumerate(rows) if r["time"] >= t_sed_freeze), None
    )
    if freeze_idx is None:
        return {"delta_tot_air": None, "ice_sed_drop_pct": None,
                "note": "ended before t_sed_freeze"}

    tot_air_at_freeze = rows[freeze_idx]["tot_air"]
    tot_air_final     = rows[-1]["tot_air"]
    delta_tot_air     = tot_air_final - tot_air_at_freeze

    # ice_sed_interf drop from relax_monitor.dat
    ice_sed_drop_pct = None
    if relax_rows:
        freeze_relax_idx = next(
            (i for i, r in enumerate(relax_rows) if r["t"] >= t_sed_freeze), None
        )
        if freeze_relax_idx is not None:
            isi0 = relax_rows[freeze_relax_idx]["ice_sed_interf"]
            isif = relax_rows[-1]["ice_sed_interf"]
            if abs(isi0) > 1e-20:
                # Positive = dropped (bad); negative = grew (unlikely but OK)
                ice_sed_drop_pct = (isi0 - isif) / abs(isi0) * 100.0

    return {
        "delta_tot_air":     delta_tot_air,
        "tot_air_at_freeze": tot_air_at_freeze,
        "tot_air_final":     tot_air_final,
        "ice_sed_drop_pct":  ice_sed_drop_pct,
    }


# ── build helper ──────────────────────────────────────────────────────────────
def _build(root: Path) -> None:
    """Run make clean && make, streaming output. Aborts on failure."""
    for cmd in (["make", "clean"], ["make"]):
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
            sys.exit(f"\n  [build] ERROR: '{label}' failed (exit {proc.returncode}). Aborting.")
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
                    sys.stdout.write(line); sys.stdout.flush()
                    f.write(line); f.flush()
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
def run_sweep(binary: str, opts_file: str,
              t_freeze_vals: list[float],
              k_prefactor_vals: list[float],
              difvap_vals: list[float],
              out_dir: str, timeout: int, t_final: float,
              eps: float,
              pass_air_drift: float,
              pass_icesd_drop: float) -> list[dict]:

    results = []
    total  = len(t_freeze_vals) * len(k_prefactor_vals) * len(difvap_vals)
    run_no = 0

    for t_freeze, k_pref, difvap in product(t_freeze_vals, k_prefactor_vals, difvap_vals):
        run_no += 1
        k_sed_pen = k_pref / (eps * eps)
        label = (f"tsf{t_freeze:.0f}_ksp{k_pref:.1e}_dvp{difvap:.1e}"
                 .replace("+", ""))

        print(f"\n{'='*72}")
        print(f"  Run {run_no}/{total}: {label}")
        print(f"    t_sed_freeze={t_freeze:.0f} s   "
              f"k_pref={k_pref:.2e}  (k_sed_pen={k_sed_pen:.3e})   "
              f"difvap_pen={difvap:.2e}")
        print(f"{'='*72}")

        extra_flags = [
            "-t_sed_freeze", str(t_freeze),
            "-k_sed_pen",    str(k_sed_pen),
            "-difvap_pen",   str(difvap),
            "-t_final",      str(t_final),
            "-flag_sed_mode", "1",
        ]
        run_dir = os.path.join(out_dir, label)
        stdout, rc = _run(binary, opts_file, extra_flags, run_dir, timeout)

        rows        = _parse_monitor_rows(stdout)
        relax_rows  = _parse_relax_monitor(run_dir)
        bounds_viol = _parse_bounds_violation(stdout)
        metrics     = _compute_metrics(rows, relax_rows, t_freeze)

        delta_tot_air    = metrics.get("delta_tot_air")
        ice_sed_drop_pct = metrics.get("ice_sed_drop_pct")

        pass_exit    = (rc == 0)
        pass_bounds  = not bounds_viol
        pass_air     = (delta_tot_air is not None
                        and abs(delta_tot_air) < pass_air_drift)
        pass_icesd   = (ice_sed_drop_pct is None   # no sediment data — skip
                        or ice_sed_drop_pct < pass_icesd_drop)
        passed       = pass_exit and pass_bounds and pass_air and pass_icesd

        row = {
            "t_sed_freeze":      t_freeze,
            "k_sed_prefactor":   k_pref,
            "k_sed_pen":         k_sed_pen,
            "difvap_pen":        difvap,
            "rc":                rc,
            "delta_tot_air":     delta_tot_air,
            "ice_sed_drop_pct":  ice_sed_drop_pct,
            "pass_exit":         pass_exit,
            "pass_bounds":       pass_bounds,
            "pass_air":          pass_air,
            "pass_icesd":        pass_icesd,
            "PASS":              passed,
        }
        results.append(row)

        dstr  = f"{delta_tot_air:.4e}"   if delta_tot_air    is not None else "N/A"
        idstr = f"{ice_sed_drop_pct:.1f}%" if ice_sed_drop_pct is not None else "N/A"
        print(f"  rc={rc}  Δtot_air={dstr}  ice_sed_drop={idstr}")
        print(f"  bounds={'OK' if pass_bounds else 'VIOLATION'}  "
              f"→ {'PASS' if passed else 'FAIL'}")

    return results


# ── plotting ──────────────────────────────────────────────────────────────────
def _make_index(vals: list) -> dict:
    return {v: i for i, v in enumerate(vals)}


def _plot(results: list[dict], out_dir: str,
          t_freeze_vals: list[float],
          k_prefactor_vals: list[float],
          difvap_vals: list[float],
          pass_air_drift: float) -> None:
    """
    Per difvap_pen: heatmap of |Δtot_air| as function of
    (t_sed_freeze × k_sed_prefactor) with PASS/FAIL overlay.
    Plus a summary heatmap of the best |Δtot_air| over all difvap values.
    """
    ti = _make_index(t_freeze_vals)
    ki = _make_index(k_prefactor_vals)

    nrows = len(t_freeze_vals)
    ncols = len(k_prefactor_vals)

    cmap = "RdYlGn_r"
    vmax = pass_air_drift * 4

    # ── one heatmap per difvap_pen ──────────────────────────────────────────
    for difvap in difvap_vals:
        sub = [r for r in results if r["difvap_pen"] == difvap]
        if not sub:
            continue

        air_grid  = np.full((nrows, ncols), np.nan)
        pass_grid = np.zeros((nrows, ncols), dtype=bool)

        for r in sub:
            i = ti.get(r["t_sed_freeze"])
            j = ki.get(r["k_sed_prefactor"])
            if i is None or j is None:
                continue
            if r["delta_tot_air"] is not None:
                air_grid[i, j] = abs(r["delta_tot_air"])
            pass_grid[i, j] = r["PASS"]

        fig, ax = plt.subplots(figsize=(max(5, ncols * 1.8), max(4, nrows * 1.2)))
        im = ax.imshow(air_grid, aspect="auto", cmap=cmap,
                       vmin=0, vmax=vmax, origin="upper")
        plt.colorbar(im, ax=ax, label="|Δtot_air|  [m in 1-D]")
        ax.axhline(-0.5, color="none")  # force frame

        ax.set_xticks(range(ncols))
        ax.set_xticklabels([f"{k:.1e}" for k in k_prefactor_vals], rotation=30, ha="right")
        ax.set_yticks(range(nrows))
        ax.set_yticklabels([f"{t:.0f}" for t in t_freeze_vals])
        ax.set_xlabel("k_sed_pen prefactor  (k_sed_pen = prefactor / ε²)", fontsize=10)
        ax.set_ylabel("t_sed_freeze  [s]", fontsize=10)
        ax.set_title(f"difvap_pen = {difvap:.1e}  —  |Δtot_air| after freeze", fontsize=11)

        for i in range(nrows):
            for j in range(ncols):
                val = air_grid[i, j]
                marker = "✓" if pass_grid[i, j] else "✗"
                txt = f"{val:.2e}\n{marker}" if not np.isnan(val) else "N/A"
                ax.text(j, i, txt, ha="center", va="center",
                        fontsize=9, color="black")

        ax.axhline(y=pass_air_drift, color="none")  # invisible, just forces layout
        plt.tight_layout()
        fname = os.path.join(out_dir, f"heatmap_dvp{difvap:.1e}.png".replace("+", ""))
        fig.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Plot: {fname}")

    # ── summary: best |Δtot_air| over all difvap_pen ──────────────────────
    best_air  = np.full((nrows, ncols), np.nan)
    any_pass  = np.zeros((nrows, ncols), dtype=bool)

    for r in results:
        i = ti.get(r["t_sed_freeze"])
        j = ki.get(r["k_sed_prefactor"])
        if i is None or j is None:
            continue
        val = r["delta_tot_air"]
        if val is not None:
            aval = abs(val)
            if np.isnan(best_air[i, j]) or aval < best_air[i, j]:
                best_air[i, j] = aval
        if r["PASS"]:
            any_pass[i, j] = True

    fig, ax = plt.subplots(figsize=(max(5, ncols * 1.8), max(4, nrows * 1.2)))
    im = ax.imshow(best_air, aspect="auto", cmap=cmap,
                   vmin=0, vmax=vmax, origin="upper")
    plt.colorbar(im, ax=ax, label="Best |Δtot_air| across all difvap_pen  [m in 1-D]")
    ax.set_xticks(range(ncols))
    ax.set_xticklabels([f"{k:.1e}" for k in k_prefactor_vals], rotation=30, ha="right")
    ax.set_yticks(range(nrows))
    ax.set_yticklabels([f"{t:.0f}" for t in t_freeze_vals])
    ax.set_xlabel("k_sed_pen prefactor", fontsize=10)
    ax.set_ylabel("t_sed_freeze  [s]", fontsize=10)
    ax.set_title("Best |Δtot_air| over all difvap_pen values", fontsize=11)
    for i in range(nrows):
        for j in range(ncols):
            val = best_air[i, j]
            marker = "✓" if any_pass[i, j] else "✗"
            txt = f"{val:.2e}\n{marker}" if not np.isnan(val) else "N/A"
            ax.text(j, i, txt, ha="center", va="center", fontsize=9, color="black")
    plt.tight_layout()
    fname = os.path.join(out_dir, "heatmap_summary.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Summary plot: {fname}")


def _plot_ice_sed(results: list[dict], out_dir: str,
                  t_freeze_vals: list[float],
                  k_prefactor_vals: list[float],
                  difvap_vals: list[float],
                  pass_icesd_drop: float) -> None:
    """Heatmap of ice_sed_drop_pct — companion to the main air drift plots."""
    ti = _make_index(t_freeze_vals)
    ki = _make_index(k_prefactor_vals)
    nrows, ncols = len(t_freeze_vals), len(k_prefactor_vals)

    best_drop = np.full((nrows, ncols), np.nan)
    any_pass  = np.zeros((nrows, ncols), dtype=bool)

    for r in results:
        i = ti.get(r["t_sed_freeze"])
        j = ki.get(r["k_sed_prefactor"])
        if i is None or j is None:
            continue
        val = r["ice_sed_drop_pct"]
        if val is not None:
            if np.isnan(best_drop[i, j]) or val < best_drop[i, j]:
                best_drop[i, j] = val
        if r["PASS"]:
            any_pass[i, j] = True

    fig, ax = plt.subplots(figsize=(max(5, ncols * 1.8), max(4, nrows * 1.2)))
    im = ax.imshow(best_drop, aspect="auto", cmap="RdYlGn_r",
                   vmin=0, vmax=pass_icesd_drop * 2, origin="upper")
    plt.colorbar(im, ax=ax, label="ice_sed_interf drop after freeze  [%]")
    ax.set_xticks(range(ncols))
    ax.set_xticklabels([f"{k:.1e}" for k in k_prefactor_vals], rotation=30, ha="right")
    ax.set_yticks(range(nrows))
    ax.set_yticklabels([f"{t:.0f}" for t in t_freeze_vals])
    ax.set_xlabel("k_sed_pen prefactor", fontsize=10)
    ax.set_ylabel("t_sed_freeze  [s]", fontsize=10)
    ax.set_title("Best ice_sed_interf drop % over all difvap_pen values", fontsize=11)
    for i in range(nrows):
        for j in range(ncols):
            val = best_drop[i, j]
            marker = "✓" if any_pass[i, j] else "✗"
            txt = f"{val:.1f}%\n{marker}" if not np.isnan(val) else "N/A"
            ax.text(j, i, txt, ha="center", va="center", fontsize=9, color="black")
    plt.tight_layout()
    fname = os.path.join(out_dir, "heatmap_icesd_drop.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Ice-sed drop plot: {fname}")


# ── CSV + summary ─────────────────────────────────────────────────────────────
def _save_csv(results: list[dict], out_dir: str) -> None:
    if not results:
        return
    path = os.path.join(out_dir, "sweep_sed_params.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    print(f"  CSV: {path}")


def _print_summary(results: list[dict], pass_air_drift: float,
                   pass_icesd_drop: float) -> None:
    print(f"\n{'─'*95}")
    print("  SUMMARY — sediment parameter sweep")
    print(f"  {'t_freeze':>9}  {'k_pref':>9}  {'difvap':>9}  "
          f"{'Δtot_air':>11}  {'icesd_drop':>11}  "
          f"{'bnd':>4}  {'air':>4}  {'ics':>4}  {'PASS':>5}")
    print(f"  {'─'*85}")
    for r in results:
        dstr  = f"{r['delta_tot_air']:.3e}"   if r["delta_tot_air"]    is not None else "N/A"
        idstr = f"{r['ice_sed_drop_pct']:.1f}%" if r["ice_sed_drop_pct"] is not None else "N/A"
        print(f"  {r['t_sed_freeze']:>9.0f}  {r['k_sed_prefactor']:>9.2e}  "
              f"{r['difvap_pen']:>9.2e}  {dstr:>11}  {idstr:>11}  "
              f"{'✓' if r['pass_bounds'] else '✗':>4}  "
              f"{'✓' if r['pass_air']   else '✗':>4}  "
              f"{'✓' if r['pass_icesd'] else '✗':>4}  "
              f"{'✓' if r['PASS']       else '✗':>5}")

    passed = [r for r in results if r["PASS"]]
    if passed:
        # Prefer smallest t_freeze that passes, then smallest |Δtot_air|
        best_t = min(passed, key=lambda r: (r["t_sed_freeze"], abs(r["delta_tot_air"] or 1e30)))
        print(f"\n  Best (smallest t_freeze that passes):")
        print(f"    t_sed_freeze={best_t['t_sed_freeze']:.0f} s  "
              f"k_pref={best_t['k_sed_prefactor']:.2e}  "
              f"difvap={best_t['difvap_pen']:.2e}  "
              f"Δtot_air={best_t['delta_tot_air']:.3e}")
        print(f"  {len(passed)}/{len(results)} combinations passed "
              f"(|Δtot_air|<{pass_air_drift:.1e}, ice_sed_drop<{pass_icesd_drop:.0f}%).")
    else:
        print("\n  No combinations passed. Consider relaxing --pass-air-drift or --pass-icesd-drop.")
    print(f"{'─'*95}")


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Sweep t_sed_freeze, k_sed_pen, and difvap_pen for the permafrost solver.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",          default=BINARY_DEFAULT)
    p.add_argument("--opts",            default=OPTS_DEFAULT,
                   help="Opts file (default: test_1D_IceSlab.opts)")
    p.add_argument("--out-dir",         default=OUTDIR_DEFAULT)
    p.add_argument("--t-freeze",        nargs="+", type=float,
                   default=T_SED_FREEZE_DEFAULT,
                   metavar="T",
                   help="t_sed_freeze values [s] (default: 10 50 100 300 900)")
    p.add_argument("--k-prefactor",     nargs="+", type=float,
                   default=K_SED_PREFACTOR_DEFAULT,
                   metavar="K",
                   help="k_sed_pen prefactor values; k_sed_pen = val/eps² "
                        "(default: 1e-5 1e-4 1e-3)")
    p.add_argument("--difvap-pen",      nargs="+", type=float,
                   default=DIFVAP_PEN_DEFAULT,
                   metavar="D",
                   help="difvap_pen values [m²/s] (default: 0 1e-5 1e-3)")
    p.add_argument("--t-final",         type=float, default=T_FINAL_DEFAULT,
                   help="Simulated end time [s] (default: 3600)")
    p.add_argument("--timeout",         type=int,   default=TIMEOUT_DEFAULT,
                   help="Wall-clock timeout per run [s] (default: 600)")
    p.add_argument("--eps",             type=float, default=EPS_DEFAULT,
                   help="Interface width ε [m] for k_sed_pen conversion")
    p.add_argument("--pass-air-drift",  type=float, default=PASS_AIR_DRIFT_DEFAULT,
                   help="Max |Δtot_air| to pass (default: 5e-6)")
    p.add_argument("--pass-icesd-drop", type=float, default=PASS_ICESD_DROP_DEFAULT,
                   help="Max ice_sed_interf drop %% to pass (default: 20)")
    p.add_argument("--skip-build",      action="store_true",
                   help="Skip 'make clean && make'")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    total = len(args.t_freeze) * len(args.k_prefactor) * len(args.difvap_pen)
    print("=" * 72)
    print("  SEDIMENT PARAMETER SWEEP")
    print(f"  binary:          {args.binary}")
    print(f"  opts:            {args.opts}")
    print(f"  out_dir:         {args.out_dir}")
    print(f"  t_sed_freeze:    {args.t_freeze}")
    print(f"  k_pref:          {args.k_prefactor}")
    print(f"  difvap_pen:      {args.difvap_pen}")
    print(f"  t_final:         {args.t_final:.4g} s")
    print(f"  timeout:         {args.timeout} s/run")
    print(f"  total runs:      {total}")
    print(f"  pass_air_drift:  {args.pass_air_drift:.1e}")
    print(f"  pass_icesd_drop: {args.pass_icesd_drop:.0f} %")
    print(f"  build:           {'skipped' if args.skip_build else 'make clean && make'}")
    print("=" * 72)

    if not args.skip_build:
        _build(ROOT)

    results = run_sweep(
        args.binary, args.opts,
        args.t_freeze, args.k_prefactor, args.difvap_pen,
        args.out_dir, args.timeout, args.t_final,
        eps=args.eps,
        pass_air_drift=args.pass_air_drift,
        pass_icesd_drop=args.pass_icesd_drop,
    )

    _print_summary(results, args.pass_air_drift, args.pass_icesd_drop)
    _save_csv(results, args.out_dir)
    _plot(results, args.out_dir,
          args.t_freeze, args.k_prefactor, args.difvap_pen,
          args.pass_air_drift)
    _plot_ice_sed(results, args.out_dir,
                  args.t_freeze, args.k_prefactor, args.difvap_pen,
                  args.pass_icesd_drop)

    print(f"\n  All output in: {args.out_dir}")


if __name__ == "__main__":
    main()
