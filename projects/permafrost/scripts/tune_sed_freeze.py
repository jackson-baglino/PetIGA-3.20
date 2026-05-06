#!/usr/bin/env python3
"""
tune_sed_freeze.py — 2-D sweep over t_sed_freeze × k_pen.

Physics motivation
------------------
When sediment is frozen (phi_sed stationary) the 2-phase ice Allen-Cahn can
generate spurious air at the ice-sediment boundary if the interface was not
adequately relaxed during the 3-phase period.  A longer 3-phase window (larger
t_sed_freeze) gives more time for the diffuse interface to settle, but delays
the freeze and slows the overall simulation.  This script finds the smallest
t_sed_freeze that avoids:

  1. Spurious air generation at the ice-sed boundary
     → measured via ice_sed_interf in relax_monitor.dat (should stay stable)
  2. Unphysical total-air growth (phase blow-up)
     → measured via tot_air in outp.txt
  3. Poor solver convergence (stagnant or shrinking time steps)
     → measured via dt history in outp.txt

Metrics per (t_sed_freeze, k_pen) run
--------------------------------------
  dt_final         Final adaptive time step [s]                ↑ larger is better
  dt_growth_frac   Fraction of steps where dt increases         > DT_GROW_THRESH
  ice_sed_chg_pct  % change in ice-sed interface after freeze   < ISED_THRESH
  air_growth_pct   % growth of total air from step 0 to end     < AIR_THRESH
  max_snes_iters   Max Newton iterations per time step (NL)     ≤ SNES_THRESH
  phase_ok         No phase bound violation / crash / abort      must be True
  PASS             All criteria met simultaneously

Usage
-----
  python scripts/tune_sed_freeze.py \\
      --binary ./permafrost \\
      --opts   inputs/tests/test_1D_IceSlab.opts \\
      --out-dir /tmp/tune_tsf

  # Custom sweep grids:
  python scripts/tune_sed_freeze.py \\
      --sweep-tsf  0 1e-4 1e-3 1e-2 0.1 1.0 \\
      --sweep-kpen 1e7 1e8 1e9 1e10 1e11 \\
      --t-final 5.0
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from itertools import product
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Pass/fail thresholds  (tuneable via CLI)
# ---------------------------------------------------------------------------
DT_GROW_THRESH  = 0.5    # fraction of steps where dt must grow
ISED_THRESH     = 5.0    # % — max allowed decrease in ice-sed interface post-freeze
AIR_THRESH      = 10.0   # % — max allowed growth of total air volume
SNES_THRESH     = 7      # max Newton iters per nonlinear solve

# ---------------------------------------------------------------------------
# Default sweep grids
# ---------------------------------------------------------------------------
SWEEP_TSF_DEFAULT  = [0.0, 1e-4, 1e-3, 1e-2, 0.1]  # t_sed_freeze [s]
SWEEP_KPEN_DEFAULT = [1e7, 1e8, 1e9, 1e10, 1e11]    # k_pen


# ---------------------------------------------------------------------------
# Monitor-table parser  (outp.txt — pipe-delimited rows from Monitor())
# ---------------------------------------------------------------------------
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
_MONITOR_ROW_RE = re.compile(
    r"^\s*(\d+)\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"   # step | t | dt
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"  # ice|air|sed
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"    # temp | rhov
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*$",    # ia_interf | trip_junc
    re.MULTILINE,
)


def _parse_outp(path: str) -> list[dict]:
    """Parse the pipe-delimited monitor table written to outp.txt."""
    if not os.path.isfile(path):
        return []
    text = Path(path).read_text()
    rows = []
    for m in _MONITOR_ROW_RE.finditer(text):
        rows.append({
            "step":      int(m.group(1)),
            "t":         float(m.group(2)),
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


def _parse_relax(path: str) -> list[dict]:
    """Parse space-delimited relax_monitor.dat (8 columns, skip # header)."""
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            try:
                rows.append({
                    "step":     int(parts[0]),
                    "t":        float(parts[1]),
                    "tot_ice":  float(parts[2]),
                    "tot_sed":  float(parts[3]),
                    "ice_air":  float(parts[4]),
                    "sed_air":  float(parts[5]),
                    "ice_sed":  float(parts[6]),
                    "tot_trip": float(parts[7]),
                })
            except ValueError:
                continue
    return rows


def _parse_max_snes(log: str) -> int:
    """Max Newton iterations from 'Nonlinear solve … iterations N' lines only."""
    iters = []
    for line in log.splitlines():
        m = re.search(r"Nonlinear solve.*?iterations\s+(\d+)", line, re.I)
        if m:
            iters.append(int(m.group(1)))
    return max(iters) if iters else -1


def _phase_ok(log: str, rc: int) -> bool:
    if rc not in (0, None):
        return False
    return re.search(r"ABORT|phase.*out.of.bounds|blow.?up", log, re.I) is None


# ---------------------------------------------------------------------------
# Run helper
# ---------------------------------------------------------------------------
def _run(binary: str, opts_file: str, extra: list[str],
         run_dir: str, timeout: int) -> tuple[str, int]:
    cmd = (
        ["mpirun", "-n", "4", binary, "-options_file", os.path.abspath(opts_file)]
        + extra
        + ["-output_path", run_dir]
    )
    env = os.environ.copy()
    env["folder"] = run_dir
    os.makedirs(run_dir, exist_ok=True)
    try:
        res = subprocess.run(
            cmd, capture_output=True, text=True, env=env, timeout=timeout,
            cwd=os.path.dirname(os.path.abspath(binary)) or ".",
        )
        combined = res.stdout + res.stderr
        # outp.txt mirrors what run_permafrost.sh does via `| tee`
        with open(os.path.join(run_dir, "outp.txt"), "w") as f:
            f.write(res.stdout)
        with open(os.path.join(run_dir, "run.log"), "w") as f:
            f.write(combined)
        return combined, res.returncode
    except subprocess.TimeoutExpired:
        print(f"    [TIMEOUT after {timeout}s]")
        return "", -2


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------
def _dt_metrics(rows: list[dict]) -> tuple[float, float]:
    """(dt_final, dt_growth_frac)."""
    if not rows:
        return float("nan"), float("nan")
    dts = [r["dt"] for r in rows]
    dt_final = dts[-1]
    if len(dts) < 2:
        return dt_final, float("nan")
    n_grow = sum(1 for a, b in zip(dts, dts[1:]) if b > a)
    return dt_final, n_grow / (len(dts) - 1)


def _ice_sed_metric(relax_rows: list[dict], tsf: float) -> float:
    """
    % decrease in ice-sed interface from the freeze point to the end of the run.
    Positive = interface shrank (air invaded — bad).
    Near zero = interface stable (good).
    """
    if not relax_rows:
        return float("nan")
    freeze_idx = next(
        (i for i, r in enumerate(relax_rows) if r["t"] >= tsf), 0
    )
    ref = relax_rows[freeze_idx]["ice_sed"]
    end = relax_rows[-1]["ice_sed"]
    if ref < 1e-30:
        return float("nan")
    return (ref - end) / ref * 100.0  # positive = shrinkage (bad)


def _air_growth_metric(rows: list[dict]) -> float:
    """% growth of total air from step 0 to end of run."""
    if len(rows) < 2:
        return float("nan")
    air0 = rows[0]["tot_air"]
    airN = rows[-1]["tot_air"]
    if air0 < 1e-30:
        return float("nan")
    return (airN - air0) / air0 * 100.0


# ---------------------------------------------------------------------------
# Single-run evaluation
# ---------------------------------------------------------------------------
def evaluate(binary, opts_file, tsf, kp, out_dir, t_final, n_relax, timeout):
    label   = f"tsf={tsf:.1e}_kp={kp:.1e}"
    run_dir = os.path.join(out_dir, label)
    print(f"\n  ── {label} ──")

    extra = [
        "-t_sed_freeze",      str(tsf),
        "-k_pen",             str(kp),
        "-t_final",           str(t_final),
        "-n_relax",           str(n_relax),
        "-Permafrost_output",  "0",  # skip sol*.dat to save disk
        "-Permafrost_monitor", "1",
    ]

    log, rc = _run(binary, opts_file, extra, run_dir, timeout)
    if rc == -2:
        return _null_row(tsf, kp, "TIMEOUT")

    outp_rows  = _parse_outp(os.path.join(run_dir, "outp.txt"))
    relax_rows = _parse_relax(os.path.join(run_dir, "relax_monitor.dat"))

    if not outp_rows:
        return _null_row(tsf, kp, "NO_DATA")

    dt_final, dt_grow_frac = _dt_metrics(outp_rows)
    ised_chg_pct           = _ice_sed_metric(relax_rows, tsf)
    air_growth_pct         = _air_growth_metric(outp_rows)
    max_snes               = _parse_max_snes(log)
    ph_ok                  = _phase_ok(log, rc)

    p_dt    = not np.isnan(dt_grow_frac)  and dt_grow_frac   >= DT_GROW_THRESH
    p_ised  = np.isnan(ised_chg_pct)      or  ised_chg_pct   <  ISED_THRESH
    p_air   = np.isnan(air_growth_pct)    or  air_growth_pct  <  AIR_THRESH
    p_snes  = max_snes >= 0               and max_snes        <= SNES_THRESH
    p_phase = ph_ok
    passed  = p_dt and p_ised and p_air and p_snes and p_phase

    row = {
        "t_sed_freeze":    tsf,
        "k_pen":           kp,
        "rc":              rc,
        "n_steps":         len(outp_rows),
        "dt_final":        dt_final,
        "dt_growth_frac":  dt_grow_frac,
        "ised_chg_pct":    ised_chg_pct,
        "air_growth_pct":  air_growth_pct,
        "max_snes_iters":  max_snes,
        "phase_ok":        ph_ok,
        "pass_dt":         p_dt,
        "pass_ised":       p_ised,
        "pass_air":        p_air,
        "pass_snes":       p_snes,
        "pass_phase":      p_phase,
        "PASS":            passed,
        "status":          "OK",
    }
    _print_row(row)
    return row


def _null_row(tsf, kp, status):
    row = {k: float("nan") for k in [
        "dt_final", "dt_growth_frac", "ised_chg_pct",
        "air_growth_pct", "max_snes_iters",
    ]}
    row.update({
        "t_sed_freeze": tsf, "k_pen": kp, "rc": -1, "n_steps": 0,
        "phase_ok": False,
        "pass_dt": False, "pass_ised": False, "pass_air": False,
        "pass_snes": False, "pass_phase": False, "PASS": False,
        "status": status,
    })
    _print_row(row)
    return row


def _fmt(v, fmt=".2e"):
    if isinstance(v, (int, np.integer)):
        return str(v)
    if isinstance(v, float) and not np.isnan(v):
        return f"{v:{fmt}}"
    return "N/A"


def _print_row(r):
    verdict = "PASS" if r["PASS"] else "FAIL"
    print(f"    dt_final={_fmt(r['dt_final'])}  grow={_fmt(r['dt_growth_frac'],'.2f')}  "
          f"ised_chg={_fmt(r['ised_chg_pct'],'.2f')}%  air_grow={_fmt(r['air_growth_pct'],'.2f')}%  "
          f"snes={r['max_snes_iters']}  phase={'OK' if r['phase_ok'] else 'FAIL'}  → {verdict}")


# ---------------------------------------------------------------------------
# ASCII table
# ---------------------------------------------------------------------------
def _print_table(results: list[dict], tsf_vals, kp_vals):
    print("\n" + "═" * 115)
    print("  SWEEP RESULTS — t_sed_freeze × k_pen")
    print("═" * 115)
    print(f"  {'tsf':>9} {'k_pen':>10} | "
          f"{'dt_final':>10} {'grow_f':>7} | "
          f"{'ised_chg%':>10} {'air_grow%':>10} | "
          f"{'snes':>5} {'phase':>6} | "
          f"{'PASS':>5}")
    print("  " + "─" * 111)
    for r in results:
        snes = int(r["max_snes_iters"]) if not np.isnan(float(r["max_snes_iters"])) else "N/A"
        print(
            f"  {r['t_sed_freeze']:>9.2e} {r['k_pen']:>10.2e} | "
            f"{_fmt(r['dt_final']):>10} {_fmt(r['dt_growth_frac'],'.2f'):>7} | "
            f"{_fmt(r['ised_chg_pct'],'.2f'):>10} {_fmt(r['air_growth_pct'],'.2f'):>10} | "
            f"{str(snes):>5} {'OK' if r['phase_ok'] else 'FAIL':>6} | "
            f"{'✓' if r['PASS'] else '✗':>5}"
        )
    print("═" * 115)

    passed = [r for r in results if r["PASS"]]
    if passed:
        best_tsf = min(r["t_sed_freeze"] for r in passed)
        print(f"\n  Minimum passing t_sed_freeze : {best_tsf:.2e} s")
        print(f"  ({len(passed)}/{len(results)} parameter combinations pass all criteria)")
    else:
        print("\n  No combinations passed — consider relaxing thresholds or "
              "increasing t_final to allow more post-freeze data.")
    print()


# ---------------------------------------------------------------------------
# 2-D heatmaps
# ---------------------------------------------------------------------------
def _heatmap(results: list[dict], tsf_vals: list, kp_vals: list, out_dir: str):
    nt = len(tsf_vals)
    nk = len(kp_vals)

    def _grid(key, default=np.nan):
        mat = np.full((nt, nk), default, dtype=float)
        for r in results:
            ti = tsf_vals.index(r["t_sed_freeze"])
            ki = kp_vals.index(r["k_pen"])
            v = r[key]
            mat[ti, ki] = float(v) if not (isinstance(v, float) and np.isnan(v)) else default
        return mat

    dt_grid   = _grid("dt_final")
    ised_grid = _grid("ised_chg_pct")
    air_grid  = _grid("air_growth_pct")
    snes_grid = _grid("max_snes_iters", default=-1.0)
    pass_grid = _grid("PASS", default=0.0)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle("Sediment freeze sweep: t_sed_freeze × k_pen", fontsize=13)

    kp_labels  = [f"{v:.0e}" for v in kp_vals]
    tsf_labels = [f"{v:.0e}" for v in tsf_vals]

    def _draw(ax, data, title, cmap, norm=None, fmt=".2e", bad="lightgray"):
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad(color=bad)
        masked = np.ma.masked_invalid(data)
        im = ax.imshow(masked, aspect="auto", cmap=cmap_obj, norm=norm, origin="upper")
        ax.set_xticks(range(nk)); ax.set_xticklabels(kp_labels, fontsize=8, rotation=30)
        ax.set_yticks(range(nt)); ax.set_yticklabels(tsf_labels, fontsize=8)
        ax.set_xlabel("k_pen", fontsize=9)
        ax.set_ylabel("t_sed_freeze [s]", fontsize=9)
        ax.set_title(title, fontsize=10)
        plt.colorbar(im, ax=ax, shrink=0.75)
        for ti in range(nt):
            for ki in range(nk):
                v = data[ti, ki]
                if np.isfinite(v):
                    ax.text(ki, ti, f"{v:{fmt}}", ha="center", va="center",
                            fontsize=7)

    # dt_final — log scale, high = good
    dt_pos = np.where(dt_grid > 0, dt_grid, np.nan)
    vdt = dt_pos[np.isfinite(dt_pos)]
    norm_dt = mcolors.LogNorm(vmin=vdt.min(), vmax=vdt.max()) if vdt.size > 1 else None
    _draw(axes[0, 0], dt_pos,
          "Final dt [s]  (↑ larger = better)", "Blues", norm=norm_dt)

    # ice-sed interface change — low = good
    _draw(axes[0, 1], np.clip(ised_grid, -5, 20),
          "Ice-sed interface Δ [%]  (↓ smaller = better)", "Reds", fmt=".1f")

    # air growth — low = good
    _draw(axes[0, 2], np.clip(air_grid, -1, 20),
          "Total air growth [%]  (↓ smaller = better)", "Oranges", fmt=".1f")

    # max SNES iters
    snes_pos = np.where(snes_grid >= 0, snes_grid, np.nan)
    _draw(axes[1, 0], snes_pos,
          "Max SNES Newton iters  (↓ smaller = better)", "Purples", fmt=".0f")

    # dt growth fraction
    _draw(axes[1, 1], _grid("dt_growth_frac"),
          "dt growth fraction  (↑ > 0.5 preferred)", "Greens", fmt=".2f")

    # PASS/FAIL
    _draw(axes[1, 2], pass_grid, "PASS (1) / FAIL (0)",
          "RdYlGn", fmt=".0f", bad="gray")

    plt.tight_layout()
    path = os.path.join(out_dir, "sweep_sed_freeze.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Heatmap saved to: {path}")


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
def _save_csv(results, out_dir):
    path = os.path.join(out_dir, "sweep_sed_freeze.csv")
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
def _parse_args():
    p = argparse.ArgumentParser(
        description="Sweep t_sed_freeze × k_pen for optimal sediment freeze timing.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",      default="./permafrost")
    p.add_argument("--opts",        default="inputs/tests/test_1D_IceSlab.opts",
                   help="Base PETSc options file")
    p.add_argument("--out-dir",     default="/tmp/tune_tsf")
    p.add_argument("--sweep-tsf",   nargs="+", type=float,
                   default=SWEEP_TSF_DEFAULT,
                   help="t_sed_freeze values [s]")
    p.add_argument("--sweep-kpen",  nargs="+", type=float,
                   default=SWEEP_KPEN_DEFAULT,
                   help="k_pen values")
    p.add_argument("--t-final",     type=float, default=2.0,
                   help="Simulation end time per run [s] (default 2.0)")
    p.add_argument("--n-relax",     type=int, default=0,
                   help="n_relax steps before physics (default 0)")
    p.add_argument("--timeout",     type=int, default=300,
                   help="Per-run timeout [s] (default 300)")
    p.add_argument("--dt-grow-thresh", type=float, default=DT_GROW_THRESH)
    p.add_argument("--ised-thresh",    type=float, default=ISED_THRESH,
                   help="Max ice-sed shrinkage %% for PASS (default 5.0)")
    p.add_argument("--air-thresh",     type=float, default=AIR_THRESH,
                   help="Max total air growth %% for PASS (default 10.0)")
    p.add_argument("--snes-thresh",    type=int,   default=SNES_THRESH)
    return p.parse_args()


def main():
    args = _parse_args()

    global DT_GROW_THRESH, ISED_THRESH, AIR_THRESH, SNES_THRESH
    DT_GROW_THRESH = args.dt_grow_thresh
    ISED_THRESH    = args.ised_thresh
    AIR_THRESH     = args.air_thresh
    SNES_THRESH    = args.snes_thresh

    os.makedirs(args.out_dir, exist_ok=True)

    tsf_vals = sorted(set(args.sweep_tsf))
    kp_vals  = sorted(set(args.sweep_kpen))
    total    = len(tsf_vals) * len(kp_vals)

    print(f"\n{'='*65}")
    print(f"  tune_sed_freeze.py")
    print(f"  t_sed_freeze: {tsf_vals}")
    print(f"  k_pen:        {kp_vals}")
    print(f"  t_final/run:  {args.t_final} s   n_relax: {args.n_relax}")
    print(f"  Total runs:   {total}")
    print(f"  Output:       {args.out_dir}")
    print(f"  Pass criteria:")
    print(f"    dt_growth_frac >= {DT_GROW_THRESH}")
    print(f"    ice-sed change <  {ISED_THRESH}%")
    print(f"    air growth     <  {AIR_THRESH}%")
    print(f"    max SNES iters <= {SNES_THRESH}")
    print(f"{'='*65}")

    results = []
    for idx, (tsf, kp) in enumerate(product(tsf_vals, kp_vals), 1):
        print(f"\n[{idx}/{total}]", end="")
        row = evaluate(
            binary=args.binary, opts_file=args.opts,
            tsf=tsf, kp=kp,
            out_dir=args.out_dir,
            t_final=args.t_final,
            n_relax=args.n_relax,
            timeout=args.timeout,
        )
        results.append(row)

    _print_table(results, tsf_vals, kp_vals)
    _save_csv(results, args.out_dir)
    _heatmap(results, tsf_vals, kp_vals, args.out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
