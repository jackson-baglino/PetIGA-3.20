#!/usr/bin/env python3
"""
tune_sed_penalty.py — Sweep t_sed_freeze and k_sed_pen with difvap_pen fixed at 0.

Restores the original vapor equation (no diffusion penalty) and isolates the
effect of the sediment penalty parameters on spurious air generation at the
ice-sediment boundary.

IC: test_1D_IceSlab.opts — ice slab [0.35 Lx, 0.65 Lx] adjacent to a
sediment slab [0.65 Lx, 0.85 Lx].  flag_sed_mode=1 is enforced: sediment
evolves freely for t < t_sed_freeze, then switches to the penalty.

Metrics (from monitoring output + relax_monitor.dat):
  - Δtot_air    = tot_air(final) − tot_air(at t_sed_freeze)
                  ≈ 0 when no spurious air appears after the freeze.
  - ice_sed_drop = relative drop in ∫(φ_ice²·φ_sed²) after freeze [%]
                  ≈ 0 when the ice-sed interface is stable.

Pass criteria:
  1. Clean exit (rc == 0).
  2. |Δtot_air|    < --pass-air-drift   (default 5e-6).
  3. ice_sed_drop  < --pass-icesd-drop  (default 20 %).

Usage
-----
  cd <project_root>
  python test/tune_sed_penalty.py

  python test/tune_sed_penalty.py \\
      --t-freeze 10 50 100 300 900 \\
      --k-prefactor 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 \\
      --t-final 3600 --timeout 600
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

ROOT = Path(__file__).parent.parent.resolve()

# ── defaults ──────────────────────────────────────────────────────────────────
BINARY_DEFAULT  = str(ROOT / "permafrost")
OPTS_DEFAULT    = str(ROOT / "inputs" / "tests" / "test_1D_IceSlab.opts")
OUTDIR_DEFAULT  = str(ROOT / "test" / "tune_sed_penalty")
TIMEOUT_DEFAULT = 600    # s wall-clock per run
T_FINAL_DEFAULT = 3600.0 # s simulated time

EPS_DEFAULT = 9.3295e-7  # m — interface width (k_sed_pen = prefactor / eps²)

T_SED_FREEZE_DEFAULT    = [10.0, 50.0, 100.0, 300.0, 900.0]
K_SED_PREFACTOR_DEFAULT = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]

PASS_AIR_DRIFT_DEFAULT  = 5e-6   # absolute |Δtot_air| [m in 1-D]
PASS_ICESD_DROP_DEFAULT = 20.0   # % drop in ice_sed_interf after freeze


# ── monitor / bounds parsers ──────────────────────────────────────────────────
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"

_MROW = re.compile(
    r"^\s+(\d+)\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|\s*(" + _NUM + r")\s*\|"
    r"\s*(" + _NUM + r")\s*$",
    re.MULTILINE,
)

_BOUNDS_LINE = re.compile(
    r"BOUNDS:\s+phi_ice\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
    r"\s+phi_sed\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
    r"\s+phi_air\s+\[(" + _NUM + r"),\s*(" + _NUM + r")\]"
)
PHASE_LO, PHASE_HI = -0.25, 1.25


def _parse_monitor_rows(stdout: str) -> list[dict]:
    rows = []
    for m in _MROW.finditer(stdout):
        rows.append({
            "step":    int(m.group(1)),
            "time":    float(m.group(2)),
            "tot_ice": float(m.group(4)),
            "tot_air": float(m.group(5)),
            "tot_sed": float(m.group(6)),
        })
    return rows


def _parse_bounds_violation(stdout: str) -> bool:
    for m in _BOUNDS_LINE.finditer(stdout):
        if any(float(m.group(i)) < PHASE_LO or float(m.group(i)) > PHASE_HI
               for i in range(1, 7)):
            return True
    return False


def _parse_relax_monitor(run_dir: str) -> list[dict]:
    rows = []
    try:
        with open(os.path.join(run_dir, "relax_monitor.dat")) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or not line:
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    rows.append({"t": float(parts[1]),
                                 "ice_sed_interf": float(parts[6])})
    except OSError:
        pass
    return rows


def _compute_metrics(rows: list[dict], relax_rows: list[dict],
                     t_sed_freeze: float) -> dict:
    if not rows:
        return {"delta_tot_air": None, "ice_sed_drop_pct": None}

    freeze_idx = next(
        (i for i, r in enumerate(rows) if r["time"] >= t_sed_freeze), None
    )
    if freeze_idx is None:
        return {"delta_tot_air": None, "ice_sed_drop_pct": None,
                "note": "ended before t_sed_freeze"}

    delta_tot_air = rows[-1]["tot_air"] - rows[freeze_idx]["tot_air"]

    ice_sed_drop_pct = None
    if relax_rows:
        fri = next((i for i, r in enumerate(relax_rows) if r["t"] >= t_sed_freeze), None)
        if fri is not None:
            isi0 = relax_rows[fri]["ice_sed_interf"]
            isif = relax_rows[-1]["ice_sed_interf"]
            if abs(isi0) > 1e-20:
                ice_sed_drop_pct = (isi0 - isif) / abs(isi0) * 100.0

    return {"delta_tot_air": delta_tot_air, "ice_sed_drop_pct": ice_sed_drop_pct}


# ── build / run helpers ───────────────────────────────────────────────────────
def _build(root: Path) -> None:
    for cmd in (["make", "clean"], ["make"]):
        label = " ".join(cmd)
        print(f"\n  [build] {label}")
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, cwd=str(root),
        )
        for line in proc.stdout:
            sys.stdout.write(f"    {line}"); sys.stdout.flush()
        proc.wait()
        if proc.returncode != 0:
            sys.exit(f"\n  [build] ERROR: '{label}' failed. Aborting.")
    print("  [build] Build successful.")


def _run(binary: str, opts_file: str, extra_flags: list[str],
         run_dir: str, timeout: int) -> tuple[str, int]:
    cmd = ["mpiexec", "-np", "4", binary,
           "-options_file", os.path.abspath(opts_file)] + extra_flags
    env = os.environ.copy()
    env["folder"] = run_dir
    os.makedirs(run_dir, exist_ok=True)
    collected: list[str] = []

    try:
        with open(os.path.join(run_dir, "outp.txt"), "w") as f:
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
                proc.kill(); t.join(timeout=5)
                print(f"\n    [TIMEOUT after {timeout}s]")
                return "".join(collected), -2
            t.join()
        return "".join(collected), proc.returncode
    except Exception as e:
        print(f"    [ERROR: {e}]")
        return "", -3


# ── sweep engine ──────────────────────────────────────────────────────────────
def run_sweep(binary: str, opts_file: str,
              t_freeze_vals: list[float],
              k_prefactor_vals: list[float],
              out_dir: str, timeout: int, t_final: float,
              eps: float,
              pass_air_drift: float,
              pass_icesd_drop: float) -> list[dict]:
    results = []
    total  = len(t_freeze_vals) * len(k_prefactor_vals)
    run_no = 0

    for t_freeze in t_freeze_vals:
        for k_pref in k_prefactor_vals:
            run_no += 1
            k_sed_pen = k_pref / (eps * eps)
            label = f"tsf{t_freeze:.0f}_ksp{k_pref:.1e}".replace("+", "")

            print(f"\n{'='*68}")
            print(f"  Run {run_no}/{total}: t_sed_freeze={t_freeze:.0f} s  "
                  f"k_pref={k_pref:.2e}  (k_sed_pen={k_sed_pen:.3e})")
            print(f"{'='*68}")

            extra_flags = [
                "-t_sed_freeze", str(t_freeze),
                "-k_sed_pen",    str(k_sed_pen),
                "-difvap_pen",   "0.0",
                "-t_final",      str(t_final),
                "-flag_sed_mode", "1",
            ]
            run_dir = os.path.join(out_dir, label)
            stdout, rc = _run(binary, opts_file, extra_flags, run_dir, timeout)

            rows       = _parse_monitor_rows(stdout)
            relax_rows = _parse_relax_monitor(run_dir)
            bounds_vio = _parse_bounds_violation(stdout)
            metrics    = _compute_metrics(rows, relax_rows, t_freeze)

            delta   = metrics.get("delta_tot_air")
            isdrop  = metrics.get("ice_sed_drop_pct")

            pass_exit   = (rc == 0)
            pass_bounds = not bounds_vio
            pass_air    = delta is not None and abs(delta) < pass_air_drift
            pass_icesd  = isdrop is None or isdrop < pass_icesd_drop
            passed      = pass_exit and pass_bounds and pass_air and pass_icesd

            results.append({
                "t_sed_freeze":      t_freeze,
                "k_sed_prefactor":   k_pref,
                "k_sed_pen":         k_sed_pen,
                "rc":                rc,
                "delta_tot_air":     delta,
                "ice_sed_drop_pct":  isdrop,
                "pass_exit":         pass_exit,
                "pass_bounds":       pass_bounds,
                "pass_air":          pass_air,
                "pass_icesd":        pass_icesd,
                "PASS":              passed,
            })

            dstr = f"{delta:.4e}" if delta is not None else "N/A"
            istr = f"{isdrop:.1f}%" if isdrop is not None else "N/A"
            print(f"  rc={rc}  Δtot_air={dstr}  ice_sed_drop={istr}  "
                  f"→ {'PASS' if passed else 'FAIL'}")

    return results


# ── plotting ──────────────────────────────────────────────────────────────────
def _plot(results: list[dict], out_dir: str,
          t_freeze_vals: list[float],
          k_prefactor_vals: list[float],
          pass_air_drift: float,
          pass_icesd_drop: float) -> None:

    nrows = len(t_freeze_vals)
    ncols = len(k_prefactor_vals)
    ti = {v: i for i, v in enumerate(t_freeze_vals)}
    ki = {v: i for i, v in enumerate(k_prefactor_vals)}

    air_grid  = np.full((nrows, ncols), np.nan)
    drop_grid = np.full((nrows, ncols), np.nan)
    pass_grid = np.zeros((nrows, ncols), dtype=bool)

    for r in results:
        i = ti.get(r["t_sed_freeze"])
        j = ki.get(r["k_sed_prefactor"])
        if i is None or j is None:
            continue
        if r["delta_tot_air"] is not None:
            air_grid[i, j]  = abs(r["delta_tot_air"])
        if r["ice_sed_drop_pct"] is not None:
            drop_grid[i, j] = r["ice_sed_drop_pct"]
        pass_grid[i, j] = r["PASS"]

    klab = [f"{k:.1e}" for k in k_prefactor_vals]
    tlab = [f"{t:.0f}" for t in t_freeze_vals]

    def _heatmap(grid, title, cbar_label, vmax, fname):
        fig, ax = plt.subplots(figsize=(max(6, ncols * 1.8), max(4, nrows * 1.2)))
        im = ax.imshow(grid, aspect="auto", cmap="RdYlGn_r",
                       vmin=0, vmax=vmax, origin="upper")
        plt.colorbar(im, ax=ax, label=cbar_label)
        ax.set_xticks(range(ncols)); ax.set_xticklabels(klab, rotation=30, ha="right")
        ax.set_yticks(range(nrows)); ax.set_yticklabels(tlab)
        ax.set_xlabel("k_sed_pen prefactor  (k_sed_pen = prefactor / ε²)", fontsize=10)
        ax.set_ylabel("t_sed_freeze  [s]", fontsize=10)
        ax.set_title(title, fontsize=11)
        for i in range(nrows):
            for j in range(ncols):
                val = grid[i, j]
                marker = "✓" if pass_grid[i, j] else "✗"
                txt = f"{val:.2e}\n{marker}" if not np.isnan(val) else "N/A"
                ax.text(j, i, txt, ha="center", va="center", fontsize=9)
        plt.tight_layout()
        fig.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Plot: {fname}")

    _heatmap(air_grid,
             "|Δtot_air| after freeze  (difvap_pen = 0)",
             "|Δtot_air|  [m in 1-D]",
             pass_air_drift * 4,
             os.path.join(out_dir, "heatmap_air_drift.png"))

    _heatmap(drop_grid,
             "ice_sed_interf drop % after freeze  (difvap_pen = 0)",
             "ice_sed_interf drop  [%]",
             pass_icesd_drop * 2,
             os.path.join(out_dir, "heatmap_icesd_drop.png"))

    # Line plot: Δtot_air vs k_pref for each t_sed_freeze
    fig, ax = plt.subplots(figsize=(8, 5))
    k_arr = np.array(k_prefactor_vals)
    for i, t_freeze in enumerate(t_freeze_vals):
        vals = air_grid[i, :]
        style = "-o" if not np.all(np.isnan(vals)) else "--"
        ax.semilogx(k_arr, vals, style, label=f"t_freeze={t_freeze:.0f} s", lw=1.8)
    ax.axhline(pass_air_drift, color="gray", ls="--", lw=1.2,
               label=f"threshold = {pass_air_drift:.1e}")
    ax.set_xlabel("k_sed_pen prefactor", fontsize=11)
    ax.set_ylabel("|Δtot_air|  [m in 1-D]", fontsize=11)
    ax.set_title("Spurious air vs. sediment penalty  (difvap_pen = 0)", fontsize=11)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fname = os.path.join(out_dir, "lineplot_air_drift.png")
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {fname}")


# ── CSV + summary ─────────────────────────────────────────────────────────────
def _save_csv(results: list[dict], out_dir: str) -> None:
    if not results:
        return
    path = os.path.join(out_dir, "sweep_sed_penalty.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader(); w.writerows(results)
    print(f"  CSV: {path}")


def _print_summary(results: list[dict],
                   pass_air_drift: float,
                   pass_icesd_drop: float) -> None:
    print(f"\n{'─'*82}")
    print("  SUMMARY — sediment penalty sweep  (difvap_pen = 0)")
    print(f"  {'t_freeze':>9}  {'k_pref':>9}  {'k_sed_pen':>12}  "
          f"{'Δtot_air':>11}  {'icesd_drop':>11}  {'PASS':>5}")
    print(f"  {'─'*72}")
    for r in results:
        dstr = f"{r['delta_tot_air']:.3e}"    if r["delta_tot_air"]    is not None else "N/A"
        istr = f"{r['ice_sed_drop_pct']:.1f}%" if r["ice_sed_drop_pct"] is not None else "N/A"
        print(f"  {r['t_sed_freeze']:>9.0f}  {r['k_sed_prefactor']:>9.2e}  "
              f"{r['k_sed_pen']:>12.3e}  {dstr:>11}  {istr:>11}  "
              f"{'✓' if r['PASS'] else '✗':>5}")

    passed = [r for r in results if r["PASS"]]
    if passed:
        best = min(passed,
                   key=lambda r: (r["t_sed_freeze"], abs(r["delta_tot_air"] or 1e30)))
        print(f"\n  Best (smallest t_freeze that passes):")
        print(f"    t_sed_freeze = {best['t_sed_freeze']:.0f} s")
        print(f"    k_sed_pen    = {best['k_sed_pen']:.3e}  "
              f"(prefactor {best['k_sed_prefactor']:.2e})")
        print(f"    Δtot_air     = {best['delta_tot_air']:.3e}")
        print(f"  {len(passed)}/{len(results)} passed "
              f"(|Δtot_air|<{pass_air_drift:.1e}, "
              f"ice_sed_drop<{pass_icesd_drop:.0f}%).")
    else:
        print("\n  No combinations passed. "
              "Consider relaxing --pass-air-drift or --pass-icesd-drop.")
    print(f"{'─'*82}")


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Sweep t_sed_freeze and k_sed_pen (difvap_pen=0).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",          default=BINARY_DEFAULT)
    p.add_argument("--opts",            default=OPTS_DEFAULT,
                   help="Opts file (default: test_1D_IceSlab.opts)")
    p.add_argument("--out-dir",         default=OUTDIR_DEFAULT)
    p.add_argument("--t-freeze",        nargs="+", type=float,
                   default=T_SED_FREEZE_DEFAULT, metavar="T",
                   help="t_sed_freeze values [s]")
    p.add_argument("--k-prefactor",     nargs="+", type=float,
                   default=K_SED_PREFACTOR_DEFAULT, metavar="K",
                   help="k_sed_pen prefactor values (k_sed_pen = val / eps²)")
    p.add_argument("--t-final",         type=float, default=T_FINAL_DEFAULT)
    p.add_argument("--timeout",         type=int,   default=TIMEOUT_DEFAULT)
    p.add_argument("--eps",             type=float, default=EPS_DEFAULT)
    p.add_argument("--pass-air-drift",  type=float, default=PASS_AIR_DRIFT_DEFAULT)
    p.add_argument("--pass-icesd-drop", type=float, default=PASS_ICESD_DROP_DEFAULT)
    p.add_argument("--skip-build",      action="store_true",
                   help="Skip 'make clean && make'")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    total = len(args.t_freeze) * len(args.k_prefactor)
    print("=" * 68)
    print("  SEDIMENT PENALTY SWEEP  (difvap_pen = 0)")
    print(f"  binary:          {args.binary}")
    print(f"  opts:            {args.opts}")
    print(f"  out_dir:         {args.out_dir}")
    print(f"  t_sed_freeze:    {args.t_freeze}")
    print(f"  k_prefactor:     {args.k_prefactor}")
    print(f"  t_final:         {args.t_final:.4g} s")
    print(f"  timeout:         {args.timeout} s/run")
    print(f"  total runs:      {total}")
    print(f"  pass_air_drift:  {args.pass_air_drift:.1e}")
    print(f"  pass_icesd_drop: {args.pass_icesd_drop:.0f} %")
    print(f"  build:           {'skipped' if args.skip_build else 'make clean && make'}")
    print("=" * 68)

    if not args.skip_build:
        _build(ROOT)

    results = run_sweep(
        args.binary, args.opts,
        args.t_freeze, args.k_prefactor,
        args.out_dir, args.timeout, args.t_final,
        eps=args.eps,
        pass_air_drift=args.pass_air_drift,
        pass_icesd_drop=args.pass_icesd_drop,
    )

    _print_summary(results, args.pass_air_drift, args.pass_icesd_drop)
    _save_csv(results, args.out_dir)
    _plot(results, args.out_dir,
          args.t_freeze, args.k_prefactor,
          args.pass_air_drift, args.pass_icesd_drop)

    print(f"\n  All output in: {args.out_dir}")


if __name__ == "__main__":
    main()
