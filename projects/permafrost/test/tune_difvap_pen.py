#!/usr/bin/env python3
"""
tune_difvap_pen.py — 3D sweep of difvap_pen × k_pen × t_sed_freeze.

For each grid point the script runs one simulation using test_2D_IceSlab.opts
(or a user-supplied opts file).  Only the three swept parameters are overridden
on the command line; t_final, physics flags, and everything else come from the
opts file unchanged.

Goal: find parameter combinations that suppress spurious air generation at the
ice–sediment interface after the sediment freeze event.

Primary metric
  delta_tot_air = tot_air(final) − tot_air(at t_sed_freeze)
  Positive values indicate air grew after the freeze = spurious air (bad).

Secondary metrics
  ice_sed_drop_pct — % drop in ice-sediment interface area after freeze
  max_snes_iters   — solver health
  phase bounds      — any phase field outside [PHASE_LO, PHASE_HI]

PASS = all four criteria simultaneously satisfied.

Usage
-----
  cd <project_root>
  python test/tune_difvap_pen.py [--skip-build]

  # 3D sweep (difvap_pen × k_pen × t_sed_freeze):
  python test/tune_difvap_pen.py --skip-build \\
      --sweep-difvap  1e-7 1e-6 1e-5 \\
      --sweep-kpen    1e7 1e8 1e9 \\
      --sweep-tsf     0.0 1e-4 1e-3 1e-2 \\
      --timeout 600

  # 2D sweep (t_sed_freeze fixed at opts value):
  python test/tune_difvap_pen.py \\
      --sweep-difvap 1e-6 1e-5 1e-4 \\
      --sweep-kpen   1e6  1e7  1e8  \\
      --timeout 600

  # Single-point smoke test:
  python test/tune_difvap_pen.py --skip-build \\
      --sweep-difvap 1e-5 --sweep-kpen 1e7 --timeout 120

  # Use a different opts file:
  python test/tune_difvap_pen.py \\
      --opts inputs/tests/test_1D_IceSlab.opts

Note on simulation length
  This script passes NO -t_final override.  For the sweep to run in reasonable
  wall-clock time AND capture spurious air after the freeze, set t_final in the
  opts file to approximately max(t_sed_freeze) + 300*delt_t.  The --timeout
  flag provides a hard wall-clock cutoff per run.
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

# ── project root (one level up from test/) ───────────────────────────────────
ROOT = Path(__file__).parent.parent.resolve()

# ── defaults ──────────────────────────────────────────────────────────────────
BINARY_DEFAULT = str(ROOT / "permafrost")
OPTS_DEFAULT   = str(ROOT / "inputs" / "tests" / "test_2D_IceSlab.opts")
OUTDIR_DEFAULT = str(ROOT / "test" / "tune_penalty2d")
TIMEOUT_DEFAULT = 600  # s wall-clock per run

# Default sweep grids
DIFVAP_SWEEP_DEFAULT = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
KPEN_SWEEP_DEFAULT   = [1e5,  1e6,  1e7,  1e8,  1e9]
TSF_SWEEP_DEFAULT    = None  # None → read t_sed_freeze from opts file (no tsf sweep)

# Pass/fail thresholds
PASS_AIR_DRIFT_DEFAULT   = 5e-6   # absolute |Δtot_air| [m in 1-D, m² in 2-D]
PASS_ICESD_DROP_DEFAULT  = 20.0   # % drop in ice-sed interface area after freeze
PASS_SNES_ITERS          = 450    # max Newton iterations per step

# Phase bounds for out-of-bounds detection
PHASE_LO = -0.25
PHASE_HI  =  1.25


# ── monitor line parser ───────────────────────────────────────────────────────
# Monitor table columns (from monitoring.c):
#   STEP | TIME | DT | TOT_ICE | TOT_AIR | TOT_SED | TEMP | TOT_RHOV | I-A INTERF | TRIPL_JUNC
#     1      2    3      4          5         6        7        8            9            10
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


def _parse_monitor_rows(stdout: str) -> list[dict]:
    """Parse every 10-column monitor table row from stdout."""
    rows = []
    for m in _MROW.finditer(stdout):
        rows.append({
            "step":     int(m.group(1)),
            "time":     float(m.group(2)),
            "tot_ice":  float(m.group(4)),
            "tot_air":  float(m.group(5)),
            "tot_sed":  float(m.group(6)),
            "tot_rhov": float(m.group(8)),
        })
    return rows


def _parse_bounds_violation(stdout: str) -> bool:
    """Return True if any phase field is reported outside [PHASE_LO, PHASE_HI]."""
    for m in _BOUNDS_LINE.finditer(stdout):
        vals = [float(m.group(i)) for i in range(1, 7)]
        if any(v < PHASE_LO or v > PHASE_HI for v in vals):
            return True
    return False


def _parse_snes_iters(stdout: str) -> int:
    """Return the maximum Newton iteration count seen in stdout."""
    iters = [int(m.group(1))
             for line in stdout.splitlines()
             for m in [re.search(r"converged.*?iterations\s+(\d+)", line, re.I)]
             if m]
    return max(iters) if iters else -1


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
                        "step":           int(parts[0]),
                        "t":              float(parts[1]),
                        "ice_sed_interf": float(parts[6]),
                    })
    except OSError:
        pass
    return rows


def _compute_metrics(rows: list[dict], relax_rows: list[dict],
                     t_sed_freeze: float) -> dict:
    """
    Compute delta_tot_air and ice_sed_interf stability after t_sed_freeze.
    Both metrics are measured relative to the first step where t >= t_sed_freeze.
    """
    if not rows:
        return {"delta_tot_air": None, "ice_sed_drop_pct": None}

    freeze_idx = next(
        (i for i, r in enumerate(rows) if r["time"] >= t_sed_freeze), None
    )
    if freeze_idx is None:
        return {
            "delta_tot_air": None,
            "ice_sed_drop_pct": None,
            "note": "simulation ended before t_sed_freeze",
        }

    tot_air_at_freeze = rows[freeze_idx]["tot_air"]
    tot_air_final     = rows[-1]["tot_air"]
    delta_tot_air     = tot_air_final - tot_air_at_freeze

    ice_sed_drop_pct = None
    if relax_rows:
        freeze_relax_idx = next(
            (i for i, r in enumerate(relax_rows) if r["t"] >= t_sed_freeze), None
        )
        if freeze_relax_idx is not None:
            isi0 = relax_rows[freeze_relax_idx]["ice_sed_interf"]
            isif = relax_rows[-1]["ice_sed_interf"]
            if abs(isi0) > 1e-20:
                ice_sed_drop_pct = (isi0 - isif) / abs(isi0) * 100.0

    return {
        "delta_tot_air":     delta_tot_air,
        "tot_air_at_freeze": tot_air_at_freeze,
        "tot_air_final":     tot_air_final,
        "ice_sed_drop_pct":  ice_sed_drop_pct,
    }


# ── opts file helpers ─────────────────────────────────────────────────────────

def _get_t_sed_freeze(opts_file: str, default: float = 300.0) -> float:
    """Parse -t_sed_freeze from an opts file; returns default if not found."""
    try:
        with open(opts_file) as f:
            for line in f:
                m = re.match(r"^\s*-t_sed_freeze\s+([\d.eE+\-]+)", line)
                if m:
                    return float(m.group(1))
    except OSError:
        pass
    return default


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


# ── post-processing scripts ───────────────────────────────────────────────────
SCRIPT_VTK      = str(ROOT / "scripts"     / "plotpermafrost.py")
SCRIPT_SCALARS  = str(ROOT / "postprocess" / "plot_scalars.py")
SCRIPT_PROFS1D  = str(ROOT / "postprocess" / "plot1D_profiles.py")
SCRIPT_PROFS2D  = str(ROOT / "postprocess" / "plot2D_snapshot.py")


def _postprocess(run_dir: str, dim: int = 2) -> None:
    """Run post-processing on a completed run directory."""
    if not os.path.isfile(os.path.join(run_dir, "igasol.dat")):
        print(f"  [post] No igasol.dat in {run_dir} — skipping.")
        return

    def _sp(cmd: list[str], desc: str) -> None:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               cwd=run_dir, timeout=300)
            if r.returncode != 0:
                print(f"    [post] WARNING: {desc} exited {r.returncode}")
            else:
                print(f"    [post] {desc} done.")
        except subprocess.TimeoutExpired:
            print(f"    [post] WARNING: {desc} timed out.")
        except FileNotFoundError as e:
            print(f"    [post] WARNING: {e}")

    os.makedirs(os.path.join(run_dir, "vtkOut"), exist_ok=True)
    _sp(["python", SCRIPT_VTK, "--vtk-dir", "vtkOut"], "VTK conversion")

    if os.path.isfile(os.path.join(run_dir, "SSA_evo.dat")):
        _sp(["python", SCRIPT_SCALARS,
             "--file", "SSA_evo.dat", "--time-unit", "s", "--save", "scalars.png"],
            "Scalar time-series")

    if dim == 1:
        _sp(["python", SCRIPT_PROFS1D, "--dir", "."],
            "1D profiles (phase + thermal + first/last)")


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

def run_sweep(binary: str, opts_file: str,
              difvap_sweep: list[float],
              kpen_sweep: list[float],
              tsf_sweep: list[float],
              out_dir: str, timeout: int,
              pass_air_drift: float,
              pass_icesd_drop: float,
              skip_postprocess: bool) -> list[dict]:

    dim          = _get_dim(opts_file)
    multi_tsf    = len(tsf_sweep) > 1

    print(f"  opts file:      {opts_file}")
    print(f"  dim={dim}  t_sed_freeze sweep: {tsf_sweep}")

    results = []
    total   = len(difvap_sweep) * len(kpen_sweep) * len(tsf_sweep)
    run_no  = 0

    for d, k, tsf in product(difvap_sweep, kpen_sweep, tsf_sweep):
        run_no += 1
        # Include tsf in directory name only when multiple values are swept
        tsf_tag = f"_t{tsf:.2e}".replace("+", "") if multi_tsf else ""
        label   = f"d{d:.2e}_k{k:.2e}{tsf_tag}".replace("+", "")
        run_dir = os.path.join(out_dir, label)

        print(f"\n{'='*70}")
        print(f"  Run {run_no}/{total}: difvap_pen={d:.2e}  k_pen={k:.2e}  t_sed_freeze={tsf:.2e}")
        print(f"{'='*70}")

        flags = ["-difvap_pen", str(d), "-k_pen", str(k), "-t_sed_freeze", str(tsf)]

        stdout, rc = _run(binary, opts_file, flags, run_dir, timeout)

        rows        = _parse_monitor_rows(stdout)
        relax_rows  = _parse_relax_monitor(run_dir)
        bounds_viol = _parse_bounds_violation(stdout)
        max_snes    = _parse_snes_iters(stdout)
        metrics     = _compute_metrics(rows, relax_rows, tsf)

        delta_tot_air    = metrics.get("delta_tot_air")
        ice_sed_drop_pct = metrics.get("ice_sed_drop_pct")

        pass_exit   = (rc == 0 or rc == -2)  # -2 = timeout; partial output still useful
        pass_bounds = not bounds_viol
        pass_air    = (delta_tot_air is not None
                       and abs(delta_tot_air) < pass_air_drift)
        pass_icesd  = (ice_sed_drop_pct is None      # no relax_monitor.dat — skip
                       or ice_sed_drop_pct < pass_icesd_drop)
        pass_iters  = (max_snes >= 0 and max_snes <= PASS_SNES_ITERS)
        passed      = pass_exit and pass_bounds and pass_air and pass_icesd and pass_iters

        row = {
            "t_sed_freeze":     tsf,
            "difvap_pen":       d,
            "k_pen":            k,
            "rc":               rc,
            "delta_tot_air":    delta_tot_air,
            "ice_sed_drop_pct": ice_sed_drop_pct,
            "max_snes_iters":   max_snes,
            "bounds_ok":        pass_bounds,
            "pass_air":         pass_air,
            "pass_interf":      pass_icesd,
            "pass_iters":       pass_iters,
            "pass_bounds":      pass_bounds,
            "PASS":             passed,
        }
        results.append(row)

        dstr  = (f"{delta_tot_air:.4e}"     if delta_tot_air    is not None else "N/A")
        idstr = (f"{ice_sed_drop_pct:.1f}%" if ice_sed_drop_pct is not None else "N/A")
        print(f"  rc={rc}  Δtot_air={dstr}  ice_sed_drop={idstr}  "
              f"SNES_max={max_snes}  bounds={'OK' if pass_bounds else 'VIOLATION'}")
        print(f"  → {'PASS' if passed else 'FAIL'}")

        if not skip_postprocess:
            _postprocess(run_dir, dim=dim)

    return results


# ── 2D heatmap plotting ───────────────────────────────────────────────────────

def _plot_slice(results_slice: list[dict], difvap_sweep: list[float],
                kpen_sweep: list[float], out_dir: str,
                pass_air_drift: float, pass_icesd_drop: float,
                title_extra: str, fname: str) -> None:
    """
    Four-panel 2D heatmap for one t_sed_freeze slice:
    difvap_pen (y) × k_pen (x).
    Panels: |Δtot_air|, ice_sed_drop_pct, max_snes_iters, PASS mask.
    """
    nd = len(difvap_sweep)
    nk = len(kpen_sweep)

    di = {v: i for i, v in enumerate(difvap_sweep)}
    ki = {v: i for i, v in enumerate(kpen_sweep)}

    air_grid   = np.full((nd, nk), np.nan)
    icesd_grid = np.full((nd, nk), np.nan)
    snes_grid  = np.full((nd, nk), np.nan)
    pass_grid  = np.zeros((nd, nk), dtype=bool)

    for r in results_slice:
        i = di.get(r["difvap_pen"])
        j = ki.get(r["k_pen"])
        if i is None or j is None:
            continue
        if r["delta_tot_air"] is not None:
            air_grid[i, j]   = abs(r["delta_tot_air"])
        if r["ice_sed_drop_pct"] is not None:
            icesd_grid[i, j] = r["ice_sed_drop_pct"]
        if r["max_snes_iters"] >= 0:
            snes_grid[i, j]  = r["max_snes_iters"]
        pass_grid[i, j] = r["PASS"]

    x_labels = [f"{k:.1e}" for k in kpen_sweep]
    y_labels  = [f"{d:.1e}" for d in difvap_sweep]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"difvap_pen × k_pen sweep — spurious air diagnostics{title_extra}",
                 fontsize=13)

    def _hm(ax, data, title, vmin, vmax, cmap, fmt):
        im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                       origin="upper", interpolation="nearest")
        plt.colorbar(im, ax=ax)
        ax.set_xticks(range(nk)); ax.set_xticklabels(x_labels, rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(nd)); ax.set_yticklabels(y_labels, fontsize=8)
        ax.set_xlabel("k_pen", fontsize=10)
        ax.set_ylabel("difvap_pen  [m²/s]", fontsize=10)
        ax.set_title(title, fontsize=11)
        for ii in range(nd):
            for jj in range(nk):
                v = data[ii, jj]
                txt = (fmt % v) if not np.isnan(v) else "N/A"
                ax.text(jj, ii, txt, ha="center", va="center", fontsize=7, color="black")

    _hm(axes[0, 0], air_grid,
        f"|Δtot_air|  (threshold = {pass_air_drift:.1e})",
        0, pass_air_drift * 5, "RdYlGn_r", "%.2e")
    _hm(axes[0, 1], icesd_grid,
        f"ice-sed interface drop [%]  (threshold = {pass_icesd_drop:.0f}%)",
        0, pass_icesd_drop * 2, "RdYlGn_r", "%.1f%%")
    _hm(axes[1, 0], snes_grid,
        f"Max SNES iters / step  (threshold = {PASS_SNES_ITERS})",
        0, PASS_SNES_ITERS * 1.2, "RdYlGn_r", "%d")

    pass_float = pass_grid.astype(float)
    im4 = axes[1, 1].imshow(pass_float, aspect="auto", cmap="RdYlGn",
                             vmin=0, vmax=1, origin="upper", interpolation="nearest")
    plt.colorbar(im4, ax=axes[1, 1])
    axes[1, 1].set_xticks(range(nk))
    axes[1, 1].set_xticklabels(x_labels, rotation=35, ha="right", fontsize=8)
    axes[1, 1].set_yticks(range(nd)); axes[1, 1].set_yticklabels(y_labels, fontsize=8)
    axes[1, 1].set_xlabel("k_pen", fontsize=10)
    axes[1, 1].set_ylabel("difvap_pen  [m²/s]", fontsize=10)
    n_pass = pass_grid.sum()
    axes[1, 1].set_title(f"PASS / FAIL  ({n_pass}/{nd * nk} pass)", fontsize=11)
    for ii in range(nd):
        for jj in range(nk):
            axes[1, 1].text(jj, ii, "✓" if pass_grid[ii, jj] else "✗",
                            ha="center", va="center", fontsize=11,
                            color="white" if pass_grid[ii, jj] else "black")

    plt.tight_layout()
    path = os.path.join(out_dir, fname)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Plot saved to: {path}")


def _plot(results: list[dict], difvap_sweep: list[float],
          kpen_sweep: list[float], tsf_sweep: list[float],
          out_dir: str, pass_air_drift: float, pass_icesd_drop: float) -> None:
    """Produce one heatmap per t_sed_freeze slice."""
    if len(tsf_sweep) == 1:
        _plot_slice(results, difvap_sweep, kpen_sweep, out_dir,
                    pass_air_drift, pass_icesd_drop,
                    title_extra="", fname="sweep_penalty2d.png")
    else:
        for tsf in tsf_sweep:
            slc = [r for r in results if r["t_sed_freeze"] == tsf]
            tag = f"{tsf:.2e}".replace("+", "")
            _plot_slice(slc, difvap_sweep, kpen_sweep, out_dir,
                        pass_air_drift, pass_icesd_drop,
                        title_extra=f"  |  t_sed_freeze = {tsf:.2e} s",
                        fname=f"sweep_penalty2d_tsf{tag}.png")


# ── CSV + summary ─────────────────────────────────────────────────────────────

def _save_csv(results: list[dict], out_dir: str) -> None:
    if not results:
        return
    path = os.path.join(out_dir, "sweep_penalty2d.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader()
        w.writerows(results)
    print(f"  CSV saved to:  {path}")


def _print_summary(results: list[dict], pass_air_drift: float,
                   pass_icesd_drop: float) -> None:
    print(f"\n{'─'*110}")
    print("  SUMMARY — difvap_pen × k_pen × t_sed_freeze sweep")
    print(f"  {'t_sed_frz':>10}  {'difvap_pen':>12}  {'k_pen':>12}  {'Δtot_air':>12}  "
          f"{'icesd_drop':>11}  {'SNES':>5}  {'bnd':>4}  {'PASS':>5}")
    print(f"  {'─'*102}")
    for r in results:
        dstr  = (f"{r['delta_tot_air']:.3e}"     if r["delta_tot_air"]    is not None else "N/A")
        idstr = (f"{r['ice_sed_drop_pct']:.1f}%" if r["ice_sed_drop_pct"] is not None else "N/A")
        snes  = r["max_snes_iters"] if r["max_snes_iters"] >= 0 else -1
        print(f"  {r['t_sed_freeze']:>10.2e}  {r['difvap_pen']:>12.2e}  {r['k_pen']:>12.3e}  "
              f"{dstr:>12}  {idstr:>11}  {snes:>5}  "
              f"{'✓' if r['bounds_ok'] else '✗':>4}  "
              f"{'✓' if r['PASS'] else '✗':>5}")

    passed = [r for r in results if r["PASS"]]
    if passed:
        best = min(passed, key=lambda r: abs(r["delta_tot_air"] or 1e30))
        print(f"\n  Best (smallest |Δtot_air| among passing):")
        print(f"    t_sed_freeze={best['t_sed_freeze']:.2e}  difvap_pen={best['difvap_pen']:.2e}"
              f"  k_pen={best['k_pen']:.3e}  Δtot_air={best['delta_tot_air']:.3e}")
        print(f"  {len(passed)}/{len(results)} combinations passed "
              f"(|Δtot_air|<{pass_air_drift:.1e}, ice_sed_drop<{pass_icesd_drop:.0f}%).")
    else:
        print("\n  No combinations passed all criteria. "
              "Consider relaxing --pass-air-drift or --pass-icesd-drop.")
    print(f"{'─'*110}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Independent 2D sweep of difvap_pen × k_pen for spurious-air suppression.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--binary",          default=BINARY_DEFAULT,
                   help="Path to the permafrost binary")
    p.add_argument("--opts",            default=OPTS_DEFAULT,
                   help="Opts file (t_final, t_sed_freeze, physics flags all come from here)")
    p.add_argument("--out-dir",         default=OUTDIR_DEFAULT,
                   help="Output directory for results")
    p.add_argument("--sweep-difvap",    nargs="+", type=float,
                   default=DIFVAP_SWEEP_DEFAULT, metavar="D",
                   help="difvap_pen values [m²/s] (default: 1e-7 1e-6 1e-5 1e-4 1e-3)")
    p.add_argument("--sweep-kpen",      nargs="+", type=float,
                   default=KPEN_SWEEP_DEFAULT, metavar="K",
                   help="k_pen values (default: 1e5 1e6 1e7 1e8 1e9)")
    p.add_argument("--sweep-tsf",       nargs="+", type=float,
                   default=TSF_SWEEP_DEFAULT, metavar="T",
                   help="t_sed_freeze values [s] (default: read single value from opts file)")
    p.add_argument("--timeout",         type=int,   default=TIMEOUT_DEFAULT,
                   help="Wall-clock timeout per run [s] (default: 600)")
    p.add_argument("--pass-air-drift",  type=float, default=PASS_AIR_DRIFT_DEFAULT,
                   help="Max |Δtot_air| to pass (default: 5e-6)")
    p.add_argument("--pass-icesd-drop", type=float, default=PASS_ICESD_DROP_DEFAULT,
                   help="Max ice-sed interface drop %% to pass (default: 20)")
    p.add_argument("--skip-build",      action="store_true",
                   help="Skip 'make clean && make' before the sweep")
    p.add_argument("--skip-postprocess", action="store_true",
                   help="Skip VTK/snapshot post-processing after each run")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Resolve t_sed_freeze sweep: CLI values take precedence; fall back to opts file
    tsf_sweep = (sorted(set(args.sweep_tsf))
                 if args.sweep_tsf is not None
                 else [_get_t_sed_freeze(args.opts)])

    total = len(args.sweep_difvap) * len(args.sweep_kpen) * len(tsf_sweep)
    print("=" * 70)
    print("  difvap_pen × k_pen × t_sed_freeze SWEEP")
    print(f"  binary:          {args.binary}")
    print(f"  opts:            {args.opts}")
    print(f"  out_dir:         {args.out_dir}")
    print(f"  difvap_pen:      {args.sweep_difvap}")
    print(f"  k_pen:           {args.sweep_kpen}")
    print(f"  t_sed_freeze:    {tsf_sweep}")
    print(f"  total runs:      {total}")
    print(f"  timeout:         {args.timeout} s/run")
    print(f"  pass_air_drift:  {args.pass_air_drift:.1e}")
    print(f"  pass_icesd_drop: {args.pass_icesd_drop:.0f} %")
    print(f"  post-processing: {'disabled' if args.skip_postprocess else 'enabled'}")
    print(f"  build:           {'skipped' if args.skip_build else 'make clean && make'}")
    print("=" * 70)

    if not args.skip_build:
        _build(ROOT)

    results = run_sweep(
        args.binary, args.opts,
        args.sweep_difvap, args.sweep_kpen, tsf_sweep,
        args.out_dir, args.timeout,
        pass_air_drift=args.pass_air_drift,
        pass_icesd_drop=args.pass_icesd_drop,
        skip_postprocess=args.skip_postprocess,
    )

    _print_summary(results, args.pass_air_drift, args.pass_icesd_drop)
    _save_csv(results, args.out_dir)
    _plot(results, args.sweep_difvap, args.sweep_kpen, tsf_sweep,
          args.out_dir, args.pass_air_drift, args.pass_icesd_drop)

    print(f"\n  All output saved to: {args.out_dir}")


if __name__ == "__main__":
    main()
