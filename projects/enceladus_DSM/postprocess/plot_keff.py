#!/usr/bin/env python3
"""plot_keff.py — effective thermal conductivity vs time and vs SSA.

    python3 plot_keff.py --dir <run> [--save-dir <dir>] [--baseline-days 1]

Writes two figures into <dir>/plots/ (or --save-dir):

    keff_time.png   k_xx, k_yy and k_iso = (k_xx + k_yy)/2 against time
    keff_ssa.png    the same three against specific surface area

WHICH CSV. An in-line run writes k_eff.csv. A -keff_replay writes
k_eff_<law>.csv beside it, and between 2026-09-23 and 2026-09-26 in-line runs
also wrote k_eff_tensor.csv, because tensor became the default. When several
exist the tensor file wins, since tensor is the current law. The law used is
read from outp.txt ("band interpolation: ...") and printed on the figure, so
an arith k_eff.csv from before 2026-09-23 is labelled as such rather than
passing as tensor.

SSA. Column 0 of SSA_evo.dat is int phi^2(1-phi)^2 dV / eps, and each unit of
interface length contributes eps/6 to that integral (see plot_ssa.py), so

    SSA = 6 * column0 / (Lx * Ly)          [1/m]

This is interface length per unit CELL area, the convention
studies/keff_sintering/analyze_pilot.py uses. The ice fraction is conserved
(phi_bar is constant to seven figures in every pilot run), so SSA per unit ice
area is this divided by a constant. The shape of k_eff(SSA) does not depend on
which convention is used.

THE FIRST DAY IS GREYED OUT. The initial condition is an analytic sum of tanh
profiles, not an equilibrated phase field, so the first hours are the field
relaxing, not sintering (CAMPAIGN.md "The baseline is t = 1 day"). Those
samples are drawn in grey on both figures and never used as a baseline.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless on HPC
import matplotlib.pyplot as plt

from pplib import load_ssa, opt_float, read_opts

INTERFACE_FACTOR = 6.0          # int phi^2(1-phi)^2 dx = eps/6; see plot_ssa.py
DAY = 86400.0

# Okabe-Ito, same order as preprocess/figstyle.py C[0], C[1]: CVD-safe.
C_XX, C_YY = "#0072B2", "#D55E00"
C_ISO = "#1a1a1a"               # the headline series, in ink
C_RELAX = "#b8b8b8"             # IC-relaxation samples


def find_keff_csv(run: Path):
    """(path, law) for the k_eff series to plot, or (None, None)."""
    law = None
    outp = run / "outp.txt"
    if outp.is_file():
        with open(outp, errors="replace") as fh:
            for line in fh:
                m = re.search(r"band interpolation:\s*(\w+)", line)
                if m:
                    law = m.group(1)
                    break
    for name, file_law in (("k_eff_tensor.csv", "tensor"),
                           ("k_eff.csv", law),
                           ("k_eff_sharp.csv", "sharp")):
        p = run / name
        if p.is_file():
            # k_eff.csv from before the law was logged: runs before
            # 2026-09-23 are arith (the only law that existed then).
            return p, (file_law or "arith, assumed (no law in outp.txt)")
    return None, None


def load(run: Path):
    kf, law = find_keff_csv(run)
    if kf is None:
        return None
    k = np.atleast_1d(np.genfromtxt(kf, delimiter=",", names=True))
    ssa = load_ssa(str(run))
    if ssa is None or len(k) == 0:
        return None
    opts = read_opts(str(run))
    lx = opt_float(opts, "-Lx", None)
    ly = opt_float(opts, "-Ly", None)
    if not lx or not ly:
        return None

    # Match on STEP: the k_eff sample and the SSA row come from the same
    # accepted step, so no interpolation is needed. Fall back to time if a
    # step is missing (e.g. a merge renumbered one file but not the other).
    ssa_step = ssa[:, 3].astype(int)
    row_of = {s: i for i, s in enumerate(ssa_step)}
    idx = np.array([row_of.get(int(s), -1) for s in k["step"]])
    miss = idx < 0
    if miss.any():
        idx[miss] = [int(np.argmin(np.abs(ssa[:, 2] - t))) for t in k["time"][miss]]
    s = INTERFACE_FACTOR * ssa[idx, 0] / (lx * ly)

    return {"t": k["time"], "kxx": k["k_00"], "kyy": k["k_11"], "kiso": k["k_iso"],
            "ssa": s, "law": law, "csv": kf.name}


def _series(ax, x, d, relax, xlabel, direct_labels=True):
    live = ~relax
    for key, col, lw, lab in (("kxx", C_XX, 1.6, r"$k_{xx}$"),
                              ("kyy", C_YY, 1.6, r"$k_{yy}$"),
                              ("kiso", C_ISO, 2.4, r"$k_\mathrm{iso}$")):
        y = d[key]
        if relax.any():
            ax.plot(x[relax], y[relax], "-", color=C_RELAX, lw=lw, zorder=1)
        ax.plot(x[live], y[live], "-", color=col, lw=lw, label=lab, zorder=2)
        # Direct label at the right-hand end of each live curve -- only where
        # that end is the late-time end, i.e. against time. Against SSA the
        # right-hand end of the live curve abuts the grey relaxation segment.
        if live.any() and direct_labels:
            ax.annotate(lab, (x[live][-1], y[live][-1]), xytext=(6, 0),
                        textcoords="offset points", va="center", fontsize=12,
                        color="#333333")
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]", fontsize=14)
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: cwd)")
    p.add_argument("--save-dir", default=None,
                   help="where the two PNGs go (default: <dir>/plots)")
    p.add_argument("--baseline-days", type=float, default=1.0,
                   help="samples before this are IC relaxation and drawn grey (default 1)")
    a = p.parse_args(argv)

    run = Path(a.dir)
    d = load(run)
    if d is None:
        print(f"  no k_eff CSV (or no SSA_evo.dat / -Lx -Ly) in {run}; nothing to plot")
        return 0                   # not an error: most runs do not pass -keff

    out = Path(a.save_dir) if a.save_dir else run / "plots"
    os.makedirs(out, exist_ok=True)
    relax = d["t"] < a.baseline_days * DAY
    tday = d["t"] / DAY
    n_live = int((~relax).sum())
    sub = (f"{len(d['t'])} samples ({d['csv']}, {d['law']} law); "
           f"grey = first {a.baseline_days:g} d, IC relaxation")

    # --- k_eff vs time ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    if relax.any():
        ax.axvspan(0, a.baseline_days, color="#f0f0f0", zorder=0, lw=0)
    _series(ax, tday, d, relax, "Time [d]")
    ax.legend(fontsize=11, loc="lower right", frameon=False)
    i0 = int(np.argmax(~relax)) if n_live else 0
    rise = (d["kiso"][-1] / d["kiso"][i0] - 1) * 100
    ax.set_title(f"Effective thermal conductivity vs time\n"
                 f"$k_\\mathrm{{iso}}$ {d['kiso'][i0]:.4f} → {d['kiso'][-1]:.4f} "
                 f"({rise:+.1f}% from t = {tday[i0]:.2f} d)", fontsize=15)
    fig.text(0.01, 0.005, sub, fontsize=9, color="#555555")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    f1 = out / "keff_time.png"
    fig.savefig(f1, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- k_eff vs SSA ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    _series(ax, d["ssa"], d, relax, r"SSA  [m$^{-1}$]  (interface length per cell area)",
            direct_labels=False)
    ax.legend(fontsize=11, loc="upper right", frameon=False)
    # Time runs right-to-left here (SSA falls as the packing sinters).
    ax.annotate("", xy=(0.12, 0.93), xytext=(0.30, 0.93), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#555555", lw=1.2))
    ax.text(0.31, 0.93, "time", transform=ax.transAxes, va="center",
            fontsize=11, color="#555555")
    ax.set_title(f"Effective thermal conductivity vs specific surface area\n"
                 f"SSA {d['ssa'][i0]:.4g} → {d['ssa'][-1]:.4g} m$^{{-1}}$ "
                 f"from t = {tday[i0]:.2f} d", fontsize=15)
    fig.text(0.01, 0.005, sub, fontsize=9, color="#555555")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    f2 = out / "keff_ssa.png"
    fig.savefig(f2, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  k_eff: {len(d['t'])} samples from {d['csv']} ({d['law']}), "
          f"{n_live} after {a.baseline_days:g} d; k_iso rise {rise:+.1f}%")
    print(f"  wrote {f1}")
    print(f"  wrote {f2}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
