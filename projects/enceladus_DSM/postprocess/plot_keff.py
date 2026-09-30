#!/usr/bin/env python3
"""plot_keff.py — effective thermal conductivity vs time and vs SSA.

    python3 plot_keff.py --dir <run> [--save-dir <dir>]

Writes four figures under <dir>/plots/keff/ (or --save-dir):

    absolute/keff_time.png     k_xx, k_yy, k_iso = (k_xx + k_yy)/2 vs time [d]
    absolute/keff_ssa.png      the same three vs SSA [1/m]
    absolute/keff_offdiag_time.png
                               k_xy and k_yx vs time, and |k_xy| relative to
                               k_xx and k_yy -- a check that they are small
    normalized/keff_time.png   k_eff / k_eff,0  vs  t / tau_sub
    normalized/keff_ssa.png    k_eff / k_eff,0  vs  SSA / SSA_0

THE OPENING SAMPLE IS t = 0, the same rule as plot_keff_snapshots.py and the
movies (pplib.opening_step): the first k_eff sample with 1 s <= t <= 1 h. By
1 s the vapour field has relaxed onto the ice geometry; the IC itself (t = 0)
still carries a uniform vapour field. Samples before it are not drawn, and
nothing after it is greyed out: the fast early rise as the packing relaxes is
data, explained in the manuscript rather than annotated on the figure.

NORMALIZATION. Each quantity is divided by its own MEASURED value at that
opening sample -- k_xx by k_xx,0, k_yy by k_yy,0, k_iso by k_iso,0, SSA by
SSA_0 -- never by an interpolated one. (Until 2026-09-30 this divided by values
interpolated to t = 11 tau_sub, labelled "b".) Every curve runs from the
opening sample to the last one; nothing is extrapolated.

LINES. k_xx dashed, k_yy dotted, k_iso solid, in one cool family (below):
normalized, the three often lie on top of one another, and the dash pattern
is what keeps all three visible.

Time has no initial value to divide by, so it is made dimensionless with
tau_sub, the solver's interface-kinetic timescale (logged in outp.txt). It is
temperature-dependent, so t/tau_sub is the natural axis for putting different
temperatures on one plot. Whether it actually collapses them is something the
plot shows; it is not assumed. With no tau_sub in outp.txt the normalized
time axis falls back to days.

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

from pplib import load_ssa, opening_step, opt_float, read_opts

INTERFACE_FACTOR = 6.0          # int phi^2(1-phi)^2 dx = eps/6; see plot_ssa.py
DAY = 86400.0

# One family, sampled from cmocean `deep` at 0.38 / 0.64 / 0.92: seafoam,
# steel blue, deep indigo. Single-run curves take a COOL map because the
# parameter sweeps own the others -- temperature is cmocean `thermal`,
# porosity an amp-to-black map -- so a fixed-parameter figure never reads as
# one point of a sweep. k_iso, the headline, is the darkest (13.4:1 on white).
# Checked (Machado 2009 CVD, OKLab dE x100): every pair >= 17.9 under protan,
# deutan and tritan, >= 19.7 normal; k_xx is 2.7:1 on white, so it is never
# identified by colour alone -- the legend names it.
C_XX, C_YY = "#55ada3", "#3e6b96"
C_ISO = "#352949"
LS_XX, LS_YY, LS_ISO = (0, (5, 2.5)), (0, (1, 1.8)), "-"   # dashed, dotted, solid


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
    # accepted step. A k_eff sample with no SSA row at its step is DROPPED,
    # not paired with the nearest-in-time row: only measured pairs are shown.
    ssa_step = ssa[:, 3].astype(int)
    row_of = {s: i for i, s in enumerate(ssa_step)}
    idx = np.array([row_of.get(int(s), -1) for s in k["step"]])
    keep = idx >= 0
    if not keep.all():
        print(f"  {int((~keep).sum())} k_eff sample(s) have no SSA_evo.dat row at "
              "their step; dropped")
    k, idx = k[keep], idx[keep]
    if len(k) == 0:
        return None
    s = INTERFACE_FACTOR * ssa[idx, 0] / (lx * ly)

    return {"t": k["time"], "step": k["step"].astype(int),
            "kxx": k["k_00"], "kyy": k["k_11"], "kiso": k["k_iso"],
            "kxy": k["k_01"], "kyx": k["k_10"],
            "ssa": s, "law": law, "csv": kf.name}


def opening_index(d):
    """Index of the opening sample: the first with 1 s <= t <= 1 h
    (pplib.opening_step), else the first with t > 0, else 0."""
    op = opening_step(d["step"], d["t"])
    if op is not None:
        return int(np.flatnonzero(d["step"] == op)[0])
    pos = np.flatnonzero(d["t"] > 0)
    return int(pos[0]) if pos.size else 0


def from_opening(d):
    """d with every sample before the opening one dropped."""
    i0 = opening_index(d)
    return {k: (v[i0:] if isinstance(v, np.ndarray) else v) for k, v in d.items()}


def read_tau_sub(run: Path):
    """tau_sub [s] from outp.txt's parameter table, or None."""
    outp = run / "outp.txt"
    if not outp.is_file():
        return None
    with open(outp, errors="replace") as fh:
        for line in fh:
            m = re.match(r"\s*tau_sub\s+([0-9.eE+-]+)\s*s", line)
            if m:
                return float(m.group(1))
    return None


SERIES = (("kxx", C_XX, LS_XX, 1.6, r"$k_{xx}$"),
          ("kyy", C_YY, LS_YY, 1.8, r"$k_{yy}$"),
          ("kiso", C_ISO, LS_ISO, 2.2, r"$k_\mathrm{iso}$"))


def _series(ax, x, ys, xlabel, ylabel, direct_labels=True):
    """Every sample from the opening one to the last, measured pairs only."""
    for key, col, ls, lw, lab in SERIES:
        y = ys[key]
        # k_iso UNDER the others: where they coincide, the dashes show on it.
        ax.plot(x, y, color=col, ls=ls, lw=lw, label=lab,
                zorder={"kiso": 2, "kyy": 3, "kxx": 4}[key],
                dash_capstyle="round")
        # Direct label at the late-time end -- against time only; against SSA
        # the right-hand end is t = 0.
        if direct_labels:
            ax.annotate(lab, (x[-1], y[-1]), xytext=(6, 0),
                        textcoords="offset points", va="center", fontsize=12,
                        color="#333333")
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def _legend(ax, loc):
    ax.legend(fontsize=11, loc=loc, frameon=False, handlelength=2.6)


def _time_arrow(ax):
    """SSA falls as the packing sinters, so time runs right-to-left."""
    ax.annotate("", xy=(0.12, 0.93), xytext=(0.30, 0.93), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", color="#555555", lw=1.2))
    ax.text(0.31, 0.93, "time", transform=ax.transAxes, va="center",
            fontsize=11, color="#555555")


def _save(fig, path, note):
    fig.text(0.01, 0.005, note, fontsize=9, color="#555555")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: cwd)")
    p.add_argument("--save-dir", default=None,
                   help="root for absolute/ and normalized/ (default: <dir>/plots/keff)")
    a = p.parse_args(argv)

    run = Path(a.dir)
    d = load(run)
    if d is None:
        print(f"  no k_eff CSV (or no SSA_evo.dat / -Lx -Ly) in {run}; nothing to plot")
        return 0                   # not an error: most runs do not pass -keff

    root = Path(a.save_dir) if a.save_dir else run / "plots" / "keff"
    out_abs, out_norm = root / "absolute", root / "normalized"
    os.makedirs(out_abs, exist_ok=True)
    os.makedirs(out_norm, exist_ok=True)

    d = from_opening(d)
    base = {key: float(d[key][0]) for key in ("kxx", "kyy", "kiso", "ssa")}
    t0 = float(d["t"][0])
    tday = d["t"] / DAY
    rise = (d["kiso"][-1] / base["kiso"] - 1) * 100
    note = (f"{len(d['t'])} samples ({d['csv']}, {d['law']} law); "
            f"opening sample t_0 = {t0:.4g} s (step {d['step'][0]})")
    k_label = r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]"
    written = []

    # --- absolute ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    _series(ax, tday, d, "Time [d]", k_label)
    _legend(ax, "lower right")
    ax.set_title(f"Effective thermal conductivity vs time\n"
                 f"$k_\\mathrm{{iso}}$ {base['kiso']:.4f} → {d['kiso'][-1]:.4f} "
                 f"({rise:+.1f}% from t$_0$ = {t0:.4g} s)", fontsize=15)
    _save(fig, out_abs / "keff_time.png", note)
    written.append(out_abs / "keff_time.png")

    fig, ax = plt.subplots(figsize=(10, 6))
    _series(ax, d["ssa"], d,
            r"SSA  [m$^{-1}$]  (interface length per cell area)", k_label,
            direct_labels=False)
    _legend(ax, "upper right")
    _time_arrow(ax)
    ax.set_title(f"Effective thermal conductivity vs specific surface area\n"
                 f"SSA {base['ssa']:.4g} → {d['ssa'][-1]:.4g} m$^{{-1}}$ "
                 f"from t$_0$ = {t0:.4g} s", fontsize=15)
    _save(fig, out_abs / "keff_ssa.png", note)
    written.append(out_abs / "keff_ssa.png")

    # --- off-diagonal, absolute, vs time --------------------------------------
    # Two panels, not one axis: k_xy is ~10x smaller than k_xx, so on a shared
    # axis it would sit on the floor. (a) the components themselves, (b) their
    # size relative to each diagonal entry, which is the question being asked.
    # k_xy and k_yx should coincide -- the homogenized tensor is symmetric --
    # so the asymmetry is printed as a check on the corrector solve.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for ax in (ax1, ax2):
        ax.grid(True, alpha=0.25, lw=0.6)
        ax.tick_params(labelsize=11)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    ax1.axhline(0.0, color="#999999", lw=0.8, ls=":")
    ax1.plot(tday, d["kxy"], "-", color=C_XX, lw=2.2, label=r"$k_{xy}$")
    ax1.plot(tday, d["kyx"], "--", color=C_YY, lw=1.4, label=r"$k_{yx}$")
    ax1.set_ylabel(r"$k_{xy}$, $k_{yx}$  [W m$^{-1}$ K$^{-1}$]", fontsize=13)
    ax1.legend(fontsize=11, loc="best", frameon=False)
    r_xx = np.abs(d["kxy"]) / d["kxx"] * 100
    r_yy = np.abs(d["kxy"]) / d["kyy"] * 100
    ax2.plot(tday, r_xx, color=C_XX, ls=LS_XX, lw=1.8, label=r"$|k_{xy}|\,/\,k_{xx}$")
    ax2.plot(tday, r_yy, color=C_YY, ls=LS_YY, lw=2.0, label=r"$|k_{xy}|\,/\,k_{yy}$")
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel("relative to diagonal  [%]", fontsize=13)
    ax2.set_xlabel("Time [d]", fontsize=14)
    ax2.legend(fontsize=11, loc="best", frameon=False)
    asym = float(np.max(np.abs(d["kxy"] - d["kyx"])))
    ax1.set_title(f"Off-diagonal conductivity vs time\n"
                  f"max $|k_{{xy}}|/k_{{yy}}$ = {r_yy.max():.2f}%,  "
                  f"max $|k_{{xy}} - k_{{yx}}|$ = {asym:.1e} W m$^{{-1}}$ K$^{{-1}}$",
                  fontsize=15)
    _save(fig, out_abs / "keff_offdiag_time.png", note)
    written.append(out_abs / "keff_offdiag_time.png")
    print(f"  off-diagonal: max |k_xy|/k_xx = {r_xx.max():.2f}%, "
          f"max |k_xy|/k_yy = {r_yy.max():.2f}%, max |k_xy - k_yx| = {asym:.2e}")

    # --- normalized ----------------------------------------------------------
    kn = {key: d[key] / base[key] for key, *_ in SERIES}
    sn = d["ssa"] / base["ssa"]
    tau = read_tau_sub(run)
    if tau:
        tn, t_label, t_ref = d["t"] / tau, r"$t\,/\,\tau_\mathrm{sub}$", \
            f"$\\tau_\\mathrm{{sub}}$ = {tau:.4g} s"
    else:
        tn, t_label, t_ref = tday, "Time [d]", "no tau_sub in outp.txt"
    k_nlabel = r"$k_\mathrm{eff}\,/\,k_{\mathrm{eff},0}$"
    bnote = f"subscript 0 = the opening sample, t_0 = {t0:.4g} s; {t_ref}"

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    # No direct labels: normalized, the three curves often coincide.
    _series(ax, tn, kn, t_label, k_nlabel, direct_labels=False)
    _legend(ax, "lower right")
    ax.set_title(f"Normalized effective thermal conductivity vs normalized time\n"
                 f"$k_\\mathrm{{iso}}/k_{{\\mathrm{{iso}},0}}$ → {kn['kiso'][-1]:.3f}  at  "
                 f"{t_label} = {tn[-1]:.4g}", fontsize=15)
    _save(fig, out_norm / "keff_time.png", note + "\n" + bnote)
    written.append(out_norm / "keff_time.png")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    ax.axvline(1.0, color="#999999", lw=0.8, ls=":")
    _series(ax, sn, kn, r"SSA$\,/\,$SSA$_0$", k_nlabel, direct_labels=False)
    _legend(ax, "upper right")
    _time_arrow(ax)
    ax.set_title(f"Normalized effective thermal conductivity vs normalized SSA\n"
                 f"SSA/SSA$_0$ → {sn[-1]:.3f},  $k_\\mathrm{{iso}}/k_{{\\mathrm{{iso}},0}}$ → "
                 f"{kn['kiso'][-1]:.3f}", fontsize=15)
    _save(fig, out_norm / "keff_ssa.png", note + "\n" + bnote)
    written.append(out_norm / "keff_ssa.png")

    print(f"  k_eff: {len(d['t'])} samples from {d['csv']} ({d['law']}); opening "
          f"sample t_0 = {t0:.4g} s; k_iso rise {rise:+.1f}%; "
          f"tau_sub = {tau if tau else 'n/a'} s")
    for w in written:
        print(f"  wrote {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
