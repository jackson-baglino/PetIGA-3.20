#!/usr/bin/env python3
"""plot_keff.py — effective thermal conductivity vs time and vs SSA.

    python3 plot_keff.py --dir <run> [--save-dir <dir>] [--baseline-days 1]

Writes four figures under <dir>/plots/keff/ (or --save-dir):

    absolute/keff_time.png     k_xx, k_yy, k_iso = (k_xx + k_yy)/2 vs time [d]
    absolute/keff_ssa.png      the same three vs SSA [1/m]
    absolute/keff_offdiag_time.png
                               k_xy and k_yx vs time, and |k_xy| relative to
                               k_xx and k_yy -- a check that they are small
    normalized/keff_time.png   k / k_b  vs  t / tau_sub
    normalized/keff_ssa.png    k / k_b  vs  SSA / SSA_b

NORMALIZATION. k and SSA are divided by their value AT the baseline time,
interpolated to t = 11 tau_sub (--baseline-tau; = 1 d at -20 C, the pilot's
reference, and the compare_keff.py convention), not at t = 0: t = 0 is the
unrelaxed initial condition (below). Each component is divided by its own
baseline value, so k_xx/k_xx,b, k_yy/k_yy,b and k_iso/k_iso,b all start at 1
and the plot shows the relative rise. Pass --baseline-days 0 to normalize by
the t = 0 values instead.

Time has no initial value to divide by, so it is made dimensionless with
tau_sub, the solver's interface-kinetic timescale (logged in outp.txt). It is
temperature-dependent, so t/tau_sub is the natural axis for putting different
temperatures on one plot. Whether it actually collapses them is something the
plot shows; it is not assumed. With no tau_sub in outp.txt the normalized
time axis falls back to t / t_b.

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
samples are drawn in grey on every figure and never used as a baseline.
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


SERIES = (("kxx", C_XX, 1.6, r"$k_{xx}$"),
          ("kyy", C_YY, 1.6, r"$k_{yy}$"),
          ("kiso", C_ISO, 2.4, r"$k_\mathrm{iso}$"))


def _series(ax, x, ys, relax, xlabel, ylabel, direct_labels=True):
    live = ~relax
    for key, col, lw, lab in SERIES:
        y = ys[key]
        if relax.any():
            # Run the grey segment through the first live sample so the two
            # join. Every point on it is a measured sample, matched to its
            # SSA row by step; nothing is interpolated or extrapolated.
            grey = relax.copy()
            if live.any():
                grey[int(np.argmax(live))] = True
            ax.plot(x[grey], y[grey], "-", color=C_RELAX, lw=lw, zorder=1,
                    label="first day, IC relaxation (measured)" if key == "kiso" else None)
        ax.plot(x[live], y[live], "-", color=col, lw=lw, label=lab, zorder=2)
        # Direct label at the right-hand end of each live curve -- only where
        # that end is the late-time end, i.e. against time. Against SSA the
        # right-hand end of the live curve abuts the grey relaxation segment.
        if live.any() and direct_labels:
            ax.annotate(lab, (x[live][-1], y[live][-1]), xytext=(6, 0),
                        textcoords="offset points", va="center", fontsize=12,
                        color="#333333")
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def _legend(ax, loc):
    """Legend with the k series first and the grey relaxation entry last."""
    h, l = ax.get_legend_handles_labels()
    order = sorted(range(len(l)), key=lambda i: "relaxation" in l[i])
    ax.legend([h[i] for i in order], [l[i] for i in order],
              fontsize=11, loc=loc, frameon=False)


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
    p.add_argument("--baseline-tau", type=float, default=86400.0 / 7822.3,
                   help="baseline as t/tau_sub (default 11.05 = 1 d at -20 C). "
                        "Samples before it are IC relaxation, drawn grey; values "
                        "interpolated to it are the normalization baseline")
    p.add_argument("--baseline-days", type=float, default=None,
                   help="baseline in days instead (also the fallback, 1 d, when "
                        "outp.txt has no tau_sub)")
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

    # Baseline in units of tau_sub when the run logs it (11 = 1 d at -20 C),
    # so every temperature is normalized at the same stage of sintering -- the
    # compare_keff.py convention. 1 d at -40 C is only 1.4 tau_sub, still
    # inside the IC relaxation, and read as a +40% rise. Days is the fallback.
    tau_b = read_tau_sub(run)
    t_base = (a.baseline_tau * tau_b if (tau_b and a.baseline_days is None)
              else (a.baseline_days if a.baseline_days is not None else 1.0) * DAY)
    b_days = t_base / DAY
    relax = d["t"] < t_base
    if relax.all():
        print(f"  run ends before the {b_days:.3g} d baseline; "
              "normalizing by the last sample")
    ib = int(np.argmax(~relax)) if not relax.all() else len(relax) - 1
    # The baseline is the value AT t_b, interpolated -- not the first sample
    # past it. With k_eff every 5 steps the first sample past 1 d lands at
    # ~1.4 d at -20 C, which understated the rise by ~3 points and disagreed
    # with compare_keff.py (2026-09-29).
    tb = t_base if not relax.all() else d["t"][-1]
    base = {key: float(np.interp(tb, d["t"], d[key])) for key in ("kxx", "kyy", "kiso", "ssa")}
    tday = d["t"] / DAY
    rise = (d["kiso"][-1] / base["kiso"] - 1) * 100
    note = (f"{len(d['t'])} samples ({d['csv']}, {d['law']} law); "
            f"grey = before {b_days:.3g} d, IC relaxation")
    k_label = r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]"
    written = []

    # --- absolute ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    if relax.any():
        ax.axvspan(0, b_days, color="#f0f0f0", zorder=0, lw=0)
    _series(ax, tday, d, relax, "Time [d]", k_label)
    _legend(ax, "lower right")
    ax.set_title(f"Effective thermal conductivity vs time\n"
                 f"$k_\\mathrm{{iso}}$ {base['kiso']:.4f} → {d['kiso'][-1]:.4f} "
                 f"({rise:+.1f}% from t = {tb / DAY:.2f} d)", fontsize=15)
    _save(fig, out_abs / "keff_time.png", note)
    written.append(out_abs / "keff_time.png")

    fig, ax = plt.subplots(figsize=(10, 6))
    _series(ax, d["ssa"], d, relax,
            r"SSA  [m$^{-1}$]  (interface length per cell area)", k_label,
            direct_labels=False)
    _legend(ax, "upper right")
    _time_arrow(ax)
    ax.set_title(f"Effective thermal conductivity vs specific surface area\n"
                 f"SSA {base['ssa']:.4g} → {d['ssa'][-1]:.4g} m$^{{-1}}$ "
                 f"from t = {tb / DAY:.2f} d", fontsize=15)
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
        if relax.any():
            ax.axvspan(0, b_days, color="#f0f0f0", zorder=0, lw=0)
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
    ax2.plot(tday, r_xx, "-", color=C_XX, lw=1.8, label=r"$|k_{xy}|\,/\,k_{xx}$")
    ax2.plot(tday, r_yy, "-", color=C_YY, lw=1.8, label=r"$|k_{xy}|\,/\,k_{yy}$")
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
        tn, t_label, t_ref = d["t"] / tb, r"$t\,/\,t_b$", "no tau_sub in outp.txt"
    k_nlabel = r"$k\,/\,k_b$"
    bnote = (f"baseline b = values interpolated to t_b = {tb / DAY:.2f} d; {t_ref}")

    fig, ax = plt.subplots(figsize=(10, 6))
    if relax.any():
        ax.axvspan(0, (tb / tau) if tau else 1.0, color="#f0f0f0", zorder=0, lw=0)
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    _series(ax, tn, kn, relax, t_label, k_nlabel)
    _legend(ax, "lower right")
    ax.set_title(f"Normalized effective thermal conductivity vs normalized time\n"
                 f"$k_\\mathrm{{iso}}/k_b$ → {kn['kiso'][-1]:.3f}  at  "
                 f"{t_label} = {tn[-1]:.4g}", fontsize=15)
    _save(fig, out_norm / "keff_time.png", note + "\n" + bnote)
    written.append(out_norm / "keff_time.png")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    ax.axvline(1.0, color="#999999", lw=0.8, ls=":")
    _series(ax, sn, kn, relax, r"SSA$\,/\,$SSA$_b$", k_nlabel, direct_labels=False)
    _legend(ax, "upper right")
    _time_arrow(ax)
    ax.set_title(f"Normalized effective thermal conductivity vs normalized SSA\n"
                 f"SSA/SSA$_b$ → {sn[-1]:.3f},  $k_\\mathrm{{iso}}/k_b$ → "
                 f"{kn['kiso'][-1]:.3f}", fontsize=15)
    _save(fig, out_norm / "keff_ssa.png", note + "\n" + bnote)
    written.append(out_norm / "keff_ssa.png")

    print(f"  k_eff: {len(d['t'])} samples from {d['csv']} ({d['law']}); baseline "
          f"t_b = {tb / DAY:.2f} d; k_iso rise {rise:+.1f}%; "
          f"tau_sub = {tau if tau else 'n/a'} s")
    for w in written:
        print(f"  wrote {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
