#!/usr/bin/env python3
"""compare_keff.py — overlay k_eff curves across temperatures and porosities.

    python3 compare_keff.py <root> [<root> ...] [--out <dir>] [--baseline-days 1]

Finds every run directory under the roots that holds a k_eff CSV, reads its
porosity, temperature and seed from the directory name
(``..._phi0.325_..._seed3_..._T-20__...``), and writes two families of figures:

    by_phi/phi<X>/   fixed porosity, one curve per temperature
    by_T/T<Y>/       fixed temperature, one curve per porosity

each with absolute/ and normalized/ versions of

    keff_time.png    k_iso against time        (normalized: k/k_b vs t/tau_sub)
    keff_ssa.png     k_iso against SSA         (normalized: k/k_b vs SSA/SSA_b)

A group is drawn only when it holds at least two values of the varied
parameter; the script says which groups it skipped and why.

SEEDS. Each curve is the mean over the seeds of that condition, and the band
around it is the seed min..max, so the spread between conditions can be read
against the spread within one. The mean is taken at common times (time plots),
and against SSA the mean k is plotted at the mean SSA of those same times --
both coordinates averaged at matched t, so no seed's SSA axis is resampled.

BASELINE, IN UNITS OF tau_sub. Not t = 0, because the first hours are the
initial condition relaxing (plot_keff.py). And not a fixed TIME either: on the
2026-09-25 warm-end batch a run reaches a given SSA exactly tau_sub(T)/tau_sub(-20)
times sooner (1.59x at -15 C, 2.47-2.50x at -10 C, every seed), so 1 d is
11 tau_sub at -20 C but 27 at -10 C -- a later stage of sintering. Normalizing
there offset the curves by the baseline, not by the physics. Every run is
therefore normalized at the SAME t/tau_sub: --baseline-days at the slowest
condition in the comparison (default 1 d there), i.e. at an earlier time for
the warmer runs. The IC relaxation is the same dynamics, so it scales the same
way. Pass --baseline-tau to set it directly.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_keff import DAY, load, read_tau_sub          # noqa: E402

# One hue per varied parameter, light -> dark with the parameter value: these
# are ordered magnitudes (sequential), not identities. The light end is kept
# at 0.45 so the palest curve still reads on white.
CMAP = {"T": "Oranges", "phi": "Blues"}
N_GRID = 400


def parse(run: Path):
    n = run.name
    phi = re.search(r"phi([0-9.]+?)_", n)
    T = re.search(r"_T(-?\d+)__", n) or re.search(r"_T(-?\d+)", n)
    seed = re.search(r"seed(\d+)", n)
    if not (phi and T and seed):
        return None
    return float(phi.group(1)), int(T.group(1)), int(seed.group(1))


def discover(roots):
    runs = []
    for root in roots:
        for kf in Path(root).rglob("k_eff*.csv"):
            d = kf.parent
            if d in {r["dir"] for r in runs}:
                continue
            meta = parse(d)
            if meta is None:
                print(f"  skip {d.name}: no phi/T/seed in the name")
                continue
            data = load(d)
            if data is None:
                continue
            data.update(dir=d, phi=meta[0], T=meta[1], seed=meta[2],
                        tau=read_tau_sub(d))
            runs.append(data)
    return runs


def condition_mean(runs, baseline_tau):
    """Seed-mean curves for one condition on a common log-spaced time grid."""
    taus = {r["tau"] for r in runs}
    tau = taus.pop() if len(taus) == 1 else None
    baseline_s = baseline_tau * tau
    t_end = min(r["t"][-1] for r in runs)
    t_first = max(r["t"][1] for r in runs)      # skip t = 0 for the log grid
    tg = np.unique(np.concatenate([[0.0], np.geomspace(t_first, t_end, N_GRID),
                                   [baseline_s]]))
    tg = tg[tg <= t_end]
    K = np.array([np.interp(tg, r["t"], r["kiso"]) for r in runs])
    S = np.array([np.interp(tg, r["t"], r["ssa"]) for r in runs])
    ib = int(np.searchsorted(tg, baseline_s))
    Kn, Sn = K / K[:, [ib]], S / S[:, [ib]]
    return dict(t=tg, ib=ib, tau=tau, n=len(runs),
                k=K.mean(0), klo=K.min(0), khi=K.max(0), s=S.mean(0),
                kn=Kn.mean(0), knlo=Kn.min(0), knhi=Kn.max(0), sn=Sn.mean(0))


def _style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def _curve(ax, x, y, lo, hi, ib, color, label, band=True):
    """Live part solid with a seed band; the relaxation part thin and faded."""
    ax.plot(x[:ib + 1], y[:ib + 1], "-", color=color, lw=1.0, alpha=0.45)
    ax.plot(x[ib:], y[ib:], "-", color=color, lw=2.2, label=label)
    if band:
        ax.fill_between(x[ib:], lo[ib:], hi[ib:], color=color, alpha=0.15, lw=0)


def draw_group(conds, vary, fixed_txt, out: Path, baseline_tau):
    """conds: {value: condition_mean dict}. vary: 'T' or 'phi'."""
    vals = sorted(conds)
    cmap = plt.get_cmap(CMAP[vary])
    col = {v: cmap(0.45 + 0.5 * i / max(1, len(vals) - 1)) for i, v in enumerate(vals)}
    lab = (lambda v: f"T = {v} °C") if vary == "T" else (lambda v: f"φ = {v:g}")
    nseed = sorted({c["n"] for c in conds.values()})
    note = (f"line = mean over {'/'.join(map(str, nseed))} seeds, band = seed min–max; "
            f"thin faded segment = before t = {baseline_tau:.3g} tau_sub (IC relaxation)")
    k_lab = r"$k_\mathrm{iso}$  [W m$^{-1}$ K$^{-1}$]"
    written = []

    def save(fig, path, extra=""):
        fig.text(0.01, 0.005, note + extra, fontsize=9, color="#555555")
        fig.tight_layout(rect=(0, 0.035, 1, 1))
        os.makedirs(path.parent, exist_ok=True)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    # absolute vs time
    fig, ax = plt.subplots(figsize=(10, 6))
    for v in vals:
        c = conds[v]
        _curve(ax, c["t"] / DAY, c["k"], c["klo"], c["khi"], c["ib"], col[v], lab(v))
    _style(ax, "Time [d]", k_lab)
    ax.legend(fontsize=11, loc="lower right", frameon=False)
    ax.set_title(f"$k_\\mathrm{{iso}}$ vs time, {fixed_txt}", fontsize=15)
    save(fig, out / "absolute" / "keff_time.png")

    # absolute vs SSA
    fig, ax = plt.subplots(figsize=(10, 6))
    for v in vals:
        c = conds[v]
        _curve(ax, c["s"], c["k"], c["klo"], c["khi"], c["ib"], col[v], lab(v))
    _style(ax, r"SSA  [m$^{-1}$]  (interface length per cell area)", k_lab)
    ax.legend(fontsize=11, loc="upper right", frameon=False)
    ax.set_title(f"$k_\\mathrm{{iso}}$ vs SSA, {fixed_txt}  (time runs right to left)",
                 fontsize=15)
    save(fig, out / "absolute" / "keff_ssa.png")

    # normalized vs t/tau_sub
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    use_tau = all(c["tau"] for c in conds.values())
    for v in vals:
        c = conds[v]
        x = c["t"] / c["tau"] if use_tau else c["t"] / DAY
        _curve(ax, x, c["kn"], c["knlo"], c["knhi"], c["ib"], col[v], lab(v))
    _style(ax, r"$t\,/\,\tau_\mathrm{sub}$" if use_tau else "Time [d]", r"$k\,/\,k_b$")
    ax.legend(fontsize=11, loc="lower right", frameon=False)
    ax.set_title(f"Normalized $k_\\mathrm{{iso}}$ vs normalized time, {fixed_txt}",
                 fontsize=15)
    save(fig, out / "normalized" / "keff_time.png",
         f"\nb = t = {baseline_tau:.3g} tau_sub for every run; each curve normalized "
         f"by its own seeds' baselines" + ("" if use_tau else "; no tau_sub, time in days"))

    # normalized vs SSA/SSA_b
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    ax.axvline(1.0, color="#999999", lw=0.8, ls=":")
    for v in vals:
        c = conds[v]
        _curve(ax, c["sn"], c["kn"], c["knlo"], c["knhi"], c["ib"], col[v], lab(v))
    _style(ax, r"SSA$\,/\,$SSA$_b$", r"$k\,/\,k_b$")
    ax.legend(fontsize=11, loc="upper right", frameon=False)
    ax.set_title(f"Normalized $k_\\mathrm{{iso}}$ vs normalized SSA, {fixed_txt}",
                 fontsize=15)
    save(fig, out / "normalized" / "keff_ssa.png",
         f"\nb = t = {baseline_tau:.3g} tau_sub for every run")
    return written


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("roots", nargs="+", type=Path)
    p.add_argument("--out", type=Path, default=None,
                   help="output directory (default: <first root>/compare)")
    p.add_argument("--baseline-days", type=float, default=1.0,
                   help="baseline time at the SLOWEST condition (largest tau_sub); "
                        "other conditions use the same t/tau_sub (default 1)")
    p.add_argument("--baseline-tau", type=float, default=None,
                   help="baseline as t/tau_sub directly; overrides --baseline-days")
    a = p.parse_args(argv)

    runs = discover(a.roots)
    if not runs:
        print("no runs with a k_eff CSV and phi/T/seed in the name; nothing to compare")
        return 0
    out = a.out or a.roots[0] / "compare"
    if any(r["tau"] is None for r in runs):
        print("  some runs have no tau_sub in outp.txt; cannot place a common baseline")
        return 1
    base_tau = a.baseline_tau if a.baseline_tau is not None else \
        a.baseline_days * DAY / max(r["tau"] for r in runs)
    print(f"  baseline at t/tau_sub = {base_tau:.3g} for every run")
    by = defaultdict(list)
    for r in runs:
        by[(r["phi"], r["T"])].append(r)
    print(f"  {len(runs)} runs in {len(by)} conditions: " +
          ", ".join(f"phi {k[0]:g} T {k[1]} ({len(v)} seeds)" for k, v in sorted(by.items())))
    conds = {k: condition_mean(v, base_tau) for k, v in by.items()}

    written = []
    for phi in sorted({k[0] for k in conds}):
        g = {T: c for (ph, T), c in conds.items() if ph == phi}
        if len(g) < 2:
            print(f"  skip by_phi/phi{phi:g}: only one temperature")
            continue
        written += draw_group(g, "T", f"φ = {phi:g}", out / "by_phi" / f"phi{phi:g}",
                              base_tau)
    for T in sorted({k[1] for k in conds}):
        g = {ph: c for (ph, TT), c in conds.items() if TT == T}
        if len(g) < 2:
            print(f"  skip by_T/T{T}: only one porosity")
            continue
        written += draw_group(g, "phi", f"T = {T} °C", out / "by_T" / f"T{T}",
                              base_tau)
    for w in written:
        print(f"  wrote {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
