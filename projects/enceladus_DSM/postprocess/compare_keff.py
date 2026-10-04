#!/usr/bin/env python3
"""compare_keff.py — overlay k_eff curves across temperatures and porosities.

    python3 compare_keff.py <root> [<root> ...] [--out <dir>]

Finds every run directory under the roots that holds a k_eff CSV, reads its
porosity, temperature and seed from the directory name
(``..._phi0.325_..._seed3_..._T-20__...``), and writes two families of figures:

    by_phi/phi<X>/   fixed porosity, one curve per temperature
    by_T/T<Y>/       fixed temperature, one curve per porosity

each with absolute/ and normalized/ versions of

    keff_time.png    k_iso against time        (normalized: k/k_0 vs t/tau_sub)
    keff_ssa.png     k_iso against SSA         (normalized: k/k_0 vs SSA/SSA_0)

A group is drawn only when it holds at least two values of the varied
parameter; the script says which groups it skipped and why.

SEEDS. Each curve is the mean over the seeds of that condition, and the band
around it is the seed min..max, so the spread between conditions can be read
against the spread within one. The mean is taken at common times (time plots),
and against SSA the mean k is plotted at the mean SSA of those same times --
both coordinates averaged at matched t, so no seed's SSA axis is resampled.
The common times run only over the span EVERY seed covers, from the latest
opening sample to the earliest last sample: seeds are interpolated between
their own samples, never extrapolated past either end.

THE OPENING SAMPLE IS t = 0 and the normalization, as in plot_keff.py and
plot_keff_snapshots.py: each run opens on its first k_eff sample with
1 s <= t <= 1 h, and is divided by its own measured values there (k_0,
SSA_0). Everything after it is drawn, the fast early relaxation included.
(Until 2026-09-30 every run was divided by values interpolated to
t = 11 tau_sub, labelled "b", and the earlier part was drawn faded.)

COLOUR. Ordered parameters, so sequential maps sampled light -> dark with
the value: temperature on cmocean `thermal` (0.15-0.85, so neither the
near-black nor the pale-yellow end is used), porosity on an amp map that runs
to black instead of white. Single-run figures use a cool map instead
(plot_keff.py), so neither reads as the other.
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
from plot_keff import DAY, from_opening, load, read_tau_sub    # noqa: E402
import cmocean                                                  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap           # noqa: E402

# Sequential, one map per varied parameter; SPAN is the stretch of it used.
_AMP = cmocean.cm.amp
CMAP = {"T": (cmocean.cm.thermal, (0.15, 0.85)),
        "phi": (LinearSegmentedColormap.from_list(
            "amp_black", ["#0b0b0b", _AMP(0.85), _AMP(0.60), _AMP(0.38)]), (0.0, 1.0))}
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
            data = from_opening(data)
            data.update(dir=d, phi=meta[0], T=meta[1], seed=meta[2],
                        tau=read_tau_sub(d))
            runs.append(data)
    return runs


def condition_mean(runs):
    """Seed-mean curves for one condition on a common log-spaced time grid,
    over the span every seed covers. Each seed is normalized by its own
    opening sample before averaging."""
    taus = {r["tau"] for r in runs}
    tau = taus.pop() if len(taus) == 1 else None
    t_start = max(r["t"][0] for r in runs)
    t_end = min(r["t"][-1] for r in runs)
    tg = np.geomspace(t_start, t_end, N_GRID)
    K = np.array([np.interp(tg, r["t"], r["kiso"]) for r in runs])
    S = np.array([np.interp(tg, r["t"], r["ssa"]) for r in runs])
    Kn = np.array([np.interp(tg, r["t"], r["kiso"] / r["kiso"][0]) for r in runs])
    Sn = np.array([np.interp(tg, r["t"], r["ssa"] / r["ssa"][0]) for r in runs])
    return dict(t=tg, tau=tau, n=len(runs),
                k=K.mean(0), klo=K.min(0), khi=K.max(0), s=S.mean(0),
                kn=Kn.mean(0), knlo=Kn.min(0), knhi=Kn.max(0), sn=Sn.mean(0))


def _style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def _curve(ax, x, y, lo, hi, color, label, band=True):
    """The seed mean, solid, over the span every seed covers, with its band."""
    ax.plot(x, y, "-", color=color, lw=2.2, label=label)
    if band:
        ax.fill_between(x, lo, hi, color=color, alpha=0.15, lw=0)


def draw_group(conds, vary, fixed_txt, out: Path, extra_note=""):
    """conds: {value: condition_mean dict}. vary: 'T' or 'phi'."""
    vals = sorted(conds)
    cmap, (a, b) = CMAP[vary]
    col = {v: cmap(a + (b - a) * i / max(1, len(vals) - 1)) for i, v in enumerate(vals)}
    lab = (lambda v: f"T = {v} °C") if vary == "T" else (lambda v: f"φ = {v:g}")
    nseed = sorted({c["n"] for c in conds.values()})
    note = (f"line = mean over {'/'.join(map(str, nseed))} seeds, band = seed min–max; "
            f"t = 0 and subscript 0 = each run's opening sample (first with t >= 1 s)"
            + extra_note)
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
        _curve(ax, c["t"] / DAY, c["k"], c["klo"], c["khi"], col[v], lab(v))
    _style(ax, "Time [d]", k_lab)
    ax.legend(fontsize=11, loc="lower right", frameon=False)
    ax.set_title(f"$k_\\mathrm{{iso}}$ vs time, {fixed_txt}", fontsize=15)
    save(fig, out / "absolute" / "keff_time.png")

    # absolute vs SSA
    fig, ax = plt.subplots(figsize=(10, 6))
    for v in vals:
        c = conds[v]
        _curve(ax, c["s"], c["k"], c["klo"], c["khi"], col[v], lab(v))
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
        _curve(ax, x, c["kn"], c["knlo"], c["knhi"], col[v], lab(v))
    _style(ax, r"$t\,/\,\tau_\mathrm{sub}$" if use_tau else "Time [d]",
           r"$k_\mathrm{iso}\,/\,k_{\mathrm{iso},0}$")
    ax.legend(fontsize=11, loc="lower right", frameon=False)
    ax.set_title(f"Normalized $k_\\mathrm{{iso}}$ vs normalized time, {fixed_txt}",
                 fontsize=15)
    save(fig, out / "normalized" / "keff_time.png",
         "" if use_tau else "; no tau_sub, time in days")

    # normalized vs PHYSICAL time: the same curves as the t/tau_sub figure,
    # but on the clock -- how much faster one condition gets there
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    for v in vals:
        c = conds[v]
        _curve(ax, c["t"] / DAY, c["kn"], c["knlo"], c["knhi"], col[v], lab(v))
    _style(ax, "Time [d]", r"$k_\mathrm{iso}\,/\,k_{\mathrm{iso},0}$")
    ax.legend(fontsize=11, loc="lower right", frameon=False)
    ax.set_title(f"Normalized $k_\\mathrm{{iso}}$ vs time, {fixed_txt}", fontsize=15)
    save(fig, out / "normalized" / "keff_time_days.png")

    # normalized vs SSA/SSA_0
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(1.0, color="#999999", lw=0.8, ls=":")
    ax.axvline(1.0, color="#999999", lw=0.8, ls=":")
    for v in vals:
        c = conds[v]
        _curve(ax, c["sn"], c["kn"], c["knlo"], c["knhi"], col[v], lab(v))
    _style(ax, r"SSA$\,/\,$SSA$_0$", r"$k_\mathrm{iso}\,/\,k_{\mathrm{iso},0}$")
    ax.legend(fontsize=11, loc="upper right", frameon=False)
    ax.set_title(f"Normalized $k_\\mathrm{{iso}}$ vs normalized SSA, {fixed_txt}",
                 fontsize=15)
    save(fig, out / "normalized" / "keff_ssa.png")
    return written


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("roots", nargs="+", type=Path)
    p.add_argument("--out", type=Path, default=None,
                   help="output directory (default: <first root>/compare)")
    a = p.parse_args(argv)

    runs = discover(a.roots)
    if not runs:
        print("no runs with a k_eff CSV and phi/T/seed in the name; nothing to compare")
        return 0
    out = a.out or a.roots[0] / "compare"
    by = defaultdict(list)
    for r in runs:
        by[(r["phi"], r["T"])].append(r)
    print(f"  {len(runs)} runs in {len(by)} conditions: " +
          ", ".join(f"phi {k[0]:g} T {k[1]} ({len(v)} seeds)" for k, v in sorted(by.items())))
    conds = {k: condition_mean(v) for k, v in by.items()}

    written = []
    for phi in sorted({k[0] for k in conds}):
        Ts = sorted(T for (ph, T) in by if ph == phi)
        if len(Ts) < 2:
            print(f"  skip by_phi/phi{phi:g}: only one temperature")
            continue
        # PAIRED: temperatures are compared on the packings they share. A
        # 5-seed mean at one T against a 1-seed curve at another mixes the
        # temperature effect with packing-to-packing scatter (~9% in k at
        # phi 0.325), and the T effect here is a time rescaling per packing.
        common = set.intersection(*({r["seed"] for r in by[(phi, T)]} for T in Ts))
        if common:
            g = {T: condition_mean([r for r in by[(phi, T)] if r["seed"] in common]) for T in Ts}
            note = f"; PAIRED: seeds {', '.join(map(str, sorted(common)))} (common to all T)"
        else:
            print(f"  WARNING by_phi/phi{phi:g}: no seed common to all T; unpaired means")
            g = {T: conds[(phi, T)] for T in Ts}
            note = "; UNPAIRED (no seed common to all T)"
        written += draw_group(g, "T", f"φ = {phi:g}", out / "by_phi" / f"phi{phi:g}", note)
    for T in sorted({k[1] for k in conds}):
        g = {ph: c for (ph, TT), c in conds.items() if TT == T}
        if len(g) < 2:
            print(f"  skip by_T/T{T}: only one porosity")
            continue
        written += draw_group(g, "phi", f"T = {T} °C", out / "by_T" / f"T{T}")
    for w in written:
        print(f"  wrote {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
