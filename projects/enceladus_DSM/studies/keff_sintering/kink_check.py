#!/usr/bin/env python3
"""Is k_eff sampled finely enough that the plotted curves have no corners?

    venv_enceladus/bin/python studies/keff_sintering/kink_check.py <campaign dir> [--limit 0.3]

A curve drawn through samples shows a corner where three consecutive samples
are not nearly collinear. For every interior sample i of every run:

    corner_i = | k_i - (line through samples i-1 and i+1 evaluated at x_i) |
               / (k_max - k_min)                       [% of the plotted range]

measured on both axes the figures use: x = time [d] and x = SSA. It is what
you would see, and it is also 4x the linear-interpolation error between
neighbouring samples (so the curve's own error is corner/4).

Reported per run: the worst corner on each axis and where it is; then the worst
per temperature. A corner above --limit (default 0.3% of the range, ~3 px on a
1000 px tall figure) is flagged. Writes kink_check.csv and kink_check.png
into <campaign>/compare/kinks/.

A REAL kink (a pore pinching off) also shows up here. It is told apart from a
sampling corner by its neighbours: sampling corners come in runs of similar
size wherever the curve bends; a real event is one isolated spike.
"""
from __future__ import annotations

import argparse, csv, re, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "postprocess"))
from plot_keff import load, read_tau_sub  # noqa: E402

PAT = re.compile(r"phi([\d.]+)_Rave50um_LR(\d+)_seed(\d+)_L[\d.]+mm_(eps[^_]+)_perxy_T(-?\d+)__")
DAY = 86400.0


def corners(x, k):
    x0, x1, x2 = x[:-2], x[1:-1], x[2:]
    lin = k[:-2] + (k[2:] - k[:-2]) * (x1 - x0) / np.where(x2 != x0, x2 - x0, 1.0)
    return np.abs(k[1:-1] - lin) / (k.max() - k.min()) * 100


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--limit", type=float, default=0.3)
    a = ap.parse_args()
    rows = []
    for kf in sorted(a.root.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m:
            continue
        r = load(kf.parent)
        if r is None or len(r["t"]) < 5:
            continue
        i0 = int(np.argmax(r["t"] >= 1.0))
        t, s, k = r["t"][i0:] / DAY, r["ssa"][i0:], r["kiso"][i0:]
        ct, cs = corners(t, k), corners(s, k)
        it, is_ = int(np.argmax(ct)), int(np.argmax(cs))
        # isolated spike (a real event) vs a run of similar corners (sampling)
        nb = np.r_[ct[max(it - 2, 0):it], ct[it + 1:it + 3]]
        rows.append(dict(phi=float(m.group(1)), LR=int(m.group(2)), seed=int(m.group(3)),
                         eps=m.group(4), T=int(m.group(5)), n=len(k),
                         corner_t_pct=float(ct[it]), at_t_d=float(t[it + 1]),
                         corner_ssa_pct=float(cs[is_]), at_ssa_t_d=float(t[is_ + 1]),
                         p95_t_pct=float(np.percentile(ct, 95)),
                         isolated=bool(len(nb) and ct[it] > 4 * np.median(nb)),
                         max_gap_d=float(np.diff(t).max())))
    out = a.root / "compare" / "kinks"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "kink_check.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    G = defaultdict(list)
    for x in rows:
        G[x["T"]].append(x)
    print(f"{len(rows)} runs; corner = % of the plotted k range; limit {a.limit}%")
    print(f"{'T':>4} {'runs':>4} {'samples':>8} {'worst (time)':>13} {'median worst':>13} "
          f"{'worst (SSA)':>12} {'max gap d':>10} {'over limit':>10}")
    for T in sorted(G):
        g = G[T]
        print(f"{T:4d} {len(g):4d} {int(np.median([x['n'] for x in g])):8d} "
              f"{max(x['corner_t_pct'] for x in g):12.2f}% {np.median([x['corner_t_pct'] for x in g]):12.2f}% "
              f"{max(x['corner_ssa_pct'] for x in g):11.2f}% {max(x['max_gap_d'] for x in g):10.2f} "
              f"{sum(x['corner_t_pct'] > a.limit for x in g):10d}")
    bad = sorted((x for x in rows if x["corner_t_pct"] > a.limit), key=lambda x: -x["corner_t_pct"])
    for x in bad[:25]:
        print(f"  phi {x['phi']} LR{x['LR']} seed {x['seed']} {x['eps']} T {x['T']:>3}: corner {x['corner_t_pct']:.2f}% "
              f"at t = {x['at_t_d']:.3f} d ({'isolated' if x['isolated'] else 'run of corners'}); n = {x['n']}")

    Ts = sorted(G)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    for j, (key, lab) in enumerate((("corner_t_pct", "worst corner, k vs time"),
                                    ("corner_ssa_pct", "worst corner, k vs SSA"))):
        for i, T in enumerate(Ts):
            v = [x[key] for x in G[T]]
            ax[j].scatter(np.full(len(v), i) + np.linspace(-0.25, 0.25, len(v)), v, s=14,
                          color="#0072B2", alpha=0.7)
        ax[j].axhline(a.limit, color="#D55E00", lw=1, ls="--", label=f"limit {a.limit}%")
        ax[j].set_xticks(range(len(Ts))); ax[j].set_xticklabels([f"{T} °C" for T in Ts])
        ax[j].set(ylabel="% of plotted range", title=f"({'ab'[j]}) {lab}"); ax[j].set_yscale("log")
        ax[j].grid(True, alpha=0.25); ax[j].legend(frameon=False)
        for sp in ("top", "right"):
            ax[j].spines[sp].set_visible(False)
    fig.suptitle("Largest corner in each run's k_eff curve (one point per run)")
    fig.savefig(out / "kink_check.png", dpi=160)
    print(f"wrote {out}/kink_check.png and .csv")


if __name__ == "__main__":
    main()
