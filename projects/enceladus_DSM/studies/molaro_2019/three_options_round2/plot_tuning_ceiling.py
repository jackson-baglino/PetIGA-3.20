#!/usr/bin/env python3
"""
plot_tuning_ceiling.py — the result that ends the tuning campaign.

Two rounds were run at T = -20 C. Within each D_v family the two rounds differ
ONLY in the wall humidity, so each family is a clean two-point sweep in the one
free boundary parameter, at fixed transport. Plotting neck-at-78-min against
grain-recession-at-78-min turns each family into a LINE, and Molaro's measured
recession (-2.93 %) is a vertical cut through it.

Where those lines cross the cut is the neck the model predicts once it is
required to also reproduce the observed mass loss -- i.e. the only comparison
that is fair, because a model that grows a good neck while losing none of the
grain is not reproducing the experiment.

The answer is that D_v x30 and D_v x100 cross at the SAME neck, 56.9 um, 75 %
of the observed growth. A factor 3.3 in effective transport buys 0.01 um. That
is a ceiling, not a slow approach to the data, and it is why further tuning of
the vapour-only model is not worth the core-hours.

Inputs are the two committed summary.csv files, so this reproduces without the
raw batches (round 1's snapshots have since been deleted).

Usage:  python studies/molaro_2019/three_options_round2/plot_tuning_ceiling.py
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
H_EQ = 1.0000234
WINDOW = 78 * 60.0

# family label -> (round-1 key, round-2 key, colour)
FAMILIES = [
    ("untuned  (nominal $D_v$)", "1 untuned",   "1 untuned h.99715",    "#B23A48"),
    ("$D_v \\times 30$",         "2a D_v x30",  "2a D_v x30 h.99928",   "#2E86AB"),
    ("$D_v \\times 100$",        "2b D_v x100", "2b D_v x100 h.99923",  "#1B4965"),
]


def load(path):
    out = {}
    for r in csv.DictReader(open(path)):
        out[r["label"]] = (float(r["humidity"]),
                           float(r["neck_w_at_78min_um"]),
                           float(r["dR_large_78min_pct"]))
    return out


def main():
    rows = [l.split(",") for l in
            (REPO / "inputs/validation/molaro2019_fig11_T-20.csv").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    m_t = np.array([float(r[0]) for r in rows]) * 60.0
    m_w = np.array([float(r[1]) for r in rows]) * 1e-6
    m_lg = np.array([float(r[4]) for r in rows]) * 1e-6 / 2.0
    c = np.polyfit(m_t, m_lg, 1)
    dR_target = 100.0 * c[0] * WINDOW / np.polyval(c, 0.0)
    w0, w1 = m_w[0] * 1e6, m_w[-1] * 1e6

    r1 = load(HERE.parent / "three_options" / "summary.csv")
    r2 = load(HERE / "summary.csv")

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.2),
                                  gridspec_kw={"width_ratios": [1.35, 1]})
    fig.subplots_adjust(wspace=0.28, left=0.07, right=0.985, top=0.87, bottom=0.13)

    # --- A. the ceiling ----------------------------------------------------
    ax.axvline(dR_target, color="k", ls="--", lw=1.4, zorder=2)
    ax.text(dR_target - 0.06, w0 + 1.0, "Molaro's measured\nrecession, −2.93 %",
            ha="right", va="bottom", fontsize=9)
    ax.axhline(w1, color="k", ls=":", lw=1.4, zorder=2)
    # ABOVE the line, in headroom opened for it: the region below is where the
    # star and the D_v x100 endpoint both sit, and the label collided with both.
    ax.text(0.99, w1 + 0.22, f"Molaro's measured neck, {w1:.2f} µm",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=9)

    crossings = []
    for label, k1, k2, col in FAMILIES:
        h1, n1, d1 = r1[k1]
        h2, n2, d2 = r2[k2]
        ax.plot([d1, d2], [n1, n2], "-o", color=col, lw=2.2, ms=8,
                mec="white", mew=1.4, label=label, zorder=4)
        # annotate each end with its wall
        for d, n, h in ((d1, n1, h1), (d2, n2, h2)):
            ax.annotate(f"1−h={H_EQ-h:.1e}", (d, n), textcoords="offset points",
                        xytext=(6, -12), fontsize=7.5, color=col)
        f = (dR_target - d1) / (d2 - d1)
        nc = n1 + f * (n2 - n1)
        crossings.append((label, nc, col))
        ax.plot([dR_target], [nc], marker="*", ms=17, color=col,
                mec="white", mew=1.2, zorder=6)

    ax.scatter([dR_target], [w1], marker="*", s=340, color="k", zorder=7,
               label="Molaro et al. (2019)")
    ax.set_xlabel("large-grain recession at t*+78 min  [%]")
    ax.set_ylabel("neck width at t*+78 min  [µm]")
    ax.set_title("A. Each $D_v$ is a line; the wall moves you along it",
                 fontsize=11, loc="left")
    ax.legend(fontsize=9, frameon=False, loc="lower left")
    ax.grid(alpha=0.25)
    ax.set_ylim(top=w1 + 2.2)   # headroom for the neck-target label
    ax.invert_xaxis()

    # --- B. what the ceiling is --------------------------------------------
    names = [c[0] for c in crossings] + ["Molaro et al."]
    vals = [c[1] for c in crossings] + [w1]
    cols = [c[2] for c in crossings] + ["k"]
    y = np.arange(len(names))[::-1]
    ax2.barh(y, [v - w0 for v in vals], left=w0, color=cols, height=0.55, zorder=3)
    for yy, v in zip(y, vals):
        share = 100.0 * (v - w0) / (w1 - w0)
        ax2.text(v + 0.4, yy, f"{v:.2f} µm   ({share:.0f} %)",
                 va="center", fontsize=9.5)
    ax2.set_yticks(y)
    ax2.set_yticklabels(names, fontsize=9.5)
    ax2.set_xlim(w0, w1 + 16)   # room for the value labels
    ax2.axvline(w1, color="k", ls=":", lw=1.4, zorder=4)
    ax2.set_xlabel("neck width at t*+78 min, at matched recession  [µm]")
    ax2.set_title("B. $\\times30$ and $\\times100$ land in the same place",
                  fontsize=11, loc="left")
    ax2.grid(alpha=0.25, axis="x")

    fig.suptitle("Transport tuning saturates: once the model must also reproduce the "
                 "measured grain recession,\nthe vapour-only neck stops at ~75 % of the "
                 "observed growth regardless of $D_v$", fontsize=12, y=0.985)
    out = HERE / "tuning_ceiling.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")
    for label, nc, _ in crossings:
        print(f"  {label:26s} → {nc:6.2f} µm   "
              f"({100.0*(nc-w0)/(w1-w0):5.1f} % of observed growth)")
    print(f"  {'Molaro et al. (2019)':26s} → {w1:6.2f} µm   (100.0 %)")


if __name__ == "__main__":
    main()
