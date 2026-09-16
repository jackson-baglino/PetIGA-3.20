#!/usr/bin/env python3
"""Bias vs domain size: is +22% a bulk property or a small-sample artifact?"""
from __future__ import annotations
import csv, sys
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[1] / "preprocess"))
import figstyle as fs                                       # noqa: E402

rows = list(csv.DictReader((HERE / "rev_bias.csv").open()))
lr = np.array([int(r["L_over_R"]) for r in rows])
bias = np.array([float(r["bias"]) for r in rows])
k45 = np.array([float(r["k_at_eps1"]) for r in rows])
an = np.array([float(r["k_anis"]) for r in rows])
sizes = sorted(set(lr))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.1))

a1.plot(lr, bias * 100, "o", color=fs.C[1], ms=7, alpha=0.55,
        markeredgecolor="white", markeredgewidth=1.0, label="individual seeds")
m = [bias[lr == s].mean() * 100 for s in sizes]
a1.plot(sizes, m, "-o", color=fs.C[0], ms=9, lw=2, label="mean of 3 seeds")
a1.set_xscale("log"); a1.set_xticks(sizes); a1.set_xticklabels(sizes)
fs.style(a1, r"domain size  $L/R_{\rm ave}$", r"$k_{\rm eff}$ bias at $\varepsilon$=45 nm  [%]",
         "(a)  the bias does NOT shrink with domain", logy=False)
a1.set_ylim(0, 36)
a1.legend(fontsize=fs.FS_LEG, frameon=False, loc="lower left")
a1.text(0.97, 0.95, "flat -> a bulk property of the\nmicrostructure, not a\nsmall-sample artifact",
        transform=a1.transAxes, ha="right", va="top", fontsize=fs.FS_NOTE, color=fs.INK)

cv = [k45[lr == s].std(ddof=1) / k45[lr == s].mean() * 100 for s in sizes]
cva = [an[lr == s].std(ddof=1) * 100 for s in sizes]
a2.plot(sizes, cv, "-o", color=fs.C[0], ms=9, lw=2,
        label=r"seed scatter in $k_{\rm eff}$")
a2.plot(sizes, cva, "-s", color=fs.C[2], ms=8, lw=2,
        label=r"seed scatter in anisotropy")
a2.axhline(5.0, color=fs.MUTED, ls="--", lw=1.2)
a2.text(sizes[0], 5.4, "5% — a usable REV", fontsize=fs.FS_NOTE, color=fs.MUTED)
a2.set_xscale("log"); a2.set_xticks(sizes); a2.set_xticklabels(sizes)
fs.style(a2, r"domain size  $L/R_{\rm ave}$", "seed-to-seed scatter  [%]",
         "(b)  but the SCATTER does", logy=False)
a2.legend(fontsize=fs.FS_LEG, frameon=False, loc="upper right")
fig.tight_layout()
fs.save(fig, HERE, "rev_bias", dpi=190)
