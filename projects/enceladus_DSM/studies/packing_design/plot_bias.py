#!/usr/bin/env python3
"""Plot the eps ladder: how much each 'problem' actually biases the answer.

Everything that changes with eps on a FIXED packing is a diffuse-band
artifact; the eps -> 0 intercept is the material. Fits are linear in eps and
the intercept is quoted with the spread between a linear and a quadratic fit,
so the extrapolation carries its own error bar rather than a false precision.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "preprocess"))
import figstyle as fs                                       # noqa: E402

HERE = Path(__file__).parent
rows = list(csv.DictReader((HERE / "bias.csv").open()))
eps = np.array([float(r["eps_m"]) for r in rows])
order = np.argsort(eps)
eps = eps[order]
g = lambda k: np.array([float(r[k]) for r in rows])[order]

PROD = 45e-9            # the eps these packings are built for


def extrap(y):
    """(intercept, uncertainty) from linear vs quadratic fits in eps."""
    lin = np.polyfit(eps, y, 1)
    quad = np.polyfit(eps, y, 2)
    a, b = np.polyval(lin, 0.0), np.polyval(quad, 0.0)
    return 0.5 * (a + b), abs(a - b)


fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.9))
panels = [
    ("k_iso", "$k_{\\rm eff}$  [W m$^{-1}$K$^{-1}$]",
     "(a)  conductivity: biased +33%", fs.C[0], False),
    ("ssa_per_m", "SSA  [m$^{-1}$]", "(b)  SSA: biased -11%", fs.C[2], False),
    ("D_xx_over_D0", "$D_{\\rm eff}/D_0$", "(c)  vapour: geometric, not band",
     fs.C[1], True),
]
summary = []
for ax, (key, ylab, title, col, logy) in zip(axes, panels):
    y = g(key)
    ax.plot(eps * 1e9, y, "o", color=col, ms=7, markeredgecolor="white",
            markeredgewidth=1.0, zorder=3, label="measured")
    if not logy:
        i0, du = extrap(y)
        xf = np.linspace(0, eps.max() * 1.05, 50)
        ax.plot(xf * 1e9, np.polyval(np.polyfit(eps, y, 1), xf), "-",
                color=fs.MUTED, lw=1.3, zorder=2, label=r"fit $\to\varepsilon=0$")
        ax.plot([0], [i0], "*", color=fs.INK, ms=14, zorder=4, clip_on=False,
                label="sharp limit")
        at = float(np.interp(PROD, eps, y))
        bias = (at - i0) / i0
        summary.append((key, at, i0, du, bias))
        ax.annotate(f"{bias:+.0%} at\n$\\varepsilon$ = 45 nm",
                    xy=(PROD * 1e9, at), xytext=(PROD * 1e9 + 8, at),
                    fontsize=fs.FS_NOTE, color=fs.INK, va="center",
                    arrowprops=dict(arrowstyle="->", color=fs.INK, lw=1.0))
    else:
        ax.plot(eps * 1e9, g("D_yy_over_D0"), "s", color=fs.C[5], ms=7,
                markeredgecolor="white", markeredgewidth=1.0, zorder=3,
                label="$D_{yy}$ (along deposition)")
        ax.set_ylim(5e-4, 1e-2)
        ax.axhline(0.325 / 2.0, color=fs.MUTED, ls="--", lw=1.2)
        ax.text(eps.max() * 1e9, 0.325 / 2.0 * 1.12,
                "a well-connected pore\nat this porosity", ha="right",
                va="bottom", fontsize=fs.FS_NOTE, color=fs.MUTED)
        ax.set_ylim(5e-4, 0.35)
    ax.axvline(PROD * 1e9, color=fs.MUTED, ls=":", lw=1.0, zorder=1)
    fs.style(ax, r"$\varepsilon$  [nm]", ylab, title, logy=logy)
    ax.set_xlim(left=0)
    ax.legend(fontsize=fs.FS_LEG, frameon=False, loc="best")

fig.tight_layout()
fs.save(fig, HERE, "bias", dpi=190)

print(f"{'quantity':>14} {'at 45nm':>10} {'eps->0':>10} {'+/-':>9} {'bias':>8}")
for k, at, i0, du, b in summary:
    print(f"{k:>14} {at:10.4g} {i0:10.4g} {du:9.2g} {b:+8.1%}")
print(f"\nD_eff/D_0 at 45 nm:  x {g('D_xx_over_D0')[np.argmin(abs(eps-PROD))]:.2e}"
      f"   y {g('D_yy_over_D0')[np.argmin(abs(eps-PROD))]:.2e}")
