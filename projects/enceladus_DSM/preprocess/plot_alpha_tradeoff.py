#!/usr/bin/env python3
"""The alpha_c decision figure: mesh, asymptotic validity and regime vs alpha_c.

WHY A SECOND FIGURE
-------------------
plot_alpha_kinetics.py sweeps TEMPERATURE with alpha_c(T) supplied by a model.
This one sweeps alpha_c ITSELF at fixed T, because alpha_c is the choice, and
after K&P Eq. 42 was enforced per channel (comp_eps, 2026-08-11) the three
things it controls no longer move together:

  A  mesh cost         Eq. 42a ceiling scales as beta_HK ~ 1/alpha_c, so a
                       large alpha_c forces an intractable grid
  B  asymptotic validity  beta_eff/beta = 1 - a1*a2*W/(D*beta) per channel --
                       the quantitative form of Eq. 42's "W << D*beta0". It
                       goes NEGATIVE (wrong-sign effective kinetics) at large
                       alpha_c
  C  transport regime  r_neck/L* with L* = beta_HK*D_v. Kuczynski's
                       r ~ t^(1/3) is derived for the attachment-limited case
                       (r/L* < 1), so this decides whether 1/3 is even the
                       right thing to expect

All three improve as alpha_c falls, and all three are plotted against the
range Braun et al. (2024) actually report.

WHAT BRAUN ACTUALLY SAYS
------------------------
Braun, Fourteau & Lowe, "A rigorous approach to the specific surface area
evolution in snow during temperature gradient metamorphism", The Cryosphere
18(4) 1653-1668 (2024), doi:10.5194/tc-18-1653-2024:

  "the best agreement ... is obtained for values in the range 1e-3 < alpha < 1e-1"
  optimal fits: 1e-1 and 1e-1.5 (series 1), 1e-2.25 (series 2)
  at mean T = -8.1 / -7.6 C under gradients of 47 / 55 K/m

Two things to keep straight when using it:
  * it is a RANGE from fitting SSA evolution, not a measurement of alpha, and
    they call it an "effective" coefficient absorbing microscale variation;
  * it is NOT a function of temperature. They state explicitly that alpha
    depends on temperature, supersaturation and crystallographic orientation,
    and then hold it constant anyway. Any alpha_c(T) we impose is our
    assumption, not theirs.

Usage:
    python preprocess/plot_alpha_tradeoff.py --out <dir>
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comp_eps as ce

C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

BRAUN_LO, BRAUN_HI = 1.0e-3, 1.0e-1          # their stated best-agreement range
BRAUN_FITS = [(1.0e-1, "series 1"), (10 ** -1.5, "series 1 late"),
              (10 ** -2.25, "series 2")]


def _ax(a, xlabel, ylabel, title, logy=True):
    a.set_xscale("log")
    if logy:
        a.set_yscale("log")
    a.set_xlabel(xlabel, fontsize=9)
    a.set_ylabel(ylabel, fontsize=9)
    a.set_title(title, fontsize=9.5, color=INK, loc="left")
    a.grid(alpha=0.25, lw=0.5, color=GRID, which="both")
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=8)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--Lx", type=float, default=2.195942e-3)
    p.add_argument("--Ly", type=float, default=6.0e-4)
    p.add_argument("--Rave", type=float, default=5.0e-4)
    p.add_argument("--rneck", type=float, default=4.5e-5)
    p.add_argument("--T0", type=float, default=-3.0)
    p.add_argument("--floor_rR", type=float, default=0.09,
                   help="target neck-resolution floor r/R = sqrt(12 eps/R)")
    p.add_argument("--safety", type=float, default=0.5)
    p.add_argument("--label", default="Demmenie 1 mm pair, T = −3 °C")
    p.add_argument("--out", type=Path, default=Path("."))
    args = p.parse_args()

    T = args.T0
    Di, Da = ce._K_I / ce._C_I, ce._K_A / ce._C_A
    Dv = ce.Dv_T(T)
    a1, a2 = ce._A1, ce._A2
    eps_floor = args.Rave * args.floor_rR ** 2 / 12.0
    n_floor = (math.ceil(math.sqrt(2) * args.Lx / eps_floor)
               * math.ceil(math.sqrt(2) * args.Ly / eps_floor))

    alpha = np.geomspace(1e-5, 3e-1, 220)
    nodes, ratios, LL = [], {"ice": [], "air": [], "vapour": []}, []
    for a in alpha:
        b = ce.beta_HK(T, a)
        eps = args.safety * min(Di * b, Da * b, Dv * b, args.Rave)
        nodes.append(math.ceil(math.sqrt(2) * args.Lx / eps)
                     * math.ceil(math.sqrt(2) * args.Ly / eps))
        # Asymptotic validity is judged at the eps we WILL run (the neck floor),
        # not at the Eq.42 ceiling -- that is the mesh that actually gets used.
        for k, D in (("ice", Di), ("air", Da), ("vapour", Dv)):
            ratios[k].append(1.0 - a1 * a2 * eps_floor / (D * b))
        LL.append(args.rneck / (b * Dv))
    nodes = np.array(nodes, float)
    LL = np.array(LL)

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.7))
    fig.suptitle(f"Choosing $\\alpha_c$: mesh, asymptotic validity and transport "
                 f"regime — {args.label}", fontsize=11.5, color=INK)

    def braun(a):
        a.axvspan(BRAUN_LO, BRAUN_HI, color=C[0], alpha=0.10, zorder=0)
        for v, _ in BRAUN_FITS:
            a.axvline(v, color=C[0], lw=1.0, ls=(0, (2, 2)), alpha=0.65, zorder=1)

    # A -- mesh -------------------------------------------------------------
    a = ax[0]; braun(a)
    a.plot(alpha, nodes / 1e6, lw=2.4, color=C[3], label="Eq. 42 ceiling")
    a.axhline(n_floor / 1e6, lw=2.0, ls=(0, (4, 2)), color=C[1],
              label=f"neck floor $r/R$={args.floor_rR:g} ({n_floor/1e6:.0f}M)")
    a.axhline(50.0, lw=1.4, ls=":", color=MUTED)
    a.annotate("intractable (>50M)", (1.1e-5, 62), fontsize=7.5, color=MUTED)
    a.annotate("Braun 2024\nrange", (math.sqrt(BRAUN_LO * BRAUN_HI), 2e3),
               fontsize=7.5, color=C[0], ha="center")
    _ax(a, "$\\alpha_c$  [-]", "mesh nodes  [millions]",
        "A.  Mesh forced by Eq. 42a  (∝ $1/\\alpha_c$)")
    a.set_ylim(nodes.min() / 1e6 / 3, 3e4)
    a.legend(fontsize=7.5, frameon=False, loc="upper left")

    # B -- asymptotic validity ----------------------------------------------
    a = ax[1]; braun(a)
    for i, (k, col) in enumerate((("ice", C[1]), ("air", C[2]), ("vapour", C[5]))):
        a.plot(alpha, ratios[k], lw=2.2, color=col, label=f"{k} channel")
    a.axhline(0.9, lw=1.4, ls=(0, (4, 2)), color=MUTED)
    a.annotate("K&P target > 0.9", (1.1e-5, 0.915), fontsize=7.5, color=MUTED)
    a.axhline(0.0, lw=1.4, color=INK, alpha=0.6)
    a.annotate("below 0 → $\\beta_{eff}$ has the WRONG SIGN", (1.1e-5, -0.85),
               fontsize=7.5, color=C[1])
    _ax(a, "$\\alpha_c$  [-]", "$\\beta_{eff}/\\beta$", logy=False,
        title="B.  Eq. 42 validity at the mesh we run ($\\epsilon$ = neck floor)")
    a.set_ylim(-1.2, 1.15)
    a.legend(fontsize=7.5, frameon=False, loc="lower right")

    # C -- transport regime ---------------------------------------------------
    a = ax[2]; braun(a)
    a.plot(alpha, LL, lw=2.4, color=C[2])
    a.axhline(1.0, lw=1.6, ls=(0, (4, 2)), color=INK, alpha=0.7)
    a.annotate("attachment-limited  ($t^{1/3}$ is derived here)", (1.1e-5, 0.25),
               fontsize=7.5, color=C[2])
    a.annotate("diffusion-limited", (1.1e-5, 12), fontsize=7.5, color=C[1])
    _ax(a, "$\\alpha_c$  [-]", "$r_{neck}/L^*$,   $L^*=\\beta_{HK}D_v$",
        "C.  Which transport regime the run is in")
    a.set_ylim(LL.min() / 3, LL.max() * 3)

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    args.out.mkdir(parents=True, exist_ok=True)
    png = args.out / "alpha_tradeoff.png"
    fig.savefig(png, dpi=150)
    print(f"plot -> {png}")

    print(f"\n  Braun 2024 range {BRAUN_LO:g}–{BRAUN_HI:g}; their fits "
          f"{', '.join(f'{v:.1e}' for v, _ in BRAUN_FITS)}")
    print(f"  {'alpha_c':>10} {'Mnodes':>9} {'beta_eff/beta(ice)':>20} "
          f"{'r/L*':>8}  regime")
    for a_ in (1e-1, 10 ** -1.5, 1e-2, 10 ** -2.25, 1e-3, 1e-4):
        b = ce.beta_HK(T, a_)
        eps = args.safety * min(Di * b, Da * b, Dv * b, args.Rave)
        n = (math.ceil(math.sqrt(2) * args.Lx / eps)
             * math.ceil(math.sqrt(2) * args.Ly / eps)) / 1e6
        r = 1.0 - a1 * a2 * eps_floor / (Di * b)
        L = args.rneck / (b * Dv)
        reg = "diffusion" if L > 3 else ("mixed" if L > 0.3 else "attachment")
        print(f"  {a_:>10.2e} {n:>9.1f} {r:>20.3f} {L:>8.2f}  {reg}")


if __name__ == "__main__":
    main()
