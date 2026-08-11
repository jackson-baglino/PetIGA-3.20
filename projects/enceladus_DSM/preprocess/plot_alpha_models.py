#!/usr/bin/env python3
"""Choose the alpha_c(T, sigma) law: the candidates and what each one costs.

WHY THIS FIGURE
---------------
alpha_c is the model's one free parameter and it should be DETERMINED by the
local state (T and rho_v), not chosen per run. Two candidate laws are on the
table, and the choice has to be made against three things at once: the
literature band, the validity of the thin-interface asymptotics, and the mesh.
This plots all three so the decision is made on evidence.

THE CANDIDATES
--------------
  ARRH    alpha_c(T) = A*exp(-Ea/RT), anchored on Braun's range endpoints.
          Ea = 63.7 kJ/mol. Note this reproduces the repo's long-standing
          values exactly (9.007e-2 at -3 C, 1.341e-2 at -20 C). NO sigma
          dependence -- it is drawn as one curve on every panel.

  LIBB2   alpha_c(T, sigma) = A*exp(-f*sigma0(T)/sigma), the Libbrecht form
          with BOTH free parameters fitted rather than A hardwired to 1.
          This is the only candidate that responds to rho_v.

WHY THE BAND IS NOT A FIT TARGET
--------------------------------
Braun, Fourteau & Lowe (The Cryosphere 18(4) 1653-1668, 2024) report best
agreement for 1e-3 < alpha < 1e-1 -- but all three of their values come from
essentially ONE temperature (-8.1 / -7.6 C), and series 1 drifts 1e-1 -> 3.2e-2
within a single sample over 160 h as the snow coarsens. Falling curvature means
falling local sigma, and under exp(-sigma0/sigma) that means falling alpha. The
spread is the SUPERSATURATION dependence, not a temperature dependence.

So the band is a consistency check on the resulting alpha(T, sigma) SURFACE,
and an Arrhenius fitted to its endpoints charges to temperature a variation
that is mostly sigma. That is the central argument for preferring LIBB2, and
panel A is where it either shows or does not.

THE SIGMA VALUES PLOTTED
------------------------
Sintering is capillarity-driven, so the local supersaturation is Gibbs-Thomson:
sigma = d0/rho_fillet with rho_fillet ~ r^2/(2R). For both geometries in play
that lands at sigma ~ 4-5e-4:

    Molaro   R=86.75 um, r=20 um -> sigma = 4.4e-4
    Demmenie R=500 um,   r=45 um -> sigma = 4.7e-4

Libbrecht's and Braun's measurements sit at snow-crystal-growth supersaturations
of 1e-2 to 1e-1, two to three decades higher. That gap is the whole difficulty,
and the panels span it deliberately.

sigma = 0 (exact saturation) is NOT plotted: the nucleation form gives alpha = 0
there, which carries no information until we decide what sigma a saturated neck
really sees.

Usage:
    python preprocess/plot_alpha_models.py --out studies/.../results/kinetics
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

# Braun et al. (2024) optimal fits, all at T ~ -8 C.
BRAUN_T = -8.0
BRAUN_FITS = [(1.0e-1, "s1 early"), (10 ** -1.5, "s1 late"), (10 ** -2.25, "s2")]


def _ax(a, xlabel, ylabel, title, logy=True):
    a.set_xlabel(xlabel, fontsize=9)
    a.set_ylabel(ylabel, fontsize=9)
    a.set_title(title, fontsize=9.5, color=INK, loc="left")
    if logy:
        a.set_yscale("log")
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
    p.add_argument("--floor_rR", type=float, default=0.09)
    p.add_argument("--safety", type=float, default=0.5)
    p.add_argument("--label", default="Demmenie 1 mm pair")
    p.add_argument("--out", type=Path, default=Path("."))
    args = p.parse_args()

    T = np.linspace(-40.0, -1.0, 200)
    eps_floor = args.Rave * args.floor_rR ** 2 / 12.0
    n_floor = (math.ceil(math.sqrt(2) * args.Lx / eps_floor)
               * math.ceil(math.sqrt(2) * args.Ly / eps_floor))

    # sigma cases: the sintering value, plus one decade either side to bracket
    # undersaturated / supersaturated conditions around it.
    sig_cases = [(1.0e-4, "undersaturated  σ = 1e-4"),
                 (4.5e-4, "sintering σ = 4.5e-4  (d₀/ρ_fillet)"),
                 (5.0e-3, "supersaturated  σ = 5e-3")]

    a_arr = np.array([ce.alpha_arrhenius(t, clamp=None) for t in T])
    a1c, a2c = ce._A1, ce._A2
    Dm, Dv_of = ce.D_heat_of("mean"), ce.Dv_T

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    fig.suptitle(f"Choosing $\\alpha_c(T,\\sigma)$ — {args.label}; "
                 f"band = Braun et al. (2024) 1e-3–1e-1",
                 fontsize=11.5, color=INK)

    # ---- A: the candidate laws -------------------------------------------
    a = axes[0]
    a.axhspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10,
              label="literature band")
    for i, (s, lbl) in enumerate(sig_cases):
        curve = np.array([ce.alpha_libbrecht2(t, s, clamp=None) for t in T])
        a.plot(T, np.maximum(curve, 1e-40), lw=2.2, color=C[2 + i],
               label=f"LIBB2, {lbl}")
    a.plot(T, a_arr, lw=2.6, color=C[1], ls=(0, (5, 2)),
           label="ARRH $\\alpha_c(T)$ (no σ dependence)")
    for v, lb in BRAUN_FITS:
        a.plot(BRAUN_T, v, "o", ms=7, color=C[0], markerfacecolor="white",
               markeredgewidth=1.8, zorder=6)
    a.annotate("Braun's 3 fits —\nall at one T", (BRAUN_T, 1.3e-1),
               fontsize=7.5, color=C[0], ha="center", va="bottom")
    _ax(a, "T [°C]", "$\\alpha_c$  [-]", "A.  Candidate laws vs the band")
    a.set_ylim(1e-8, 1e0)
    a.legend(fontsize=7, frameon=False, loc="lower right")

    # ---- B: does the asymptotic expansion hold? --------------------------
    a = axes[1]
    a.axhline(0.9, lw=1.4, ls=(0, (4, 2)), color=MUTED)
    a.annotate("K&P target > 0.9", (-39, 0.915), fontsize=7.5, color=MUTED)
    a.axhline(0.0, lw=1.4, color=INK, alpha=0.6)
    for i, (s, lbl) in enumerate(sig_cases):
        r = [1 - a1c * a2c * eps_floor
             / (Dm * ce.beta_HK(t, ce.alpha_libbrecht2(t, s))) for t in T]
        a.plot(T, r, lw=2.2, color=C[2 + i], label=f"LIBB2, σ={s:g}")
    a.plot(T, [1 - a1c * a2c * eps_floor / (Dm * ce.beta_HK(t, al))
               for t, al in zip(T, np.clip(a_arr, ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI))],
           lw=2.6, ls=(0, (5, 2)), color=C[1], label="ARRH")
    _ax(a, "T [°C]", "$\\beta'/\\beta_0'$", logy=False,
        title="B.  Thin-interface validity at the neck-floor $\\epsilon$ (mean channel)")
    a.set_ylim(-0.6, 1.1)
    a.legend(fontsize=7.5, frameon=False, loc="lower right")

    # ---- C: mesh, and which regime -----------------------------------------
    a = axes[2]
    for i, (s, lbl) in enumerate(sig_cases):
        n = []
        for t in T:
            al = ce.alpha_libbrecht2(t, s)
            e = args.safety * min(Dm * ce.beta_HK(t, al), Dv_of(t) * ce.beta_HK(t, al))
            n.append(math.ceil(math.sqrt(2) * args.Lx / e)
                     * math.ceil(math.sqrt(2) * args.Ly / e))
        a.plot(T, np.array(n) / 1e6, lw=2.2, color=C[2 + i], label=f"σ={s:g}")
    a.axhline(n_floor / 1e6, lw=2.0, ls=(0, (4, 2)), color=C[1])
    a.annotate(f"neck floor {n_floor/1e6:.0f}M — binds whenever it is on top",
               (-39, n_floor / 1e6 * 1.4), fontsize=7.5, color=C[1])
    _ax(a, "T [°C]", "mesh nodes  [millions]",
        "C.  Mesh from Eq. 42 vs the neck-resolution floor")
    a.legend(fontsize=7.5, frameon=False, loc="lower left")

    fig.tight_layout(rect=(0, 0, 1, 0.92))
    args.out.mkdir(parents=True, exist_ok=True)
    png = args.out / "alpha_models.png"
    fig.savefig(png, dpi=150)
    print(f"plot -> {png}")

    # Console table. Panel A plots LIBB2 UNCLAMPED (the model's own behaviour);
    # panels B/C and the regime columns use the CLAMPED value, which is what a
    # run would actually see. Both are shown so the clamp's effect is visible
    # rather than silently folded in.
    print(f"\n  {'T':>5} {'σ':>9} {'LIBB2 raw':>10} {'clamped':>9} {'ARRH':>10} "
          f"{'β′/β₀':>8} {'r/L*':>8}")
    for t in (-3.0, -20.0):
        for s, _ in sig_cases:
            raw = ce.alpha_libbrecht2(t, s, clamp=None)
            al = ce.alpha_libbrecht2(t, s)
            b = ce.beta_HK(t, al)
            print(f"  {t:>5.0f} {s:>9.1e} {raw:>10.3e} {al:>9.3e} "
                  f"{ce.alpha_arrhenius(t):>10.3e} "
                  f"{1 - a1c*a2c*eps_floor/(Dm*b):>8.3f} "
                  f"{args.rneck/(b*Dv_of(t)):>8.2f}")

    A_fit, f_fit, _ = ce.libbrecht2_params()
    print(f"\n  LIBB2 fit: A = {A_fit:.4e}, f = {f_fit:.4e}")
    if A_fit >= ce.ALPHA_LIT_HI:
        print(f"  WARNING: A >= the clamp ceiling {ce.ALPHA_LIT_HI:g}, so alpha")
        print(f"  saturates at A for sigma >~ 1e-3 and its sigma-response is")
        print(f"  confined to sigma < 1e-3. This is a consequence of anchoring")
        print(f"  the fit at the BAND ENDPOINTS while also clamping to the band —")
        print(f"  the band was meant to be a consistency check, not a fit target.")
    print("\n  r/L* < 1 => attachment-limited, the regime Kuczynski's t^(1/3) is "
          "derived in.")


if __name__ == "__main__":
    main()
