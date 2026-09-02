#!/usr/bin/env python3
"""Why Libbrecht (2017) sigma0(T) cannot be used to set alpha_c in this solver.

WHAT THIS FIGURE SET ARGUES
---------------------------
Libbrecht, Annu. Rev. Mater. Res. 47, 271 (2017), reports a critical
supersaturation sigma0(T) for ice crystal growth, together with the nucleation
form for the attachment (condensation) coefficient

    alpha_c(T, sigma) = exp(-sigma0(T)/sigma).                            (L1)

Taking (L1) at face value is attractive: it makes alpha_c -- the model's one
genuinely free parameter -- a function of state rather than a fitted constant.
This script shows, using the SAME code path the mesh sizer uses
(comp_eps.compute_eps), that it does not survive contact with the problem.

The chain is short and every link is plotted:

    sigma0(T)  ->  alpha_c = exp(-sigma0/sigma)  ->  beta_sub ~ 1/alpha_c
               ->  eps <= d0/(beta_sub*v_n)      ->  Nx = ceil(sqrt(2)*Lx/eps)

Two facts drive the whole result.

1. THE SUPERSATURATION GAP. Libbrecht's chamber experiments are run at
   sigma ~ 1e-2 to 1e-1. Our sintering problem runs one to two decades lower:
   the Molaro wall undersaturation is 1-h = 2.7e-3, and the Gibbs-Thomson
   supersaturation at a neck fillet is d0/rho_fillet ~ 4.5e-4. Because (L1) is
   an exponential in sigma0/sigma, moving down two decades in sigma does not
   reduce alpha_c by two decades -- it annihilates it. At -20 C and
   sigma = 4.5e-4, (L1) returns alpha_c = 1e-30 (the code's own underflow
   floor), i.e. no sintering at all.

2. K&P Eq. (45) TURNS THAT INTO MESH. beta_sub ~ 1/alpha_c, and the kinetic
   interface-width bound eps <= d0/(beta_sub*v_n) is inversely proportional to
   beta_sub. A vanishing alpha_c is therefore not merely "slow physics we can
   wait out": it is a demand for an interface width below the size of a water
   molecule, and a mesh with more nodes than the domain has molecules.

THE MEASURED v_n IS THE LOAD-BEARING INPUT
------------------------------------------
Eq. (45) is only a constraint if the front actually moves. We take v_n from
Molaro et al. (2019) Fig. 11 themselves -- the integrated neck-growth velocity
3.416e-9 m/s in studies/molaro_2019/alpha_c_estimate.csv. So the argument is
not "assume a velocity and watch the mesh explode"; it is "the experiment we
are replicating moves its interface at THIS speed, and Libbrecht's alpha_c
says that is impossible to resolve."

Feeding v_n back self-consistently from the kinetics instead (--vn_feature in
comp_eps.py) makes Eq. (45) collapse to eps <= R_feat and hides the problem;
that is circular, not a rebuttal. See studies/libbrecht_kinetics/README.md.

WHAT THE PLOTS ARE
------------------
  fig1_sigma0_vs_T.png       sigma0(T): Libbrecht's data as this repo codes it
  fig2_alpha_vs_sigma0.png   alpha_c vs sigma0, at fixed sigma
  fig3_alpha_vs_T.png        alpha_c(T) at each sigma
  fig4_beta_sub_vs_T.png     beta_sub(T), with its defining equation
  fig5_eps_vs_T.png          eps(T) from comp_eps, shaded by binding bound
  fig6_Nx_vs_T.png           Nx(T), against what a machine can actually hold
  fig0_overview.png          all six on one slide
  libbrecht_constraints.csv  the numbers behind them

Usage:
    python preprocess/plot_libbrecht_constraints.py \
        --out studies/libbrecht_kinetics
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

# Okabe-Ito, the repo standard (plot_alpha_kinetics.py, fit_neck_growth.py).
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

# Colour per binding bound, shared by the fig5/fig6 shading.
BOUND_COLOR = {"B-HEAT": C[0], "B-VAPOR": C[5], "B-KINETIC": C[1], "B-CURV": C[2]}
BOUND_LABEL = {
    "B-HEAT":    r"B-HEAT binds: $\epsilon\leq s\,D_{heat}\beta_{HK}$ (K&P 43a/b)",
    "B-VAPOR":   r"B-VAPOR binds: $\epsilon\leq s\,D_v\beta_{HK}$ (K&P 43c)",
    "B-KINETIC": r"B-KINETIC binds: $\epsilon\leq s\,d_0/(\beta_{sub}v_n)$ (K&P 45)",
    "B-CURV":    r"B-CURV binds: $\epsilon\leq 0.05\,R_{ave}$ (geometric)",
}

# The supersaturations in play. Value, colour, label, and whether it is a
# condition WE impose or one LIBBRECHT measured at.
SIGMA_CASES = [
    (1.0e-1, C[2], "Libbrecht chamber, $\\sigma = 10^{-1}$", "lab"),
    (1.0e-2, C[5], "Libbrecht chamber, $\\sigma = 10^{-2}$", "lab"),
    (2.7e-3, C[3], "Molaro wall BC, $\\sigma = 2.7\\times10^{-3}$", "ours"),
    (4.5e-4, C[1], "neck fillet $d_0/\\rho$, $\\sigma = 4.5\\times10^{-4}$", "ours"),
]

# a = (m/rho_ice)^(1/3): the molecular length already implied by the constants
# in comp_eps.py. Below it a continuum phase field means nothing.
A_MOL = (ce._M_H2O / ce._RHO_ICE) ** (1.0 / 3.0)

ALPHA_FLOOR = 1.0e-30       # the floor hardwired in comp_eps.alpha_libbrecht


def _style(a, xlabel, ylabel, title, logy=True):
    a.set_xlabel(xlabel, fontsize=10)
    a.set_ylabel(ylabel, fontsize=10)
    if title:
        a.set_title(title, fontsize=11, color=INK, loc="left")
    if logy:
        a.set_yscale("log")
    a.grid(alpha=0.25, lw=0.5, color=GRID, which="both")
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=9)


def _ls(kind):
    """Ours = solid, Libbrecht's own conditions = dashed. Never colour alone."""
    return "-" if kind == "ours" else (0, (5, 2))


# =========================================================================
# Panels.  Each takes a bare Axes so it can be drawn standalone OR into the
# overview grid -- one code path, so the slide and the summary cannot drift.
# =========================================================================

def panel_sigma0(a, T, ctx):
    """1. sigma0(T) -- Libbrecht's own plot, as comp_eps.sigma0 codes it."""
    # Own T grid: the table's warmest point sits at T -> 0, so stopping at the
    # shared Tmax would leave that marker floating off the end of the curve.
    Tf = np.linspace(-40.0, -1.0e-4, 600)
    a.plot(Tf, [ce.sigma0(t) for t in Tf], lw=2.4, color=C[0],
           label="$\\sigma_0(T)$ as coded (log-log interpolation in $|T|$)")
    a.plot(ce._SIG0_T[:9], ce._SIG0_S[:9], "o", ms=7, color=C[0],
           markerfacecolor="white", markeredgewidth=1.8, zorder=5,
           label="Libbrecht (2017) table, 9 points in range")
    a.annotate("non-monotonic kink at $-6/-7$ °C\n(a digitisation artifact,\n"
               "not physics)", (-6.6, 6.6e-3), xytext=(-13.5, 3.0e-3),
               fontsize=8.5, color=C[1], ha="center", va="center",
               arrowprops=dict(arrowstyle="->", color=C[1], lw=1.1))
    a.text(-39.4, 2.35e-3,
           "$\\sigma_0$ spans a factor 27 across this range.\n"
           "$\\alpha_c = \\exp(-\\sigma_0/\\sigma)$, so that factor\n"
           "lands in the EXPONENT.",
           fontsize=8.5, color=INK, ha="left", va="bottom")
    a.set_xlim(-41, 0.5)
    a.set_ylim(2.0e-3, 1.0)
    _style(a, "T [°C]", "$\\sigma_0$  (critical supersaturation) [-]",
           "1.  Libbrecht's data: $\\sigma_0(T)$")
    a.legend(fontsize=8.5, frameon=False, loc="upper right")


def panel_alpha_vs_sigma0(a, T, ctx):
    """2. alpha_c vs sigma0 at fixed sigma. Pure exp(-sigma0/sigma)."""
    s0 = np.logspace(-3.6, -0.7, 500)
    a.axhspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10, zorder=0)
    for sig, col, lbl, kind in SIGMA_CASES:
        y = np.where(s0 / sig > 69.0775, ALPHA_FLOOR, np.exp(-s0 / sig))
        a.plot(s0, np.maximum(y, ALPHA_FLOOR), lw=2.4, color=col,
               ls=_ls(kind), label=lbl)
    a.axhline(ALPHA_FLOOR, lw=1.4, color=INK, alpha=0.6,
              label="$10^{-30}$ underflow floor (material_properties.c)")
    # Where the table's own sigma0 values sit, so x carries physical marks.
    for Tm in (-2.0, -20.0, -40.0):
        s0m = ce.sigma0(Tm)
        a.axvline(s0m, lw=1.0, ls=(0, (1, 3)), color=MUTED, zorder=0)
        a.annotate(f"$T$ = {Tm:.0f} °C", (s0m, 6.0e2), rotation=90, fontsize=8,
                   color=MUTED, ha="right", va="top")
    a.annotate("literature band for $\\alpha_c$",
               (2.9e-4, 3.0e-2), fontsize=8, color=C[0], va="center")
    a.set_xscale("log")
    a.set_ylim(1e-34, 1e3)
    _style(a, "$\\sigma_0(T)$  [-]", "$\\alpha_c$  [-]",
           "2.  $\\alpha_c = \\exp(-\\sigma_0/\\sigma)$ against $\\sigma_0$")
    a.legend(fontsize=8, frameon=False, loc="lower left")


def panel_alpha_vs_T(a, T, ctx):
    """3. alpha_c(T) at each sigma -- the supersaturation gap, in T."""
    a.axhspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10, zorder=0)
    a.annotate("literature band for $\\alpha_c$", (-39.5, 5.0e-3), fontsize=8,
               color=C[0], va="center")
    for sig, col, lbl, kind in SIGMA_CASES:
        a.plot(T, np.maximum(ctx["alpha"][sig], ALPHA_FLOOR), lw=2.4, color=col,
               ls=_ls(kind), label=lbl)
    a.axhline(ALPHA_FLOOR, lw=1.4, color=INK, alpha=0.6,
              label="$10^{-30}$ underflow floor: no sintering at all")
    a.axhline(ctx["alpha_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"$\\alpha_c$ = {ctx['alpha_run']:g}, the constant we run")
    a.set_ylim(1e-34, 1e3)
    _style(a, "T [°C]", "$\\alpha_c$  [-]",
           "3.  $\\alpha_c(T)$: two decades down in $\\sigma$, "
           "thirty down in $\\alpha_c$")
    a.legend(fontsize=8, frameon=False, loc="lower right")


def panel_beta_vs_T(a, T, ctx):
    """4. beta_sub(T), with the equation that produces it."""
    a.axhspan(2.0e4, 2.0e6, color=C[2], alpha=0.14, zorder=0,
              label="Moure & Fu (2024) Table S1: $2\\times10^4$–"
                    "$2\\times10^6$ s/m")
    for sig, col, lbl, kind in SIGMA_CASES:
        a.plot(T, ctx["beta"][sig], lw=2.4, color=col, ls=_ls(kind), label=lbl)
    a.axhline(ctx["beta_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"$\\beta_{{sub}}$ = {ctx['beta_run']:.2e} s/m at "
                    f"$\\alpha_c$ = {ctx['alpha_run']:g} (what we run)")
    a.text(0.030, 0.965,
           r"$\beta_{sub}(T)\;=\;\dfrac{\beta_{HK}}{\rho_{vs}(T)/\rho_i}"
           r"\;=\;\dfrac{\rho_i}{\rho_{vs}(T)}\;\dfrac{1}{\alpha_c}\;"
           r"\sqrt{\dfrac{2\pi m}{k_B T}}\,,$"
           "\n" r"$\qquad\alpha_c\;=\;\exp\left[-\sigma_0(T)/\sigma\right]$"
           "\n" r"(K&P $\beta_0\equiv$ M&F $\beta_{sub}$; $\beta_{HK}$ is the "
           r"scaled Hertz–Knudsen form)",
           transform=a.transAxes, fontsize=9.5, color=INK, va="top", ha="left",
           bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=GRID, lw=0.9),
           zorder=6)
    a.set_ylim(1e2, 1e46)
    _style(a, "T [°C]", "$\\beta_{sub}$  [s/m]",
           "4.  Attachment resistance $\\beta_{sub} \\propto 1/\\alpha_c$")
    a.legend(fontsize=8, frameon=False, loc="upper right")


def panel_eps_vs_T(a, T, ctx):
    """5. eps(T) straight out of comp_eps.compute_eps, shaded by binding bound."""
    prim = ctx["primary_sigma"]
    _shade_binding(a, T, ctx["binding"][prim])
    for sig, col, lbl, kind in SIGMA_CASES:
        a.plot(T, ctx["eps"][sig], lw=3.0 if sig == prim else 2.0, color=col,
               ls=_ls(kind),
               label=lbl + ("  ← shaded" if sig == prim else ""))
    a.axhline(ctx["eps_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"$\\epsilon$ = {ctx['eps_run']*1e6:.3g} µm, the production run")
    a.axhline(A_MOL, lw=1.8, color=INK, alpha=0.8,
              label=f"one water molecule, $a=(m/\\rho_i)^{{1/3}}$ = "
                    f"{A_MOL*1e10:.2f} Å")
    a.set_ylim(1e-38, 1e-2)
    _style(a, "T [°C]", "$\\epsilon$  [m]",
           "5.  Interface width from the comp_eps.py bounds")
    a.legend(fontsize=7.6, frameon=False, loc="lower right")


def panel_Nx_vs_T(a, T, ctx):
    """6. Nx(T) = ceil(sqrt(2)*Lx/eps), against what a machine can hold."""
    prim = ctx["primary_sigma"]
    _shade_binding(a, T, ctx["binding"][prim])
    for sig, col, lbl, kind in SIGMA_CASES:
        a.plot(T, ctx["Nx"][sig], lw=3.0 if sig == prim else 2.0, color=col,
               ls=_ls(kind), label=lbl)
    a.axhline(ctx["Nx_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"$N_x$ = {ctx['Nx_run']}, the production mesh "
                    f"({ctx['nodes_run']/1e6:.1f} M nodes in 2D)")
    n_atom = ctx["Lx"] / A_MOL
    a.axhline(n_atom, lw=1.8, color=INK, alpha=0.8,
              label=f"$L_x/a$ = {n_atom:.2g}: one node per water molecule")
    a.set_ylim(1e2, 1e36)
    _style(a, "T [°C]", "$N_x$  [nodes across $L_x$]",
           "6.  The mesh Libbrecht's $\\alpha_c$ demands")
    a.legend(fontsize=7.6, frameon=False, loc="upper right",
             bbox_to_anchor=(1.0, 0.93))


def _shade_binding(a, T, binding):
    """Shade contiguous T-runs sharing a binding bound; name each region once."""
    seen = set()
    i = 0
    while i < len(T):
        j = i
        while j + 1 < len(T) and binding[j + 1] == binding[i]:
            j += 1
        name = binding[i]
        a.axvspan(T[i], T[j], color=BOUND_COLOR.get(name, MUTED), alpha=0.12,
                  lw=0, zorder=0,
                  label=BOUND_LABEL[name] if name not in seen else None)
        seen.add(name)
        # Name the region inside it, so the shading is never colour-alone.
        a.annotate(name, (0.5 * (T[i] + T[j]), 0.975),
                   xycoords=("data", "axes fraction"), fontsize=9,
                   color=BOUND_COLOR.get(name, MUTED), ha="center", va="top",
                   fontweight="bold")
        i = j + 1


PANELS = [
    ("fig1_sigma0_vs_T",     panel_sigma0),
    ("fig2_alpha_vs_sigma0", panel_alpha_vs_sigma0),
    ("fig3_alpha_vs_T",      panel_alpha_vs_T),
    ("fig4_beta_sub_vs_T",   panel_beta_vs_T),
    ("fig5_eps_vs_T",        panel_eps_vs_T),
    ("fig6_Nx_vs_T",         panel_Nx_vs_T),
]


# =========================================================================

def build_context(T, args):
    """Evaluate the whole chain, for every sigma, through comp_eps itself."""
    ctx = dict(primary_sigma=args.primary_sigma, Lx=args.Lx, Ly=args.Ly,
               alpha_run=args.alpha_run)
    ctx["sig0"] = np.array([ce.sigma0(t) for t in T])
    for k in ("alpha", "beta", "eps", "Nx", "nodes", "binding"):
        ctx[k] = {}

    for sig, _c, _l, _k in SIGMA_CASES:
        al, be, ep, nx, nn, bd = [], [], [], [], [], []
        for t in T:
            a_c = ce.alpha_libbrecht(t, sig)
            r = ce.compute_eps(Lx=args.Lx, Ly=args.Ly, Lz=0.0, Rave=args.Rave,
                               T0_C=t, alpha_c=a_c, v_n=args.vn,
                               safety=args.safety, eps_over_R=args.eps_over_R)
            al.append(a_c); be.append(r["beta_uns"]); ep.append(r["eps"])
            nx.append(float(r["Nx"])); nn.append(float(r["Nx"]) * float(r["Ny"]))
            bd.append(r["binding"])
        ctx["alpha"][sig] = np.array(al); ctx["beta"][sig] = np.array(be)
        ctx["eps"][sig] = np.array(ep);   ctx["Nx"][sig] = np.array(nx)
        ctx["nodes"][sig] = np.array(nn); ctx["binding"][sig] = bd

    # The run we actually do: a constant alpha_c, sized at the run temperature.
    run = ce.compute_eps(Lx=args.Lx, Ly=args.Ly, Lz=0.0, Rave=args.Rave,
                         T0_C=args.T_run, alpha_c=args.alpha_run, v_n=args.vn,
                         safety=args.safety, eps_over_R=args.eps_over_R)
    ctx["beta_run"] = run["beta_uns"]
    ctx["eps_run"] = args.eps_run if args.eps_run > 0 else run["eps"]
    ctx["Nx_run"] = math.ceil(math.sqrt(2.0) * args.Lx / ctx["eps_run"])
    ctx["nodes_run"] = ctx["Nx_run"] * math.ceil(math.sqrt(2.0) * args.Ly
                                                 / ctx["eps_run"])
    return ctx


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Geometry: Molaro dom2, the production axisymmetric grain pair.
    p.add_argument("--Lx", type=float, default=4.50e-4)
    p.add_argument("--Ly", type=float, default=2.25e-4)
    p.add_argument("--Rave", type=float, default=8.675e-5)
    p.add_argument("--vn", type=float, default=3.415598e-9,
                   help="MEASURED front velocity [m/s]; default is the "
                        "integrated Molaro Fig. 11 neck rate from "
                        "studies/molaro_2019/alpha_c_estimate.csv")
    p.add_argument("--safety", type=float, default=0.5)
    p.add_argument("--eps_over_R", type=float, default=0.05)
    p.add_argument("--primary_sigma", type=float, default=2.7e-3,
                   help="the sigma whose binding bound is shaded in figs 5-6")
    p.add_argument("--alpha_run", type=float, default=0.1,
                   help="the constant alpha_c the production runs use")
    p.add_argument("--T_run", type=float, default=-20.0)
    p.add_argument("--eps_run", type=float, default=1.18e-7,
                   help="eps of the production geometry [m]; 0 = size it here")
    p.add_argument("--Tmin", type=float, default=-40.0)
    p.add_argument("--Tmax", type=float, default=-1.0)
    p.add_argument("--nT", type=int, default=400)
    p.add_argument("--label", default="Molaro grain pair, 450 × 225 µm")
    p.add_argument("--out", type=Path, default=Path("studies/libbrecht_kinetics"))
    args = p.parse_args()

    T = np.linspace(args.Tmin, args.Tmax, args.nT)
    ctx = build_context(T, args)
    args.out.mkdir(parents=True, exist_ok=True)

    # --- standalone, slide-sized figures ---------------------------------
    for name, fn in PANELS:
        fig, a = plt.subplots(figsize=(8.2, 5.8))
        fn(a, T, ctx)
        fig.tight_layout()
        fig.savefig(args.out / f"{name}.png", dpi=200)
        plt.close(fig)
        print(f"plot -> {args.out / (name + '.png')}")

    # --- one-slide overview ----------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(20.0, 11.4))
    fig.suptitle("Libbrecht (2017) $\\sigma_0(T)$ as the source of $\\alpha_c$: "
                 f"the chain down to $\\epsilon$ and the mesh — {args.label}, "
                 f"measured $v_n$ = {args.vn:.3g} m/s",
                 fontsize=14, color=INK)
    for ax, (_n, fn) in zip(axes.ravel(), PANELS):
        fn(ax, T, ctx)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(args.out / "fig0_overview.png", dpi=170)
    plt.close(fig)
    print(f"plot -> {args.out / 'fig0_overview.png'}")

    # --- the numbers ------------------------------------------------------
    csv = args.out / "libbrecht_constraints.csv"
    with open(csv, "w") as fh:
        fh.write("T_C,sigma0,sigma,alpha_c,beta_sub_s_per_m,eps_m,binding,"
                 "Nx,nodes_2D\n")
        for i in range(0, len(T), max(1, len(T) // 80)):
            for sig, _c, _l, _k in SIGMA_CASES:
                fh.write(f"{T[i]:.3f},{ctx['sig0'][i]:.6e},{sig:.3e},"
                         f"{ctx['alpha'][sig][i]:.6e},{ctx['beta'][sig][i]:.6e},"
                         f"{ctx['eps'][sig][i]:.6e},{ctx['binding'][sig][i]},"
                         f"{ctx['Nx'][sig][i]:.6e},{ctx['nodes'][sig][i]:.6e}\n")
    print(f"csv  -> {csv}")

    # --- the headline numbers, for the slide text -------------------------
    print(f"\n  Production run:  alpha_c = {args.alpha_run:g}, "
          f"beta_sub = {ctx['beta_run']:.3e} s/m, eps = {ctx['eps_run']:.3e} m, "
          f"Nx = {ctx['Nx_run']}, {ctx['nodes_run']/1e6:.1f} M nodes")
    print(f"  Molecular limit: a = {A_MOL:.3e} m, Lx/a = {args.Lx/A_MOL:.3e}\n")
    print(f"  {'T':>5} {'sigma':>9} {'alpha_c':>11} {'beta_sub':>11} "
          f"{'eps [m]':>11} {'bind':>10} {'Nx':>11} {'nodes':>11}")
    for t in (-2.0, -10.0, -20.0, -40.0):
        i = int(np.argmin(abs(T - t)))
        for sig, _c, _l, _k in SIGMA_CASES:
            print(f"  {T[i]:>5.1f} {sig:>9.1e} {ctx['alpha'][sig][i]:>11.3e} "
                  f"{ctx['beta'][sig][i]:>11.3e} {ctx['eps'][sig][i]:>11.3e} "
                  f"{ctx['binding'][sig][i]:>10s} {ctx['Nx'][sig][i]:>11.3e} "
                  f"{ctx['nodes'][sig][i]:>11.3e}")
        print()


if __name__ == "__main__":
    main()
