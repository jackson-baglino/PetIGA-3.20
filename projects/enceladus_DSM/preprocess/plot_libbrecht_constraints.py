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

The case is made ENTIRELY INSIDE LIBBRECHT'S OWN CALIBRATION RANGE. Only the
two supersaturations the chamber experiments were run at, sigma = 1e-1 and
1e-2, are plotted. Nothing here depends on the conditions of any one of our
simulations, so nothing here can be answered with "your boundary condition is
wrong".

1. THE LAW CANNOT SIT INSIDE THE LITERATURE BAND AT ANY ONE sigma. The band
   the literature supports is 1e-3 < alpha_c < 1e-1 (Libbrecht 2017; Braun,
   Fourteau & Lowe 2024). Evaluate (L1) at Libbrecht's OWN upper chamber
   condition, sigma = 1e-1, and it returns alpha_c = 0.33 to 0.96 across
   -40..-1 C -- three to ten times ABOVE the band's ceiling. Evaluate it at the
   lower one, sigma = 1e-2, and it drops BELOW the band's floor colder than
   -29.7 C. That is structural, not a tuning problem: sigma0 spans a factor 27
   over -2..-40 C and it sits in the exponent, so no single reference sigma
   holds alpha_c inside one decade across the range.

2. K&P Eq. (45) TURNS THAT INTO MESH. beta_sub ~ 1/alpha_c, and the kinetic
   interface-width bound eps <= d0/(beta_sub*v_n) is inversely proportional to
   beta_sub. So a falling alpha_c is not merely "slow physics we can wait out".
   At sigma = 1e-2 and -40 C -- again, Libbrecht's own condition -- (L1) gives
   alpha_c = 1.7e-5, beta_sub = 4.0e9 s/m, eps = 0.40 angstrom (an EIGHTH of a
   water molecule) and Nx = 1.6e7, eleven times more nodes across the domain
   than it has molecules across it. Below about -25 C the mesh has already
   passed Nx = 1e4.

WHAT IS DELIBERATELY NOT PLOTTED
--------------------------------
Sintering is capillarity-driven, so its supersaturation is d0/rho_fillet --
lower again than either chamber value, for any micron-scale fillet. That only
makes the result worse, and earlier versions of this figure plotted it for a
particular geometry. Those curves are gone. They were specific to a single
digital experiment, and at those sigmas (L1) returns an alpha_c decades below
the literature band, so the beta_sub they imply is one we would never run --
drawing it invites the reading that we do.

A REFERENCE CASE IS UNAVOIDABLE, AND IS ONLY A SCALE
----------------------------------------------------
Eq. (45) needs a front velocity and Nx needs a domain, so some concrete case
has to be named. The defaults are a 450 x 225 um grain pair and the measured
integrated neck rate 3.416e-9 m/s from studies/molaro_2019/. They set the
SCALE of the vertical axes; they do not carry the argument, which is about
where alpha_c(T) lands relative to a band that has nothing to do with them.

Feeding v_n back self-consistently from the kinetics (--vn_feature in
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
  fig0_master.png            all six, 3 x 2, one shared legend
  fig7_legend.png            that legend alone, for laying panels out by hand
  libbrecht_constraints.csv  the numbers behind them

Panels are --panel_width wide (3 in by default, drawn compact); the master and
the legend strip are --width (10 in). 10 in is a hard ceiling -- these go on
PowerPoint slides -- and no type is below 10 pt.

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
import figstyle as fs
from figstyle import (C, INK, MUTED, GRID, BOUND_COLOR, BOUND_LABEL,
                      GRP_BIND, GRP_REF, FS_TITLE, FS_LABEL, FS_TICK, FS_LEG,
                      FS_NOTE, FS_C_TITLE, FS_C_LABEL, FS_C_TICK, SLIDE_W_IN,
                      ASPECT, COMPACT_BELOW_IN)
from figstyle import (style as _style, legend as _legend, compact as _compact,
                      lbl as _lbl, mark_offscale as _mark_offscale,
                      shade_binding as _shade_binding, sig_tex as _sig_tex)

# LIBBRECHT'S OWN CHAMBER CONDITIONS, AND ONLY THOSE.
#
# Curves for the supersaturations of one particular simulation used to sit here
# -- an imposed wall undersaturation and a neck-fillet d0/rho. They are gone for
# two reasons. They are specific to a single digital experiment, and the case
# being made is general. And at those sigmas the Libbrecht form returns an
# alpha_c decades below the literature band, so the beta_sub it implies is one
# we would never run -- plotting it invites the reading that we do.
#
# What remains is the honest general statement: evaluate Libbrecht's own law at
# Libbrecht's own conditions and follow sigma0(T) down in temperature.
#
# sigma = 1e-3 is one plain decade below the chamber range -- no experiment
# attached to it, which is the point. It shows the trend continuing without
# resting the case on conditions anyone can dispute; the legend keeps it in its
# own group so nobody reads it as something Libbrecht measured.
#
# Value, colour, label, group key (which also picks the linestyle).
SIGMA_CASES = [
    (1.0e-1, C[2], "$\\sigma = 10^{-1}$", "chamber_hi"),
    (1.0e-2, C[5], "$\\sigma = 10^{-2}$", "chamber_lo"),
    (1.0e-3, C[3], "$\\sigma = 10^{-3}$", "extrap"),
]
_LS = {"chamber_hi": "-",
       "chamber_lo": (0, (5, 2)),
       "extrap": (0, (5, 1.5, 1, 1.5))}

# Legend group headings.
GRP_LAB = "Libbrecht's chamber"
GRP_EXTRAP = "Extrapolated"
SIG_LAB = [l for _v, _c, l, k in SIGMA_CASES if k.startswith("chamber")]
SIG_EXTRAP = [l for _v, _c, l, k in SIGMA_CASES if k == "extrap"]

# a = (m/rho_ice)^(1/3): the molecular length already implied by the constants
# in comp_eps.py. Below it a continuum phase field means nothing.
A_MOL = (ce._M_H2O / ce._RHO_ICE) ** (1.0 / 3.0)

ALPHA_FLOOR = 1.0e-30       # the floor hardwired in comp_eps.alpha_libbrecht
ALPHA_CEIL = 1.0            # every impinging molecule sticks

# -------------------------------------------------------------------------
# The display window.
#
# Plotted full-range, these quantities span sixty decades and the usable
# band is a hairline. So every axis is clipped to the LITERATURE band opened
# up by two decades on each side -- roughly "everything we could plausibly
# run, plus a margin". Outside that, a curve is marked where it leaves the
# frame and the extreme it reaches is stated, so nothing is hidden; the full
# numbers stay in libbrecht_constraints.csv and the study README.
# -------------------------------------------------------------------------
VIEW_ALPHA = (ce.ALPHA_LIT_LO * 1e-2, 2.5)        # band 1e-3..1e-1, +2 decades
MF_BETA_LO, MF_BETA_HI = 2.0e4, 2.0e6             # M&F (2024) Table S1
VIEW_BETA = (1.0e3, MF_BETA_HI * 1e2)   # trimmed: nothing runs below 1e3

# Mesh window, in Nx. NX_USABLE is what we can actually run; the frame is
# opened up ~2.3 decades above it so the molecular limit Lx/a stays visible.
NX_USABLE = 1.0e4
VIEW_NX = (3.0e2, 1.0e6)


def master_groups(ctx):
    """The legend the master grid carries, and the one fig7_legend repeats.

    Panels 3-6 each draw a dotted black line for the value they actually run;
    compact mode labels all four "what we run", which is what lets one entry
    stand for the set.
    """
    return [
        (GRP_LAB, SIG_LAB),
        (GRP_EXTRAP, SIG_EXTRAP),
        (GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"]]),
        (GRP_REF, ["$\\alpha_c = 1$  physical ceiling",
                   "$N_x = 10^4$  practical ceiling",
                   "what we run", "M&F (2024) Table S1"]),
    ]


def _ls(kind):
    """Line style, so the sigmas are never told apart by colour alone."""
    return _LS[kind]


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
           label="$\\sigma_0(T)$ as coded")
    a.plot(ce._SIG0_T[:9], ce._SIG0_S[:9], "o", ms=7, color=C[0],
           markerfacecolor="white", markeredgewidth=1.8, zorder=5,
           label="table, 9 points in range")
    if not _compact():
        _panel1_notes(a)
    a.set_xlim(-41, 0.5)
    a.set_ylim(2.0e-3, 6.0e-1)
    _style(a, "T [°C]", "$\\sigma_0$  (critical supersaturation) [-]",
           "1.  Libbrecht's data: $\\sigma_0(T)$",
           short="1.  Libbrecht's $\\sigma_0(T)$", short_y="$\\sigma_0$  [-]")
    _legend(a, [("Libbrecht (2017)", ["$\\sigma_0(T)$ as coded",
                                      "table, 9 points in range"])])


def _panel1_notes(a):
    a.annotate("non-monotonic kink at $-6/-7$ °C\n(a digitisation artifact,\n"
               "not physics)", (-6.6, 6.6e-3), xytext=(-14.5, 3.0e-3),
               fontsize=FS_NOTE, color=C[1], ha="center", va="center",
               arrowprops=dict(arrowstyle="->", color=C[1], lw=1.1))
    # Upper left: the only corner with no curve, no legend and no arrow in it.
    a.text(-39.4, 5.0e-1,
           "$\\sigma_0$ spans a factor 27 across this range.\n"
           "$\\alpha_c = \\exp(-\\sigma_0/\\sigma)$, so that factor\n"
           "lands in the EXPONENT.",
           fontsize=FS_NOTE, color=INK, ha="left", va="top")


def panel_alpha_vs_sigma0(a, T, ctx):
    """2. alpha_c vs sigma0 at fixed sigma. Pure exp(-sigma0/sigma)."""
    s0 = np.logspace(-3.6, -0.7, 500)
    a.axhspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10, zorder=0)
    a.set_xscale("log")
    a.set_ylim(*VIEW_ALPHA)
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        y = np.where(s0 / sig > 69.0775, ALPHA_FLOOR, np.exp(-s0 / sig))
        y = np.maximum(y, ALPHA_FLOOR)
        a.plot(s0, y, lw=2.4, color=col, ls=_ls(kind), label=lbl)
        _mark_offscale(a, s0, y, col, row=i % 2)
    a.axhline(ALPHA_CEIL, lw=1.6, color=INK, alpha=0.7,
              label="$\\alpha_c = 1$  physical ceiling")
    # Where the table's own sigma0 values sit, so x carries physical marks.
    for Tm in ([] if _compact() else (-2.0, -20.0, -40.0)):
        s0m = ce.sigma0(Tm)
        a.axvline(s0m, lw=1.0, ls=(0, (1, 3)), color=MUTED, zorder=0)
        a.annotate(f"$T$ = {Tm:.0f} °C", (s0m, ALPHA_CEIL * 0.55), rotation=90,
                   fontsize=FS_NOTE, color=MUTED, ha="right", va="top", zorder=6,
                   bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                             alpha=0.85))
    if not _compact():
        a.annotate("literature band for $\\alpha_c$", (2.9e-4, 3.0e-2),
                   fontsize=FS_NOTE, color=C[0], va="center")
    _style(a, "$\\sigma_0(T)$  [-]", "$\\alpha_c$  [-]",
           "2.  $\\alpha_c = \\exp(-\\sigma_0/\\sigma)$ against $\\sigma_0$",
           short="2.  $\\alpha_c$ against $\\sigma_0$")
    _legend(a, [(GRP_LAB, SIG_LAB), (GRP_EXTRAP, SIG_EXTRAP),
                (GRP_REF, ["$\\alpha_c = 1$  physical ceiling"])])


def panel_alpha_vs_T(a, T, ctx):
    """3. alpha_c(T) at each sigma -- the supersaturation gap, in T."""
    a.axhspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10, zorder=0)
    a.set_ylim(*VIEW_ALPHA)
    if not _compact():
        a.annotate("literature band for $\\alpha_c$", (-39.5, 4.0e-3),
                   fontsize=FS_NOTE, color=C[0], va="center")
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        y = np.maximum(ctx["alpha"][sig], ALPHA_FLOOR)
        a.plot(T, y, lw=2.4, color=col, ls=_ls(kind), label=lbl)
        _mark_offscale(a, T, y, col, row=i % 2)
    a.axhline(ALPHA_CEIL, lw=1.6, color=INK, alpha=0.7,
              label="$\\alpha_c = 1$  physical ceiling")
    run = _lbl(f"$\\alpha_c$ = {ctx['alpha_run']:g}  what we run", "what we run")
    a.axhline(ctx["alpha_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "T [°C]", "$\\alpha_c$  [-]",
           "3.  $\\alpha_c(T)$: two decades down in $\\sigma$, "
           "thirty down in $\\alpha_c$", short="3.  $\\alpha_c(T)$")
    _legend(a, [(GRP_LAB, SIG_LAB), (GRP_EXTRAP, SIG_EXTRAP),
                (GRP_REF, ["$\\alpha_c = 1$  physical ceiling", run])])


def panel_beta_vs_T(a, T, ctx):
    """4. beta_sub(T). The defining equation rides in the title, not the frame,
    so it cannot land on the legend or on a curve."""
    a.axhspan(MF_BETA_LO, MF_BETA_HI, color=C[2], alpha=0.14, zorder=0,
              label="M&F (2024) Table S1")
    a.set_ylim(*VIEW_BETA)
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        a.plot(T, ctx["beta"][sig], lw=2.4, color=col, ls=_ls(kind), label=lbl)
        _mark_offscale(a, T, ctx["beta"][sig], col, row=i % 2)
    run = _lbl(f"$\\beta_{{sub}}$ = {ctx['beta_run']:.1e} s/m  what we run",
               "what we run")
    a.axhline(ctx["beta_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "T [°C]", "$\\beta_{sub}$  [s/m]",
           "4.  Attachment resistance $\\beta_{sub} \\propto 1/\\alpha_c$\n\n"
           r"$\beta_{sub}(T)=\dfrac{\beta_{HK}}{\rho_{vs}(T)/\rho_i}"
           r"=\dfrac{\rho_i}{\rho_{vs}(T)}\;\dfrac{1}{\alpha_c}\;"
           r"\sqrt{\dfrac{2\pi m}{k_B T}}\,,\qquad"
           r"\alpha_c=\exp\left[-\sigma_0(T)/\sigma\right]$",
           title_size=13, pad=14,
           short="4.  $\\beta_{sub}(T) \\propto 1/\\alpha_c$")
    _legend(a, [(GRP_LAB, SIG_LAB), (GRP_EXTRAP, SIG_EXTRAP),
                (GRP_REF, ["M&F (2024) Table S1", run])])


def panel_eps_vs_T(a, T, ctx):
    """5. eps(T) straight out of comp_eps.compute_eps, shaded by binding bound."""
    prim = ctx["primary_sigma"]
    a.set_ylim(*ctx["view_eps"])
    _shade_binding(a, T, ctx["binding"][prim], y_frac=0.975)
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        a.plot(T, ctx["eps"][sig], lw=3.4 if sig == prim else 2.4, color=col,
               ls=_ls(kind),
               label=lbl)
        _mark_offscale(a, T, ctx["eps"][sig], col, row=i % 2)
    ceil = _lbl(f"$N_x = 10^4$ ceiling, $\\epsilon$ = "
                f"{ctx['eps_usable']*1e9:.0f} nm", "$N_x = 10^4$  practical ceiling")
    a.axhline(ctx["eps_usable"], lw=1.6, color=MUTED, label=ceil)
    run = _lbl(f"production run, $\\epsilon$ = {ctx['eps_run']*1e6:.3g} µm",
               "what we run")
    a.axhline(ctx["eps_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "T [°C]", "$\\epsilon$  [m]",
           "5.  Interface width from the comp_eps.py bounds "
           f"(shading follows $\\sigma = {_sig_tex(prim)}$)",
           short="5.  Interface width $\\epsilon(T)$")
    _legend(a, [(GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"]]),
                (GRP_LAB, SIG_LAB), (GRP_EXTRAP, SIG_EXTRAP),
                (GRP_REF, [ceil, run])])


def panel_Nx_vs_T(a, T, ctx):
    """6. Nx(T) = ceil(sqrt(2)*Lx/eps), against what a machine can hold."""
    prim = ctx["primary_sigma"]
    a.set_ylim(*VIEW_NX)
    _shade_binding(a, T, ctx["binding"][prim], y_frac=0.03)
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        a.plot(T, ctx["Nx"][sig], lw=3.4 if sig == prim else 2.4, color=col,
               ls=_ls(kind), label=lbl)
        _mark_offscale(a, T, ctx["Nx"][sig], col, fmt="{:.0e}", row=i % 2)
    a.axhline(NX_USABLE, lw=1.6, color=MUTED,
              label="$N_x = 10^4$  practical ceiling")
    run = _lbl(f"production mesh, $N_x$ = {ctx['Nx_run']}", "what we run")
    a.axhline(ctx["Nx_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "T [°C]", "$N_x$  [nodes across $L_x$]",
           "6.  The mesh Libbrecht's $\\alpha_c$ demands",
           short="6.  Mesh $N_x(T)$", short_y="$N_x$")
    _legend(a, [(GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"]]),
                (GRP_LAB, SIG_LAB), (GRP_EXTRAP, SIG_EXTRAP),
                (GRP_REF, ["$N_x = 10^4$  practical ceiling", run])])


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
    # eps window is the exact reciprocal of VIEW_NX, so panels 5 and 6 are
    # the same plot on two axes and can be read against each other.
    ctx["view_eps"] = (math.sqrt(2.0) * args.Lx / VIEW_NX[1],
                       math.sqrt(2.0) * args.Lx / VIEW_NX[0])
    ctx["eps_usable"] = math.sqrt(2.0) * args.Lx / NX_USABLE
    ctx["nodes_run"] = ctx["Nx_run"] * math.ceil(math.sqrt(2.0) * args.Ly
                                                 / ctx["eps_run"])
    return ctx


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Geometry: Molaro dom2, the production axisymmetric grain pair.
    # Reference-case geometry: the Molaro dom2 production domain,
    # inputs/geometry/molaro/molaro_2D_L450x225um_..._dom2.opts. Lx enters only
    # as Nx = ceil(sqrt(2)*Lx/eps), i.e. LINEARLY -- it sets where the vertical
    # axis of panel 6 sits, not the shape of any curve on it.
    p.add_argument("--Lx", type=float, default=4.50e-4, help="domain x [m]")
    p.add_argument("--Ly", type=float, default=2.25e-4, help="domain y [m]")
    p.add_argument("--Rave", type=float, default=8.675e-5,
                   help="representative grain radius [m], for B-CURV")
    p.add_argument("--vn", type=float, default=3.415598e-9,
                   help="MEASURED front velocity [m/s]; default is the "
                        "integrated Molaro Fig. 11 neck rate from "
                        "studies/molaro_2019/alpha_c_estimate.csv")
    p.add_argument("--safety", type=float, default=0.5)
    p.add_argument("--eps_over_R", type=float, default=0.05)
    p.add_argument("--primary_sigma", type=float, default=1.0e-2,
                   help="the sigma whose binding bound is shaded in figs 5-6; "
                        "must be one of the plotted SIGMA_CASES")
    p.add_argument("--alpha_run", type=float, default=0.1,
                   help="the constant alpha_c the production runs use")
    p.add_argument("--T_run", type=float, default=-20.0)
    p.add_argument("--eps_run", type=float, default=1.18e-7,
                   help="eps of the production geometry [m]; 0 = size it here")
    p.add_argument("--Tmin", type=float, default=-40.0)
    p.add_argument("--Tmax", type=float, default=-1.0)
    p.add_argument("--nT", type=int, default=400)
    p.add_argument("--width", type=float, default=SLIDE_W_IN,
                   help="width [in] of the master and the legend strip. Hard "
                        "ceiling: these go on slides, so nothing is written "
                        "wider than this")
    p.add_argument("--panel_width", type=float, default=3.0,
                   help="width [in] of each individual panel. Below "
                        "COMPACT_BELOW_IN they are drawn in compact mode -- "
                        "short titles, no prose annotations, no per-panel "
                        "legend (fig7_legend.png carries it for the set). "
                        "Pass 10 to get the fully annotated panels back")
    p.add_argument("--dpi", type=int, default=200,
                   help="200 keeps text crisp on a projector. PowerPoint that "
                        "ignores the PNG's DPI tag will insert it at width*dpi/96 "
                        "in -- set the width box back to --width, or pass "
                        "--dpi 96 for drop-in sizing at some loss of sharpness")
    p.add_argument("--panel_legend", choices=("separate", "inline", "none"),
                   default="separate",
                   help="a 3 in panel has no room for a legend inside the "
                        "frame. 'separate' writes <panel>_legend.png next to "
                        "each one; 'inline' makes the panel taller and puts "
                        "the legend under the axes in the same image; 'none' "
                        "leaves only the combined fig7_legend.png")
    p.add_argument("--label", default="reference case: 450 × 225 µm grain pair")
    p.add_argument("--out", type=Path, default=Path("studies/libbrecht_kinetics"))
    args = p.parse_args()

    T = np.linspace(args.Tmin, args.Tmax, args.nT)
    ctx = build_context(T, args)
    args.out.mkdir(parents=True, exist_ok=True)

    fs.emit(args.out, PANELS, T, ctx, width=args.width,
            panel_width=args.panel_width, dpi=args.dpi,
            panel_legend=args.panel_legend, master_groups=master_groups,
            suptitle="Libbrecht (2017) $\\sigma_0(T)$ as the source of "
                     "$\\alpha_c$: the chain down to $\\epsilon$ and the mesh\n"
                     f"{args.label}, measured $v_n$ = {args.vn:.3g} m/s")

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
