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
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comp_eps as ce

# Okabe-Ito, the repo standard (plot_alpha_kinetics.py, fit_neck_growth.py).
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

# Colour per binding bound, shared by the fig5/fig6 shading.
BOUND_COLOR = {"B-HEAT": C[0], "B-VAPOR": C[5], "B-KINETIC": C[1], "B-CURV": C[2]}
# Short, because these ride in a legend under the axes. The formulas they
# used to carry live in docs/tex/constraints_iguanatex.txt.
BOUND_LABEL = {
    "B-HEAT":    "B-HEAT  (K&P 43a/b)",
    "B-VAPOR":   "B-VAPOR  (K&P 43c)",
    "B-KINETIC": "B-KINETIC  (K&P 45)",
    "B-CURV":    "B-CURV  (geometric)",
}

# The supersaturations in play. Value, colour, label, and whether it is a
# condition WE impose or one LIBBRECHT measured at.
SIGMA_CASES = [
    (1.0e-1, C[2], "$\\sigma = 10^{-1}$", "lab"),
    (1.0e-2, C[5], "$\\sigma = 10^{-2}$", "lab"),
    (2.7e-3, C[3], "$\\sigma = 2.7\\times10^{-3}$   wall BC", "ours"),
    (4.5e-4, C[1], "$\\sigma = 4.5\\times10^{-4}$   neck fillet", "ours"),
]

# Legend group headings. Which sigma a curve belongs to is the first thing a
# reader needs, and it is exactly what the entry label cannot say briefly --
# so the heading says it once for the pair.
GRP_LAB = "Libbrecht's chamber"
GRP_OURS = "Our conditions"
GRP_BIND = "Which bound binds"
GRP_REF = "Reference"
SIG_LAB = [l for _v, _c, l, k in SIGMA_CASES if k == "lab"]
SIG_OURS = [l for _v, _c, l, k in SIGMA_CASES if k == "ours"]

# a = (m/rho_ice)^(1/3): the molecular length already implied by the constants
# in comp_eps.py. Below it a continuum phase field means nothing.
A_MOL = (ce._M_H2O / ce._RHO_ICE) ** (1.0 / 3.0)

ALPHA_FLOOR = 1.0e-30       # the floor hardwired in comp_eps.alpha_libbrecht
ALPHA_CEIL = 1.0            # every impinging molecule sticks

# Type scale. These panels are meant to fill most of a slide, so they are
# sized for a projector rather than for a screen an arm's length away. One
# scheme serves both the standalone figures and the overview grid, because
# FIGSIZE and OVERVIEW_SIZE keep the per-panel area roughly equal.
# Nothing here is below 10 pt.
FS_TITLE, FS_LABEL, FS_TICK, FS_LEG, FS_NOTE = 15, 13, 11.5, 10, 10
# 10 x 5 of plot, plus an inch for the legend strip underneath it.
FIGSIZE = (10.0, 5.8)
# One standalone panel per grid cell: anything smaller and the legends -- which
# sit UNDER the axes -- overflow into the neighbouring column.
OVERVIEW_SIZE = (3 * FIGSIZE[0], 2 * FIGSIZE[1])

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
VIEW_NX = (4.0e2, 1.0e6)


def _mark_offscale(a, x, y, color, fmt="{:.1e}", row=0):
    """Flag where a curve leaves the frame, and say how far off scale it goes.

    Clipping to the usable window is the point of the tight limits, but a
    curve that simply vanishes at the frame reads as missing data rather
    than as a number too large to run. One marker on the edge plus the
    extreme value fixes that without restoring sixty decades of axis.
    """
    lo, hi = a.get_ylim()
    for out, edge, mk, va, extreme, arrow in (
            (y > hi, hi, "^", "top", np.nanmax, "↑"),
            (y < lo, lo, "v", "bottom", np.nanmin, "↓")):
        if not out.any() or out.all():
            continue
        # Mark the run that actually contains the extreme, at the end of that
        # run touching the frame -- a curve can leave and re-enter (the sigma0
        # kink does exactly that), and a marker parked in the middle of some
        # other excursion points at the wrong place.
        k = int(np.nanargmax(y) if va == "top" else np.nanargmin(y))
        lo_i = hi_i = k
        while lo_i > 0 and out[lo_i - 1]:
            lo_i -= 1
        while hi_i < len(out) - 1 and out[hi_i + 1]:
            hi_i += 1
        j = hi_i if hi_i < len(out) - 1 else (lo_i if lo_i > 0 else k)
        pad = 8.0 + 15.0 * row
        a.plot([x[j]], [edge], mk, ms=9, color=color, clip_on=False, zorder=7)
        a.annotate(f"{arrow} {fmt.format(extreme(y))}", (x[j], edge),
                   xytext=(0, -pad if va == "top" else pad),
                   textcoords="offset points", fontsize=FS_NOTE, color=color,
                   ha="center", va=va, zorder=7,
                   bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                             alpha=0.85))


def _style(a, xlabel, ylabel, title, logy=True, title_size=FS_TITLE, pad=10):
    a.set_xlabel(xlabel, fontsize=FS_LABEL)
    a.set_ylabel(ylabel, fontsize=FS_LABEL)
    if title:
        a.set_title(title, fontsize=title_size, color=INK, loc="left", pad=pad)
    if logy:
        a.set_yscale("log")
    a.grid(alpha=0.25, lw=0.5, color=GRID, which="both")
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=FS_TICK)


def _blank():
    """An invisible handle, so a group heading can occupy a legend row."""
    return Line2D([], [], linestyle="none", marker="none")


def _legend(a, groups, y=-0.165):
    """One legend under the axes, laid out as titled columns.

    Eight entries in a flat strip is a lookup problem: the reader has to scan
    the whole row to find the curve they want. matplotlib fills a multi-column
    legend column-major, so padding every group to the same length puts each
    group in its own column under its own bold heading -- and the headings are
    what carry the provenance the short entry labels leave out.

    `groups` is [(heading, [label, ...]), ...]; labels not actually plotted in
    this panel are skipped, so one group list can serve several panels.
    """
    have = {lbl: h for h, lbl in zip(*a.get_legend_handles_labels())}
    cols = [(t, [k for k in keys if k in have]) for t, keys in groups]
    cols = [c for c in cols if c[1]]
    rows = max(len(keys) for _t, keys in cols)
    handles, labels, headings = [], [], set()
    for title, keys in cols:
        handles.append(_blank()); labels.append(title); headings.add(title)
        for k in keys:
            handles.append(have[k]); labels.append(k)
        handles += [_blank()] * (rows - len(keys))
        labels += [""] * (rows - len(keys))
    leg = a.legend(handles, labels, ncol=len(cols), loc="upper center",
                   bbox_to_anchor=(0.5, y), fontsize=FS_LEG, frameon=False,
                   handlelength=2.2, handletextpad=0.7, columnspacing=1.4,
                   labelspacing=0.45, borderaxespad=0.0, alignment="left")
    for t in leg.get_texts():
        if t.get_text() in headings:
            t.set_fontweight("bold")
    return leg


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
           label="$\\sigma_0(T)$ as coded")
    a.plot(ce._SIG0_T[:9], ce._SIG0_S[:9], "o", ms=7, color=C[0],
           markerfacecolor="white", markeredgewidth=1.8, zorder=5,
           label="table, 9 points in range")
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
    a.set_xlim(-41, 0.5)
    a.set_ylim(2.0e-3, 6.0e-1)
    _style(a, "T [°C]", "$\\sigma_0$  (critical supersaturation) [-]",
           "1.  Libbrecht's data: $\\sigma_0(T)$")
    _legend(a, [("Libbrecht (2017)", ["$\\sigma_0(T)$ as coded",
                                      "table, 9 points in range"])])


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
    for Tm in (-2.0, -20.0, -40.0):
        s0m = ce.sigma0(Tm)
        a.axvline(s0m, lw=1.0, ls=(0, (1, 3)), color=MUTED, zorder=0)
        a.annotate(f"$T$ = {Tm:.0f} °C", (s0m, ALPHA_CEIL * 0.55), rotation=90,
                   fontsize=FS_NOTE, color=MUTED, ha="right", va="top", zorder=6,
                   bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                             alpha=0.85))
    a.annotate("literature band for $\\alpha_c$", (2.9e-4, 3.0e-2),
               fontsize=FS_NOTE, color=C[0], va="center")
    _style(a, "$\\sigma_0(T)$  [-]", "$\\alpha_c$  [-]",
           "2.  $\\alpha_c = \\exp(-\\sigma_0/\\sigma)$ against $\\sigma_0$")
    _legend(a, [(GRP_LAB, SIG_LAB), (GRP_OURS, SIG_OURS),
                (GRP_REF, ["$\\alpha_c = 1$  physical ceiling"])])


def panel_alpha_vs_T(a, T, ctx):
    """3. alpha_c(T) at each sigma -- the supersaturation gap, in T."""
    a.axhspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10, zorder=0)
    a.set_ylim(*VIEW_ALPHA)
    a.annotate("literature band for $\\alpha_c$", (-39.5, 4.0e-3),
               fontsize=FS_NOTE, color=C[0], va="center")
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        y = np.maximum(ctx["alpha"][sig], ALPHA_FLOOR)
        a.plot(T, y, lw=2.4, color=col, ls=_ls(kind), label=lbl)
        _mark_offscale(a, T, y, col, row=i % 2)
    a.axhline(ALPHA_CEIL, lw=1.6, color=INK, alpha=0.7,
              label="$\\alpha_c = 1$  physical ceiling")
    a.axhline(ctx["alpha_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"$\\alpha_c$ = {ctx['alpha_run']:g}  what we run")
    _style(a, "T [°C]", "$\\alpha_c$  [-]",
           "3.  $\\alpha_c(T)$: two decades down in $\\sigma$, "
           "thirty down in $\\alpha_c$")
    _legend(a, [(GRP_LAB, SIG_LAB), (GRP_OURS, SIG_OURS),
                (GRP_REF, ["$\\alpha_c = 1$  physical ceiling",
                           f"$\\alpha_c$ = {ctx['alpha_run']:g}  what we run"])])


def panel_beta_vs_T(a, T, ctx):
    """4. beta_sub(T). The defining equation rides in the title, not the frame,
    so it cannot land on the legend or on a curve."""
    a.axhspan(MF_BETA_LO, MF_BETA_HI, color=C[2], alpha=0.14, zorder=0,
              label="M&F (2024) Table S1")
    a.set_ylim(*VIEW_BETA)
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        a.plot(T, ctx["beta"][sig], lw=2.4, color=col, ls=_ls(kind), label=lbl)
        _mark_offscale(a, T, ctx["beta"][sig], col, row=i % 2)
    a.axhline(ctx["beta_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"$\\beta_{{sub}}$ = {ctx['beta_run']:.1e} s/m  what we run")
    _style(a, "T [°C]", "$\\beta_{sub}$  [s/m]",
           "4.  Attachment resistance $\\beta_{sub} \\propto 1/\\alpha_c$\n\n"
           r"$\beta_{sub}(T)=\dfrac{\beta_{HK}}{\rho_{vs}(T)/\rho_i}"
           r"=\dfrac{\rho_i}{\rho_{vs}(T)}\;\dfrac{1}{\alpha_c}\;"
           r"\sqrt{\dfrac{2\pi m}{k_B T}}\,,\qquad"
           r"\alpha_c=\exp\left[-\sigma_0(T)/\sigma\right]$",
           title_size=13, pad=14)
    _legend(a, [(GRP_LAB, SIG_LAB), (GRP_OURS, SIG_OURS),
                (GRP_REF, ["M&F (2024) Table S1",
                           f"$\\beta_{{sub}}$ = {ctx['beta_run']:.1e} s/m  "
                           f"what we run"])])


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
    a.axhline(ctx["eps_usable"], lw=1.6, color=MUTED,
              label=f"$N_x = 10^4$ ceiling, $\\epsilon$ = "
                    f"{ctx['eps_usable']*1e9:.0f} nm")
    a.axhline(ctx["eps_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"production run, $\\epsilon$ = {ctx['eps_run']*1e6:.3g} µm")
    _style(a, "T [°C]", "$\\epsilon$  [m]",
           "5.  Interface width from the comp_eps.py bounds "
           "(shading follows $\\sigma = 2.7\\times10^{-3}$)")
    _legend(a, [(GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"]]),
                (GRP_LAB, SIG_LAB), (GRP_OURS, SIG_OURS),
                (GRP_REF, [f"$N_x = 10^4$ ceiling, $\\epsilon$ = "
                           f"{ctx['eps_usable']*1e9:.0f} nm",
                           f"production run, $\\epsilon$ = "
                           f"{ctx['eps_run']*1e6:.3g} µm"])])


def panel_Nx_vs_T(a, T, ctx):
    """6. Nx(T) = ceil(sqrt(2)*Lx/eps), against what a machine can hold."""
    prim = ctx["primary_sigma"]
    a.set_ylim(*VIEW_NX)
    _shade_binding(a, T, ctx["binding"][prim], y_frac=0.975)
    for i, (sig, col, lbl, kind) in enumerate(SIGMA_CASES):
        a.plot(T, ctx["Nx"][sig], lw=3.4 if sig == prim else 2.4, color=col,
               ls=_ls(kind), label=lbl)
        _mark_offscale(a, T, ctx["Nx"][sig], col, fmt="{:.0e}", row=i % 2)
    a.axhline(NX_USABLE, lw=1.6, color=MUTED,
              label="$N_x = 10^4$  practical ceiling")
    a.axhline(ctx["Nx_run"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"production mesh, $N_x$ = {ctx['Nx_run']}")
    _style(a, "T [°C]", "$N_x$  [nodes across $L_x$]",
           "6.  The mesh Libbrecht's $\\alpha_c$ demands")
    _legend(a, [(GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"]]),
                (GRP_LAB, SIG_LAB), (GRP_OURS, SIG_OURS),
                (GRP_REF, ["$N_x = 10^4$  practical ceiling",
                           f"production mesh, $N_x$ = {ctx['Nx_run']}"])])


def _shade_binding(a, T, binding, y_frac=0.975):
    """Shade contiguous T-runs sharing a binding bound; name each region once.

    y_frac lets the caller drop the region names clear of whatever else lives
    at the top of that panel -- panel 5 parks its legend there.
    """
    y, va = y_frac, "top"
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
        a.annotate(name, (0.5 * (T[i] + T[j]), y),
                   xycoords=("data", "axes fraction"), fontsize=FS_NOTE + 1,
                   color=BOUND_COLOR.get(name, MUTED), ha="center", va=va,
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
        fig, a = plt.subplots(figsize=FIGSIZE, layout="constrained")
        fn(a, T, ctx)
        fig.savefig(args.out / f"{name}.png", dpi=200)
        plt.close(fig)
        print(f"plot -> {args.out / (name + '.png')}")

    # --- one-slide overview ----------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=OVERVIEW_SIZE, layout="constrained")
    fig.suptitle("Libbrecht (2017) $\\sigma_0(T)$ as the source of $\\alpha_c$: "
                 f"the chain down to $\\epsilon$ and the mesh — {args.label}, "
                 f"measured $v_n$ = {args.vn:.3g} m/s",
                 fontsize=19, color=INK)
    for ax, (_n, fn) in zip(axes.ravel(), PANELS):
        fn(ax, T, ctx)
    fig.savefig(args.out / "fig0_overview.png", dpi=140)
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
