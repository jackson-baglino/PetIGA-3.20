#!/usr/bin/env python3
"""What a CONSTANT alpha_c costs, for the three values the literature supports.

WHAT THIS FIGURE SET IS FOR
---------------------------
Its companion, plot_libbrecht_constraints.py, argues that Libbrecht's
sigma0(T) cannot set alpha_c: no single reference sigma keeps
alpha_c = exp(-sigma0/sigma) inside the literature band across -40..-1 C, and
following it down in temperature demands an interface width below the size of
a water molecule.

This set assumes the conclusion of that one and asks the practical question
instead. If alpha_c is simply a CONSTANT chosen inside the band the literature
supports -- 1e-3 < alpha_c < 1e-1 (Libbrecht 2017; Braun, Fourteau & Lowe,
The Cryosphere 18, 1653, 2024) -- what does each choice cost?

    alpha_c  ->  beta_sub ~ 1/alpha_c  ->  the K&P bounds  ->  eps  ->  Nx

No sigma0 table, no nucleation law, no extrapolation. Just three numbers the
literature already endorses, pushed through comp_eps.compute_eps().

THE RESULT: alpha_c IS NOT MONOTONE IN COST
-------------------------------------------
It is tempting to read "alpha_c is a free parameter, pick a small one and the
kinetics are slow but cheap". The bounds say otherwise, because two of them
move in OPPOSITE directions:

    B-HEAT / B-VAPOR   eps <= s*D*beta_HK        ~ 1/alpha_c   loosens as alpha_c falls
    B-KINETIC          eps <= s*d0/(beta_sub*vn) ~ alpha_c     TIGHTENS as alpha_c falls

so Nx(alpha_c) is V-shaped and there is a genuine optimum, near alpha_c ~ 0.02
to 0.11 over -40..-5 C for the reference case. Both ends are expensive: at
alpha_c = 1e-4 the mesh is 2.6e6 across x at -40 C, and at alpha_c = 1 it is
2.1e4 -- above the practical ceiling either way. Panels 4 and 5 are that
statement; panel 4 shows the mechanism, panel 5 shows where the window is.

The high end has a second cost the mesh does not show. Raising alpha_c shrinks
beta_HK, which shrinks the binding eps, which is what makes the thin-interface
corrections a larger share of tau_sub. Panel 6 plots that share directly:
below f_kin ~ 0.5 comp_eps flags the sharp-interface analogy as stretched.

WHAT THE PLOTS ARE
------------------
  fig1_beta_sub_vs_T.png   beta_sub(T) at each alpha_c, vs M&F's Table S1 range
  fig2_eps_vs_T.png        eps(T), shaded by which bound binds
  fig3_Nx_vs_T.png         Nx(T), vs the practical ceiling and the run we do
  fig4_bounds_vs_alpha.png the four bounds against alpha_c -- why there is a V
  fig5_Nx_vs_alpha.png     Nx against alpha_c at three temperatures
  fig6_fkin_vs_alpha.png   is the sharp-interface mapping still meaningful?
  fig0_master.png          all six, 3 x 2, one shared legend
  fig7_legend.png          that legend alone
  alpha_c_sizing.csv       the numbers behind them

Layout, type scale and the 10 in width ceiling come from figstyle.py, shared
with the Libbrecht set so the two cannot drift.

Usage:
    python preprocess/plot_alpha_constant_constraints.py \
        --out studies/alpha_c_sizing
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comp_eps as ce
import figstyle as fs
from figstyle import (C, INK, MUTED, GRID, BOUND_COLOR, BOUND_LABEL,
                      GRP_BIND, GRP_REF, FS_NOTE, FS_C_TICK, SLIDE_W_IN)
from figstyle import (style as _style, legend as _legend, compact as _compact,
                      lbl as _lbl, mark_offscale as _mark_offscale,
                      shade_binding as _shade_binding)

# The three values, and only these. Same visual language as the Libbrecht set:
# green solid, blue dashed, pink dash-dot, in descending order of alpha_c.
ALPHA_CASES = [
    (1.0e-1, C[2], "$\\alpha_c = 10^{-1}$", "-"),
    (1.0e-2, C[5], "$\\alpha_c = 10^{-2}$", (0, (5, 2))),
    (1.0e-3, C[3], "$\\alpha_c = 10^{-3}$", (0, (5, 1.5, 1, 1.5))),
]
GRP_ALPHA = "Literature band"
ALPHA_LABELS = [l for _v, _c, l, _s in ALPHA_CASES]

# Temperatures the alpha_c sweeps (panels 4-6) are evaluated at. Temperature is
# an ORDERED variable, so it gets a sequential ramp rather than three more
# categorical hues -- which also keeps it clear of the bound colours, since the
# master legend carries both families at once.
T_RAMP = ["#3B0F70", "#B63679", "#FE9F6D"]      # magma, cold -> warm
SWEEP_T = [(-40.0, T_RAMP[0], "$T = -40$ °C", "-"),
           (-20.0, T_RAMP[1], "$T = -20$ °C", (0, (5, 2))),
           (-5.0, T_RAMP[2], "$T = -5$ °C", (0, (1, 1.6)))]
GRP_SWEEP = "Temperature"
SWEEP_LABELS = [l for _v, _c, l, _s in SWEEP_T]

# a = (m/rho_ice)^(1/3): the molecular length implied by comp_eps's constants.
A_MOL = (ce._M_H2O / ce._RHO_ICE) ** (1.0 / 3.0)

MF_BETA_LO, MF_BETA_HI = 2.0e4, 2.0e6      # M&F (2024) Table S1
NX_USABLE = 1.0e4                          # the practical ceiling

# Display windows, sized to the curves rather than to sixty decades of nothing.
VIEW_BETA = (1.0e3, 2.0e8)
VIEW_NX_T = (3.0e2, 5.0e5)                 # panels 2-3 (eps is the reciprocal)
VIEW_NX_A = (3.0e2, 5.0e6)                 # panel 5, which spans more alpha_c
VIEW_BOUND = (1.0e-9, 1.0e-3)              # panel 4
ALPHA_SWEEP = np.logspace(-4.0, 0.0, 220)  # panel 4-6 x-axis


def _ls(style):
    return style


# =========================================================================
# Panels
# =========================================================================

def _alpha_band(a, horizontal=True):
    """Shade the band the literature supports. On panels 4-6 that is a band in
    x, and it is the only region a choice is allowed to come from."""
    span = a.axvspan if not horizontal else a.axhspan
    span(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10, zorder=0,
         label="literature band, $10^{-3}$–$10^{-1}$")


def panel_beta(a, T, ctx):
    """1. beta_sub(T) = (rho_i/rho_vs)(1/alpha_c) sqrt(2 pi m / kT)."""
    a.axhspan(MF_BETA_LO, MF_BETA_HI, color=C[2], alpha=0.14, zorder=0,
              label="M&F (2024) Table S1")
    a.set_ylim(*VIEW_BETA)
    for i, (al, col, lbl, st) in enumerate(ALPHA_CASES):
        a.plot(T, ctx["beta"][al], lw=2.4, color=col, ls=_ls(st), label=lbl)
        _mark_offscale(a, T, ctx["beta"][al], col, row=i % 2)
    run = _lbl(f"$\\beta_{{sub}}$ = {ctx['beta_run']:.1e} s/m  what we run",
               "what we run")
    a.axhline(ctx["beta_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "T [°C]", "$\\beta_{sub}$  [s/m]",
           "1.  Attachment resistance $\\beta_{sub} \\propto 1/\\alpha_c$",
           short="1.  $\\beta_{sub}(T)$")
    _legend(a, [(GRP_ALPHA, ALPHA_LABELS),
                (GRP_REF, ["M&F (2024) Table S1", run])])


def panel_eps(a, T, ctx):
    """2. eps(T) straight out of comp_eps.compute_eps."""
    prim = ctx["primary_alpha"]
    a.set_ylim(*ctx["view_eps"])
    _shade_binding(a, T, ctx["binding"][prim], y_frac=0.975)
    for i, (al, col, lbl, st) in enumerate(ALPHA_CASES):
        a.plot(T, ctx["eps"][al], lw=3.4 if al == prim else 2.4, color=col,
               ls=_ls(st), label=lbl)
        _mark_offscale(a, T, ctx["eps"][al], col, row=i % 2)
    ceil = _lbl(f"$N_x = 10^4$ ceiling, $\\epsilon$ = "
                f"{ctx['eps_usable']*1e9:.0f} nm", "$N_x = 10^4$  practical ceiling")
    a.axhline(ctx["eps_usable"], lw=1.6, color=MUTED, label=ceil)
    run = _lbl(f"production run, $\\epsilon$ = {ctx['eps_run']*1e6:.3g} µm",
               "what we run")
    a.axhline(ctx["eps_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "T [°C]", "$\\epsilon$  [m]",
           "2.  Interface width, shaded by which bound binds at "
           "$\\alpha_c = 10^{-2}$",
           short="2.  Interface width $\\epsilon(T)$")
    _legend(a, [(GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"]]),
                (GRP_ALPHA, ALPHA_LABELS), (GRP_REF, [ceil, run])])


def panel_Nx(a, T, ctx):
    """3. Nx(T) = ceil(sqrt(2)*Lx/eps)."""
    prim = ctx["primary_alpha"]
    a.set_ylim(*VIEW_NX_T)
    _shade_binding(a, T, ctx["binding"][prim], y_frac=0.03)
    for i, (al, col, lbl, st) in enumerate(ALPHA_CASES):
        a.plot(T, ctx["Nx"][al], lw=3.4 if al == prim else 2.4, color=col,
               ls=_ls(st), label=lbl)
        _mark_offscale(a, T, ctx["Nx"][al], col, fmt="{:.0e}", row=i % 2)
    a.axhline(NX_USABLE, lw=1.6, color=MUTED,
              label="$N_x = 10^4$  practical ceiling")
    run = _lbl(f"production mesh, $N_x$ = {ctx['Nx_run']}", "what we run")
    a.axhline(ctx["Nx_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "T [°C]", "$N_x$  [nodes across $L_x$]",
           "3.  The mesh each constant $\\alpha_c$ demands",
           short="3.  Mesh $N_x(T)$", short_y="$N_x$")
    _legend(a, [(GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"]]),
                (GRP_ALPHA, ALPHA_LABELS),
                (GRP_REF, ["$N_x = 10^4$  practical ceiling", run])])


def panel_bounds(a, T, ctx):
    """4. The four bounds against alpha_c. Two of them move in opposite
    directions, which is the whole reason there is an optimum."""
    Tb = ctx["T_bounds"]
    a.set_xscale("log")
    a.set_ylim(*VIEW_BOUND)
    _alpha_band(a, horizontal=False)
    for key in ("B-HEAT", "B-VAPOR", "B-KINETIC"):
        a.plot(ALPHA_SWEEP, ctx["bounds"][key], lw=2.2, color=BOUND_COLOR[key],
               ls="-" if key == "B-KINETIC" else (0, (5, 2)),
               label=BOUND_LABEL[key])
    a.axhline(ctx["bounds"]["B-CURV"], lw=1.8, ls=(0, (1, 2)),
              color=BOUND_COLOR["B-CURV"], label=BOUND_LABEL["B-CURV"])
    a.plot(ALPHA_SWEEP, ctx["bounds"]["eps"], lw=3.4, color=INK,
           label="$\\epsilon$ = the smallest of them")
    _style(a, "$\\alpha_c$  [-]", "bound on $\\epsilon$  [m]",
           f"4.  The four bounds against $\\alpha_c$ at $T$ = {Tb:.0f} °C",
           short=f"4.  Bounds vs $\\alpha_c$ ({Tb:.0f} °C)",
           short_y="$\\epsilon$ bound  [m]")
    _legend(a, [("K&P bounds", [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"],
                                BOUND_LABEL["B-VAPOR"], BOUND_LABEL["B-CURV"]]),
                (GRP_REF, ["$\\epsilon$ = the smallest of them",
                           "literature band, $10^{-3}$–$10^{-1}$"])])


def panel_Nx_alpha(a, T, ctx):
    """5. Nx against alpha_c: the design chart. Read the mesh off the axis."""
    a.set_xscale("log")
    a.set_ylim(*VIEW_NX_A)
    _alpha_band(a, horizontal=False)
    for i, (t, col, lbl, st) in enumerate(SWEEP_T):
        y = ctx["Nx_sweep"][t]
        a.plot(ALPHA_SWEEP, y, lw=2.4, color=col, ls=_ls(st), label=lbl)
        _mark_offscale(a, ALPHA_SWEEP, y, col, fmt="{:.0e}", row=i % 2)
        j = int(np.argmin(y))
        a.plot([ALPHA_SWEEP[j]], [y[j]], "v", ms=8, color=col, zorder=6)
    a.axhline(NX_USABLE, lw=1.6, color=MUTED,
              label="$N_x = 10^4$  practical ceiling")
    if not _compact():
        a.annotate("▼ the cheapest $\\alpha_c$ at each $T$", (2.0e-2, 6.0e2),
                   fontsize=FS_NOTE, color=INK, ha="center")
    _style(a, "$\\alpha_c$  [-]", "$N_x$  [nodes across $L_x$]",
           "5.  Mesh against $\\alpha_c$: both ends are expensive",
           short="5.  $N_x$ vs $\\alpha_c$", short_y="$N_x$")
    _legend(a, [(GRP_SWEEP, SWEEP_LABELS),
                (GRP_REF, ["$N_x = 10^4$  practical ceiling",
                           "literature band, $10^{-3}$–$10^{-1}$"])])


def panel_fkin(a, T, ctx):
    """6. The kinetic share of the tau_sub bracket -- whether the
    sharp-interface mapping still means anything at this alpha_c."""
    a.set_xscale("log")
    a.set_ylim(0.0, 1.08)
    _alpha_band(a, horizontal=False)
    for i, (t, col, lbl, st) in enumerate(SWEEP_T):
        a.plot(ALPHA_SWEEP, ctx["fkin_sweep"][t], lw=2.4, color=col,
               ls=_ls(st), label=lbl)
    a.axhline(0.5, lw=1.8, ls=(0, (1, 2)), color=INK,
              label="$f_{kin} = 0.5$  comp_eps flags below this")
    _style(a, "$\\alpha_c$  [-]",
           "$f_{kin}$  (kinetic share of $\\tau_{sub}$)",
           "6.  Is the sharp-interface mapping still meaningful?",
           logy=False, short="6.  Validity $f_{kin}(\\alpha_c)$",
           short_y="$f_{kin}$")
    _legend(a, [(GRP_SWEEP, SWEEP_LABELS),
                (GRP_REF, ["$f_{kin} = 0.5$  comp_eps flags below this",
                           "literature band, $10^{-3}$–$10^{-1}$"])])


PANELS = [
    ("fig1_beta_sub_vs_T",   panel_beta),
    ("fig2_eps_vs_T",        panel_eps),
    ("fig3_Nx_vs_T",         panel_Nx),
    ("fig4_bounds_vs_alpha", panel_bounds),
    ("fig5_Nx_vs_alpha",     panel_Nx_alpha),
    ("fig6_fkin_vs_alpha",   panel_fkin),
]


def master_groups(ctx):
    return [
        (GRP_ALPHA, ALPHA_LABELS),
        (GRP_SWEEP, SWEEP_LABELS),
        (GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"],
                    BOUND_LABEL["B-VAPOR"], BOUND_LABEL["B-CURV"]]),
        (GRP_REF, ["$\\epsilon$ = the smallest of them",
                   "$N_x = 10^4$  practical ceiling",
                   "what we run", "M&F (2024) Table S1"]),
    ]


# =========================================================================

def build_context(T, args):
    ctx = dict(primary_alpha=args.primary_alpha, Lx=args.Lx, Ly=args.Ly,
               T_bounds=args.T_bounds)
    kw = dict(Lx=args.Lx, Ly=args.Ly, Lz=0.0, Rave=args.Rave, v_n=args.vn,
              safety=args.safety, eps_over_R=args.eps_over_R)

    for k in ("beta", "eps", "Nx", "nodes", "binding", "fkin"):
        ctx[k] = {}
    for al, _c, _l, _s in ALPHA_CASES:
        r = [ce.compute_eps(T0_C=t, alpha_c=al, **kw) for t in T]
        ctx["beta"][al] = np.array([x["beta_uns"] for x in r])
        ctx["eps"][al] = np.array([x["eps"] for x in r])
        ctx["Nx"][al] = np.array([float(x["Nx"]) for x in r])
        ctx["nodes"][al] = np.array([float(x["Nx"]) * float(x["Ny"]) for x in r])
        ctx["binding"][al] = [x["binding"] for x in r]
        ctx["fkin"][al] = np.array([x["kinetic_frac"] for x in r])

    # The alpha_c sweeps, at a few fixed temperatures.
    ctx["Nx_sweep"], ctx["fkin_sweep"] = {}, {}
    for t, _c, _l, _s in SWEEP_T:
        r = [ce.compute_eps(T0_C=t, alpha_c=a, **kw) for a in ALPHA_SWEEP]
        ctx["Nx_sweep"][t] = np.array([float(x["Nx"]) for x in r])
        ctx["fkin_sweep"][t] = np.array([x["kinetic_frac"] for x in r])

    # The four bounds themselves, at one temperature, as alpha_c varies.
    Tb = args.T_bounds
    bhk = np.array([ce.beta_HK(Tb, a) for a in ALPHA_SWEEP])
    rho_rat = ce.rho_vs_sat(Tb) / ce._RHO_ICE
    d0, Dv = ce.capillary_length(Tb), ce.Dv_T(Tb)
    s = args.safety
    ctx["bounds"] = {
        "B-HEAT":  s * ce.D_heat_of("mean") * bhk,
        "B-VAPOR": s * Dv * bhk,
        "B-KINETIC": s * d0 / ((bhk / rho_rat) * args.vn),
        "B-CURV":  args.eps_over_R * args.Rave,
    }
    ctx["bounds"]["eps"] = np.minimum(
        np.minimum(ctx["bounds"]["B-HEAT"], ctx["bounds"]["B-VAPOR"]),
        np.minimum(ctx["bounds"]["B-KINETIC"], ctx["bounds"]["B-CURV"]))

    run = ce.compute_eps(T0_C=args.T_run, alpha_c=args.alpha_run, **kw)
    ctx["beta_run"] = run["beta_uns"]
    ctx["eps_run"] = args.eps_run if args.eps_run > 0 else run["eps"]
    ctx["Nx_run"] = math.ceil(math.sqrt(2.0) * args.Lx / ctx["eps_run"])
    ctx["nodes_run"] = ctx["Nx_run"] * math.ceil(math.sqrt(2.0) * args.Ly
                                                 / ctx["eps_run"])
    # eps window is the exact reciprocal of the Nx one, so panels 2 and 3 can
    # be read against each other.
    ctx["view_eps"] = (math.sqrt(2.0) * args.Lx / VIEW_NX_T[1],
                       math.sqrt(2.0) * args.Lx / VIEW_NX_T[0])
    ctx["eps_usable"] = math.sqrt(2.0) * args.Lx / NX_USABLE
    return ctx


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Reference-case geometry: the Molaro dom2 production domain. Lx enters
    # only as Nx = ceil(sqrt(2)*Lx/eps), i.e. LINEARLY.
    p.add_argument("--Lx", type=float, default=4.50e-4, help="domain x [m]")
    p.add_argument("--Ly", type=float, default=2.25e-4, help="domain y [m]")
    p.add_argument("--Rave", type=float, default=8.675e-5,
                   help="representative grain radius [m], for B-CURV")
    p.add_argument("--vn", type=float, default=3.415598e-9,
                   help="MEASURED front velocity [m/s]; K&P Eq. 45 is only a "
                        "constraint if the front actually moves")
    p.add_argument("--safety", type=float, default=0.5)
    p.add_argument("--eps_over_R", type=float, default=0.05)
    p.add_argument("--primary_alpha", type=float, default=1.0e-2,
                   help="the alpha_c whose binding bound is shaded in figs 2-3")
    p.add_argument("--T_bounds", type=float, default=-20.0,
                   help="temperature the alpha_c sweep in fig 4 is taken at")
    p.add_argument("--alpha_run", type=float, default=0.1)
    p.add_argument("--T_run", type=float, default=-20.0)
    p.add_argument("--eps_run", type=float, default=1.18e-7)
    p.add_argument("--Tmin", type=float, default=-40.0)
    p.add_argument("--Tmax", type=float, default=-1.0)
    p.add_argument("--nT", type=int, default=400)
    p.add_argument("--width", type=float, default=SLIDE_W_IN)
    p.add_argument("--panel_width", type=float, default=3.0)
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--panel_legend", choices=("separate", "inline", "none"),
                   default="separate")
    p.add_argument("--label", default="reference case: 450 × 225 µm grain pair")
    p.add_argument("--out", type=Path, default=Path("studies/alpha_c_sizing"))
    args = p.parse_args()

    T = np.linspace(args.Tmin, args.Tmax, args.nT)
    ctx = build_context(T, args)
    args.out.mkdir(parents=True, exist_ok=True)

    fs.emit(args.out, PANELS, T, ctx, width=args.width,
            panel_width=args.panel_width, dpi=args.dpi,
            panel_legend=args.panel_legend, master_groups=master_groups,
            suptitle="A constant $\\alpha_c$ inside the literature band: what "
                     "each choice costs in $\\beta_{sub}$, $\\epsilon$ and mesh\n"
                     f"{args.label}, measured $v_n$ = {args.vn:.3g} m/s")

    # --- the numbers ------------------------------------------------------
    csv = args.out / "alpha_c_sizing.csv"
    with open(csv, "w") as fh:
        fh.write("T_C,alpha_c,beta_sub_s_per_m,eps_m,binding,Nx,nodes_2D,f_kin\n")
        for i in range(0, len(T), max(1, len(T) // 80)):
            for al, _c, _l, _s in ALPHA_CASES:
                fh.write(f"{T[i]:.3f},{al:.3e},{ctx['beta'][al][i]:.6e},"
                         f"{ctx['eps'][al][i]:.6e},{ctx['binding'][al][i]},"
                         f"{ctx['Nx'][al][i]:.6e},{ctx['nodes'][al][i]:.6e},"
                         f"{ctx['fkin'][al][i]:.4f}\n")
    print(f"csv  -> {csv}")

    print(f"\n  Production run:  alpha_c = {args.alpha_run:g}, "
          f"beta_sub = {ctx['beta_run']:.3e} s/m, eps = {ctx['eps_run']:.3e} m, "
          f"Nx = {ctx['Nx_run']}, {ctx['nodes_run']/1e6:.1f} M nodes")
    print(f"  Practical ceiling Nx = {NX_USABLE:.0g} is eps = "
          f"{ctx['eps_usable']*1e9:.0f} nm\n")
    print(f"  {'T':>5} {'alpha_c':>9} {'beta_sub':>11} {'eps [m]':>11} "
          f"{'bind':>10} {'Nx':>9} {'nodes':>11} {'f_kin':>6}")
    for t in (-40.0, -20.0, -5.0):
        i = int(np.argmin(abs(T - t)))
        for al, _c, _l, _s in ALPHA_CASES:
            print(f"  {T[i]:>5.1f} {al:>9.1e} {ctx['beta'][al][i]:>11.3e} "
                  f"{ctx['eps'][al][i]:>11.3e} {ctx['binding'][al][i]:>10s} "
                  f"{ctx['Nx'][al][i]:>9.0f} {ctx['nodes'][al][i]:>11.3e} "
                  f"{ctx['fkin'][al][i]:>6.2f}")
        print()
    print("  Cheapest alpha_c, and what it costs:")
    for t, _c, _l, _s in SWEEP_T:
        y = ctx["Nx_sweep"][t]
        j = int(np.argmin(y))
        print(f"  T = {t:6.1f} C : alpha_c = {ALPHA_SWEEP[j]:.3g} -> "
              f"Nx = {y[j]:.0f}   (Nx = {y[0]:.3g} at 1e-4, {y[-1]:.3g} at 1)")


if __name__ == "__main__":
    main()
