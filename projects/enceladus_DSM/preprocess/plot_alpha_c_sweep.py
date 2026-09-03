#!/usr/bin/env python3
"""alpha_c as a knob: sweep it over [1e-4, 1] and read off what each value costs.

WHAT THIS FIGURE SET IS FOR
---------------------------
Its companion, plot_libbrecht_constraints.py, argues that Libbrecht's
sigma0(T) cannot set alpha_c. This one takes alpha_c as what it is in the
solver -- a constant we choose -- and puts it on the x-axis, so every panel
answers "if I pick THIS alpha_c, what do I get?".

    alpha_c  ->  beta_sub ~ 1/alpha_c  ->  the K&P bounds  ->  eps  ->  Nx

The band the literature supports, 1e-3 < alpha_c < 1e-1 (Libbrecht 2017;
Braun, Fourteau & Lowe, The Cryosphere 18, 1653, 2024), is shaded on every
panel: it is the only region a choice is allowed to come from. The sweep runs
a decade past it on each side so the shape of the trade is visible.

THE RESULT: alpha_c IS NOT MONOTONE IN COST
-------------------------------------------
It is tempting to read "alpha_c is free, so pick a small one -- the kinetics
are slow but the mesh is cheap". The bounds say otherwise, because two of them
move in OPPOSITE directions:

    B-HEAT / B-VAPOR   eps <= s*D*beta_HK        ~ 1/alpha_c   loosens as alpha_c falls
    B-KINETIC          eps <= s*d0/(beta_sub*vn) ~ alpha_c     TIGHTENS as alpha_c falls

so eps(alpha_c) is a tent and Nx(alpha_c) is a V, with a genuine optimum near
alpha_c = 0.02 (at -5 C) to 0.115 (at -40 C). Both ends are expensive: at
alpha_c = 1e-4 the mesh is 2.6e6 across x at -40 C, and at alpha_c = 1 it is
2.1e4. Panel 2 is the mechanism, panels 3 and 4 are the consequence.

WHY beta_sub CARRIES A TEMPERATURE DEPENDENCE AT ALL
----------------------------------------------------
Worth being explicit, because there are two coefficients and only one of them
really moves with T:

    beta_HK(T)  = (1/alpha_c) sqrt(2 pi m / k_B T)          SCALED   (K&P beta')
    beta_sub(T) = beta_HK / (rho_vs(T)/rho_i)               UNSCALED (K&P beta_0)

beta_HK is the physical Hertz-Knudsen coefficient and its only T-dependence is
the sqrt(T) in the mean thermal speed -- worth 7.4% over -40..-1 C, which is
nothing. All of the spread in beta_sub is the rho_vs(T)/rho_i factor that
converts between the two, and rho_vs falls 44x over that same range. So panel
1's three curves are not three different kinetics; they are one kinetics seen
through a saturation vapour density that collapses as it gets cold. The solver
is passed -beta_sub0 (unscaled) and rescales internally, which is why it is
the one plotted.

WHAT THE PLOTS ARE
------------------
  fig1_beta_vs_alpha.png   beta_sub and beta_HK against alpha_c
  fig2_bounds_vs_alpha.png the four K&P bounds -- why there is a tent
  fig3_eps_vs_alpha.png    the interface width that survives
  fig4_Nx_vs_alpha.png     the mesh it demands
  fig5_fkin_vs_alpha.png   is the sharp-interface mapping still meaningful?
  fig6_Lstar_vs_alpha.png  when does alpha_c stop mattering at all?
  fig0_master.png          all six, 3 x 2, one shared legend
  fig7_legend.png          that legend alone
  alpha_c_sweep.csv        the numbers behind them

Layout, type scale and the 10 in width ceiling come from figstyle.py, shared
with the Libbrecht set so the two cannot drift.

Usage:
    python preprocess/plot_alpha_c_sweep.py --out studies/alpha_c_sizing
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

# Temperature is an ORDERED variable, so it gets a sequential ramp rather than
# categorical hues -- which also keeps it clear of the bound colours, since the
# master legend carries both families at once.
T_RAMP = ["#3B0F70", "#B63679", "#FE9F6D"]          # magma, cold -> warm
SWEEP_T = [(-40.0, T_RAMP[0], "$T = -40$ °C", "-"),
           (-20.0, T_RAMP[1], "$T = -20$ °C", (0, (5, 2))),
           (-5.0, T_RAMP[2], "$T = -5$ °C", (0, (1, 1.6)))]
GRP_SWEEP = "Temperature"
SWEEP_LABELS = [l for _v, _c, l, _s in SWEEP_T]

LIT_LABEL = "literature band, $10^{-3}$–$10^{-1}$"
MF_BETA_LO, MF_BETA_HI = 2.0e4, 2.0e6               # M&F (2024) Table S1
NX_USABLE = 1.0e4                                   # the practical ceiling

# One decade past the literature band on each side, so the shape of the trade
# is visible without inviting a choice from out there.
ALPHA = np.logspace(-4.0, 0.0, 300)

VIEW_BETA = (1.0e-3, 1.0e9)   # wide enough to hold beta_HK AND beta_sub
VIEW_BOUND = (1.0e-9, 1.0e-3)
VIEW_EPS = (1.0e-10, 3.0e-6)
VIEW_NX = (3.0e2, 5.0e6)
VIEW_LSTAR = (5.0e-8, 5.0e-3)


def _band(a):
    """Shade the band the literature supports. Every panel carries it: it is
    the only region a choice of alpha_c is allowed to come from."""
    a.axvspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10,
              zorder=0, label=LIT_LABEL)


def _mark_binding(a, binding, y_frac, color_by=BOUND_COLOR):
    """Where the binding bound changes, and what binds on each side.

    Deliberately NOT a filled span: the literature band is already a
    translucent fill across these axes, and a second one on top of it mixes
    into a grey that means neither. A rule plus the two names carries the same
    information without touching the band.
    """
    va = "top" if y_frac > 0.5 else "bottom"
    i = 0
    while i < len(ALPHA):
        j = i
        while j + 1 < len(ALPHA) and binding[j + 1] == binding[i]:
            j += 1
        if j + 1 < len(ALPHA):
            a.axvline(ALPHA[j], lw=1.4, ls=(0, (3, 2)), color=MUTED, zorder=1)
        xm = 10.0 ** (0.5 * (math.log10(ALPHA[i]) + math.log10(ALPHA[j])))
        a.annotate(binding[i], (xm, y_frac), xycoords=("data", "axes fraction"),
                   fontsize=FS_C_TICK if _compact() else FS_NOTE + 1,
                   color=color_by.get(binding[i], MUTED), ha="center", va=va,
                   fontweight="bold", zorder=7,
                   bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                             alpha=0.85))
        i = j + 1


def _logx(a):
    a.set_xscale("log")
    a.set_xlim(ALPHA[0], ALPHA[-1])


# =========================================================================
# Panels
# =========================================================================

def panel_beta(a, T, ctx):
    """1. beta_sub and beta_HK against alpha_c.

    Both are drawn because the difference between them is the answer to "why
    does beta_sub depend on temperature at all": beta_HK barely does, and the
    whole spread is the rho_vs(T)/rho_i factor between the two.
    """
    a.axhspan(MF_BETA_LO, MF_BETA_HI, color=C[2], alpha=0.14, zorder=0,
              label="M&F (2024) Table S1")
    _band(a)
    _logx(a)
    a.set_ylim(*VIEW_BETA)
    # Drawn as a line, not a band: its whole -40..-1 C spread is 8%, which is
    # sub-pixel on an axis this tall -- and that invisibility IS the point.
    a.plot(ALPHA, 0.5 * (ctx["bhk_lo"] + ctx["bhk_hi"]), lw=2.2, color=MUTED,
           label="$\\beta_{HK}$ (scaled) — 8% wide over $-40$…$-1$ °C")
    for i, (t, col, lbl, st) in enumerate(SWEEP_T):
        a.plot(ALPHA, ctx["beta"][t], lw=2.4, color=col, ls=st, label=lbl)
        _mark_offscale(a, ALPHA, ctx["beta"][t], col, row=i % 2)
    run = _lbl(f"$\\beta_{{sub}}$ = {ctx['beta_run']:.1e} s/m  what we run",
               "what we run")
    a.axhline(ctx["beta_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    if not _compact():
        a.annotate("the whole spread is $\\rho_{vs}(T)/\\rho_i$:\n"
                   "$\\rho_{vs}$ falls 44× over this range,\n"
                   "$\\beta_{HK}$ moves 7%",
                   (1.3e-4, 1.0e2), fontsize=FS_NOTE, color=INK, va="bottom")
    _style(a, "$\\alpha_c$  [-]", "$\\beta$  [s/m]",
           "1.  $\\beta_{sub} \\propto 1/\\alpha_c$, and where its $T$"
           " dependence comes from",
           short="1.  $\\beta_{sub}(\\alpha_c)$", short_y="$\\beta$  [s/m]")
    _legend(a, [(GRP_SWEEP, SWEEP_LABELS),
                (GRP_REF, ["$\\beta_{HK}$ (scaled) — 8% wide over $-40$…$-1$ °C",
                           "M&F (2024) Table S1", run, LIT_LABEL])])


def panel_bounds(a, T, ctx):
    """2. The four K&P bounds. Two move in opposite directions; that crossing
    is the entire reason there is an optimum."""
    Tb = ctx["T_bounds"]
    _band(a)
    _logx(a)
    a.set_ylim(*VIEW_BOUND)
    for key in ("B-HEAT", "B-VAPOR", "B-KINETIC"):
        a.plot(ALPHA, ctx["bounds"][key], lw=2.2, color=BOUND_COLOR[key],
               ls="-" if key == "B-KINETIC" else (0, (5, 2)),
               label=BOUND_LABEL[key])
    a.axhline(ctx["bounds"]["B-CURV"], lw=1.8, ls=(0, (1, 2)),
              color=BOUND_COLOR["B-CURV"], label=BOUND_LABEL["B-CURV"])
    a.plot(ALPHA, ctx["bounds"]["eps"], lw=3.4, color=INK,
           label="$\\epsilon$ = the smallest of them")
    _style(a, "$\\alpha_c$  [-]", "bound on $\\epsilon$  [m]",
           f"2.  The four bounds at $T$ = {Tb:.0f} °C: B-KINETIC rises, "
           f"B-HEAT falls", short=f"2.  Bounds ({Tb:.0f} °C)",
           short_y="$\\epsilon$ bound  [m]")
    _legend(a, [("K&P bounds", [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"],
                                BOUND_LABEL["B-VAPOR"], BOUND_LABEL["B-CURV"]]),
                (GRP_REF, ["$\\epsilon$ = the smallest of them", LIT_LABEL])])


def panel_eps(a, T, ctx):
    """3. The interface width that survives all four bounds."""
    prim = ctx["primary_T"]
    _logx(a)
    a.set_ylim(*VIEW_EPS)
    _band(a)
    _mark_binding(a, ctx["binding"][prim], 0.03)
    for i, (t, col, lbl, st) in enumerate(SWEEP_T):
        a.plot(ALPHA, ctx["eps"][t], lw=3.4 if t == prim else 2.4, color=col,
               ls=st, label=lbl)
        _mark_offscale(a, ALPHA, ctx["eps"][t], col, row=i % 2)
    ceil = _lbl(f"$N_x = 10^4$ ceiling, $\\epsilon$ = "
                f"{ctx['eps_usable']*1e9:.0f} nm",
                "$N_x = 10^4$  practical ceiling")
    a.axhline(ctx["eps_usable"], lw=1.6, color=MUTED, label=ceil)
    run = _lbl(f"production run, $\\epsilon$ = {ctx['eps_run']*1e6:.3g} µm",
               "what we run")
    a.axhline(ctx["eps_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    _style(a, "$\\alpha_c$  [-]", "$\\epsilon$  [m]",
           f"3.  Interface width (shading: which bound binds at "
           f"{prim:.0f} °C)", short="3.  Interface width $\\epsilon$")
    _legend(a, [(GRP_SWEEP, SWEEP_LABELS),
                (GRP_REF, [ceil, run, LIT_LABEL])])


def panel_Nx(a, T, ctx):
    """4. The mesh that follows, Nx = ceil(sqrt(2)*Lx/eps)."""
    prim = ctx["primary_T"]
    _logx(a)
    a.set_ylim(*VIEW_NX)
    _band(a)
    _mark_binding(a, ctx["binding"][prim], 0.975)
    for i, (t, col, lbl, st) in enumerate(SWEEP_T):
        y = ctx["Nx"][t]
        a.plot(ALPHA, y, lw=3.4 if t == prim else 2.4, color=col, ls=st,
               label=lbl)
        _mark_offscale(a, ALPHA, y, col, fmt="{:.0e}", row=i % 2)
        j = int(np.argmin(y))
        a.plot([ALPHA[j]], [y[j]], "v", ms=8, color=col, zorder=6)
    a.axhline(NX_USABLE, lw=1.6, color=MUTED,
              label="$N_x = 10^4$  practical ceiling")
    run = _lbl(f"production mesh, $N_x$ = {ctx['Nx_run']}", "what we run")
    a.axhline(ctx["Nx_run"], lw=1.8, ls=(0, (1, 2)), color=INK, label=run)
    if not _compact():
        a.annotate("▼ cheapest $\\alpha_c$", (4.0e-2, 6.0e2), fontsize=FS_NOTE,
                   color=INK, ha="center")
    _style(a, "$\\alpha_c$  [-]", "$N_x$  [nodes across $L_x$]",
           "4.  The mesh it demands: both ends are expensive",
           short="4.  Mesh $N_x$", short_y="$N_x$")
    _legend(a, [(GRP_SWEEP, SWEEP_LABELS),
                (GRP_REF, ["$N_x = 10^4$  practical ceiling", run, LIT_LABEL])])


def panel_fkin(a, T, ctx):
    """5. The kinetic share of the tau_sub bracket -- whether the
    sharp-interface mapping still means anything at this alpha_c."""
    _band(a)
    _logx(a)
    a.set_ylim(0.0, 1.08)
    for t, col, lbl, st in SWEEP_T:
        a.plot(ALPHA, ctx["fkin"][t], lw=2.4, color=col, ls=st, label=lbl)
    a.axhline(0.5, lw=1.8, ls=(0, (1, 2)), color=INK,
              label="$f_{kin} = 0.5$  comp_eps flags below this")
    _style(a, "$\\alpha_c$  [-]", "$f_{kin}$  (kinetic share of $\\tau_{sub}$)",
           "5.  Is the sharp-interface mapping still meaningful?", logy=False,
           short="5.  Validity $f_{kin}$", short_y="$f_{kin}$")
    _legend(a, [(GRP_SWEEP, SWEEP_LABELS),
                (GRP_REF, ["$f_{kin} = 0.5$  comp_eps flags below this",
                           LIT_LABEL])])


def panel_Lstar(a, T, ctx):
    """6. L* = beta_HK*D_v, the attachment/diffusion crossover length.

    Features smaller than L* are attachment-limited, so alpha_c controls them;
    features larger are vapour-diffusion-limited, and alpha_c stops mattering.
    That is why raising alpha_c eventually buys nothing.
    """
    _band(a)
    _logx(a)
    a.set_ylim(*VIEW_LSTAR)
    for t, col, lbl, st in SWEEP_T:
        a.plot(ALPHA, ctx["Lstar"][t], lw=2.4, color=col, ls=st, label=lbl)
    a.axhline(ctx["Rave"], lw=1.6, color=MUTED,
              label=f"grain radius {ctx['Rave']*1e6:.0f} µm")
    a.axhline(ctx["rneck"], lw=1.8, ls=(0, (1, 2)), color=INK,
              label=f"neck radius {ctx['rneck']*1e6:.0f} µm")
    _style(a, "$\\alpha_c$  [-]", "$L^* = \\beta_{HK} D_v$  [m]",
           "6.  Above $L^*$ transport limits, and $\\alpha_c$ stops mattering",
           short="6.  Regime $L^*(\\alpha_c)$", short_y="$L^*$  [m]")
    _legend(a, [(GRP_SWEEP, SWEEP_LABELS),
                (GRP_REF, [f"grain radius {ctx['Rave']*1e6:.0f} µm",
                           f"neck radius {ctx['rneck']*1e6:.0f} µm",
                           LIT_LABEL])])


PANELS = [
    ("fig1_beta_vs_alpha",   panel_beta),
    ("fig2_bounds_vs_alpha", panel_bounds),
    ("fig3_eps_vs_alpha",    panel_eps),
    ("fig4_Nx_vs_alpha",     panel_Nx),
    ("fig5_fkin_vs_alpha",   panel_fkin),
    ("fig6_Lstar_vs_alpha",  panel_Lstar),
]


def master_groups(ctx):
    return [
        (GRP_SWEEP, SWEEP_LABELS),
        (GRP_BIND, [BOUND_LABEL["B-KINETIC"], BOUND_LABEL["B-HEAT"],
                    BOUND_LABEL["B-VAPOR"], BOUND_LABEL["B-CURV"]]),
        (GRP_REF, ["$\\epsilon$ = the smallest of them",
                   "$N_x = 10^4$  practical ceiling",
                   "what we run", LIT_LABEL]),
    ]


# =========================================================================

def build_context(T, args):
    ctx = dict(primary_T=args.primary_T, T_bounds=args.T_bounds, Lx=args.Lx,
               Rave=args.Rave, rneck=args.rneck)
    kw = dict(Lx=args.Lx, Ly=args.Ly, Lz=0.0, Rave=args.Rave, v_n=args.vn,
              safety=args.safety, eps_over_R=args.eps_over_R)

    for k in ("beta", "eps", "Nx", "nodes", "binding", "fkin", "Lstar"):
        ctx[k] = {}
    for t, _c, _l, _s in SWEEP_T:
        r = [ce.compute_eps(T0_C=t, alpha_c=a, **kw) for a in ALPHA]
        ctx["beta"][t] = np.array([x["beta_uns"] for x in r])
        ctx["eps"][t] = np.array([x["eps"] for x in r])
        ctx["Nx"][t] = np.array([float(x["Nx"]) for x in r])
        ctx["nodes"][t] = np.array([float(x["Nx"]) * float(x["Ny"]) for x in r])
        ctx["binding"][t] = [x["binding"] for x in r]
        ctx["fkin"][t] = np.array([x["kinetic_frac"] for x in r])
        ctx["Lstar"][t] = np.array([ce.beta_HK(t, a) * ce.Dv_T(t) for a in ALPHA])

    # beta_HK across the whole temperature range, to show how little it moves.
    Tedge = np.linspace(args.Tmin, args.Tmax, 40)
    bhk = np.array([[ce.beta_HK(t, a) for a in ALPHA] for t in Tedge])
    ctx["bhk_lo"], ctx["bhk_hi"] = bhk.min(axis=0), bhk.max(axis=0)

    # The four bounds themselves, at one temperature, as alpha_c varies.
    Tb = args.T_bounds
    b = np.array([ce.beta_HK(Tb, a) for a in ALPHA])
    rho_rat = ce.rho_vs_sat(Tb) / ce._RHO_ICE
    d0, Dv, s = ce.capillary_length(Tb), ce.Dv_T(Tb), args.safety
    ctx["bounds"] = {
        "B-HEAT":    s * ce.D_heat_of("mean") * b,
        "B-VAPOR":   s * Dv * b,
        "B-KINETIC": s * d0 / ((b / rho_rat) * args.vn),
        "B-CURV":    args.eps_over_R * args.Rave,
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
                   help="representative grain radius [m], for B-CURV and L*")
    p.add_argument("--rneck", type=float, default=1.4e-5,
                   help="neck radius [m], for the L* regime read-out")
    p.add_argument("--vn", type=float, default=3.415598e-9,
                   help="MEASURED front velocity [m/s]; K&P Eq. 45 is only a "
                        "constraint if the front actually moves")
    p.add_argument("--safety", type=float, default=0.5)
    p.add_argument("--eps_over_R", type=float, default=0.05)
    p.add_argument("--primary_T", type=float, default=-20.0,
                   help="the temperature whose binding bound is shaded")
    p.add_argument("--T_bounds", type=float, default=-20.0,
                   help="temperature the bound sweep in fig 2 is taken at")
    p.add_argument("--alpha_run", type=float, default=0.1)
    p.add_argument("--T_run", type=float, default=-20.0)
    p.add_argument("--eps_run", type=float, default=1.18e-7)
    p.add_argument("--Tmin", type=float, default=-40.0)
    p.add_argument("--Tmax", type=float, default=-1.0)
    p.add_argument("--width", type=float, default=SLIDE_W_IN)
    p.add_argument("--panel_width", type=float, default=3.0)
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--panel_legend", choices=("separate", "inline", "none"),
                   default="separate")
    p.add_argument("--label", default="reference case: 450 × 225 µm grain pair")
    p.add_argument("--out", type=Path, default=Path("studies/alpha_c_sizing"))
    args = p.parse_args()

    ctx = build_context(None, args)
    args.out.mkdir(parents=True, exist_ok=True)

    fs.emit(args.out, PANELS, None, ctx, width=args.width,
            panel_width=args.panel_width, dpi=args.dpi,
            panel_legend=args.panel_legend, master_groups=master_groups,
            suptitle="$\\alpha_c$ as a knob: what each value costs in "
                     "$\\beta_{sub}$, $\\epsilon$ and mesh\n"
                     f"{args.label}, measured $v_n$ = {args.vn:.3g} m/s")

    # --- the numbers ------------------------------------------------------
    csv = args.out / "alpha_c_sweep.csv"
    with open(csv, "w") as fh:
        fh.write("alpha_c,T_C,beta_sub_s_per_m,eps_m,binding,Nx,nodes_2D,"
                 "f_kin,Lstar_m\n")
        for i in range(0, len(ALPHA), 4):
            for t, _c, _l, _s in SWEEP_T:
                fh.write(f"{ALPHA[i]:.6e},{t:.1f},{ctx['beta'][t][i]:.6e},"
                         f"{ctx['eps'][t][i]:.6e},{ctx['binding'][t][i]},"
                         f"{ctx['Nx'][t][i]:.6e},{ctx['nodes'][t][i]:.6e},"
                         f"{ctx['fkin'][t][i]:.4f},{ctx['Lstar'][t][i]:.6e}\n")
    print(f"csv  -> {csv}")

    print(f"\n  Production run:  alpha_c = {args.alpha_run:g}, "
          f"beta_sub = {ctx['beta_run']:.3e} s/m, eps = {ctx['eps_run']:.3e} m, "
          f"Nx = {ctx['Nx_run']}, {ctx['nodes_run']/1e6:.1f} M nodes")
    print(f"  Practical ceiling Nx = {NX_USABLE:.0g} is eps = "
          f"{ctx['eps_usable']*1e9:.0f} nm")
    print(f"  beta_HK moves {100*(ctx['bhk_hi'][0]/ctx['bhk_lo'][0] - 1):.1f}% "
          f"over {args.Tmin:g}..{args.Tmax:g} C; rho_vs falls "
          f"{ce.rho_vs_sat(args.Tmax)/ce.rho_vs_sat(args.Tmin):.0f}x\n")

    for t, _c, _l, _s in SWEEP_T:
        y, e, f = ctx["Nx"][t], ctx["eps"][t], ctx["fkin"][t]
        j = int(np.argmin(y))
        sw = [ALPHA[i] for i in range(1, len(ALPHA))
              if ctx["binding"][t][i] != ctx["binding"][t][i - 1]]
        ok = np.flatnonzero(y <= NX_USABLE)
        print(f"  T = {t:6.1f} C : cheapest alpha_c = {ALPHA[j]:.3g} -> "
              f"Nx = {y[j]:.0f}, eps = {e[j]*1e9:.0f} nm")
        print(f"                 Nx = {y[0]:.3g} at 1e-4, {y[-1]:.3g} at 1; "
              f"bound switches at alpha_c = {['%.3g' % x for x in sw]}")
        print(f"                 Nx <= 1e4 for alpha_c in "
              f"[{ALPHA[ok[0]]:.3g}, {ALPHA[ok[-1]]:.3g}], "
              f"f_kin there {f[ok].min():.2f}-{f[ok].max():.2f}")


if __name__ == "__main__":
    main()
