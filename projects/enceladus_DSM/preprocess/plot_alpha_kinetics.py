#!/usr/bin/env python3
"""How alpha_c(T) propagates into eps, the mesh, the timestep, and the REGIME.

WHY THIS EXISTS
---------------
alpha_c is the model's one genuinely free parameter, and it is not a detail
that sits off to the side: it sets beta_sub, beta_sub sets the binding K&P
bound on eps, and eps sets the mesh. It also decides which transport regime
the simulation is in, through the crossover length

    L* = beta_HK * D_v          (numerically identical to comp_eps Eq. 43c)

Below L* attachment kinetics limit the neck flux; above it vapour diffusion
does. That distinction is not cosmetic here, because the Kuczynski
evaporation-condensation exponent r ~ t^(1/3) is DERIVED for the
attachment-limited case. Run diffusion-limited and there is no reason to
expect 1/3 at all.

So a single number chosen for convenience silently selects the physics, the
cost, and whether the headline result is even well posed. These panels make
that chain visible in one figure.

THE alpha_c MODEL BEING PLOTTED
-------------------------------
A smooth two-anchor Arrhenius (comp_eps.alpha_arrhenius) pinned to the
alpha_c band that Libbrecht (2017) and Braun et al. (2024) support, roughly
[1e-4, 1e-3]. It is NOT a fit to the Libbrecht sigma0 table, for reasons
panel B makes visible: that table is 10 sparse points, non-monotonic around
-6/-7 C, and its functional form alpha = exp(-sigma0/sigma) spans a factor
27 in sigma0 across -2..-40 C, so ln(alpha) spans a factor 27 too. No single
reference sigma can hold that inside one decade. Libbrecht therefore fixes
the SIGN of the temperature dependence; the literature bounds fix the
magnitude.

Usage:
    python preprocess/plot_alpha_kinetics.py --out studies/.../verification
    python preprocess/plot_alpha_kinetics.py --Lx 3.87e-4 --Ly 1.21e-4 \
        --Rave 8.675e-5 --label "Molaro pair"
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

# Okabe-Ito, the repo standard (plot_neck_convergence.py, fit_neck_growth.py).
# A published CVD-safe set; assigned in fixed order, never cycled.
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

# Target neck-resolution floor r/R = sqrt(12*eps/R); 0.09 is the value the
# Demmenie geometry is built around (see studies/sinter_exponent/PLAN.md).
FLOOR_TARGET = 0.09


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
    p.add_argument("--rneck", type=float, default=4.5e-5,
                   help="neck radius used for the L/L* regime read-out [m]")
    p.add_argument("--safety", type=float, default=0.5)
    p.add_argument("--vn_feature", type=float, default=None,
                   help="R_feat for the Eq.(45) bound [m]; default Rave/50")
    p.add_argument("--label", default="Demmenie 1 mm pair")
    p.add_argument("--Tmin", type=float, default=-40.0)
    p.add_argument("--Tmax", type=float, default=-1.0)
    p.add_argument("--out", type=Path, default=Path("."))
    args = p.parse_args()

    rfeat = args.vn_feature if args.vn_feature else args.Rave / 50.0
    T = np.linspace(args.Tmin, args.Tmax, 160)

    # --- alpha_c models -----------------------------------------------------
    a_arr = np.array([ce.alpha_arrhenius(t) for t in T])
    a_unclamped = np.array([ce.alpha_arrhenius(t, clamp=None) for t in T])
    sig0 = np.array([ce.sigma0(t) for t in T])
    # Libbrecht family at several characteristic supersaturations.
    sig_family = [1.0e-2, 3.0e-3, 1.0e-3, 4.5e-4]
    a_lib = {s: np.array([ce.alpha_libbrecht(t, s) for t in T]) for s in sig_family}

    # --- propagate each alpha_c through the sizer ---------------------------
    eps, nodes, tau, beta_sub, Lstar, floor, kin = ([] for _ in range(7))
    bounds = {"Eq.43a heat": [], "Eq.43c vapour": [], "Eq.45 kinetic": []}
    for t, a in zip(T, a_arr):
        bhk = ce.beta_HK(t, a)
        chi = ce.rho_vs_sat(t) / ce._RHO_ICE
        bsub = bhk / chi
        vn = ce.capillary_length(t) / (bsub * rfeat)
        r = ce.compute_eps(Lx=args.Lx, Ly=args.Ly, Lz=0.0, Rave=args.Rave,
                           T0_C=t, alpha_c=a, v_n=vn, safety=args.safety)
        eps.append(r["eps"]); tau.append(r["tau_sub"]); beta_sub.append(bsub)
        nodes.append(r["Nx"] * max(r["Ny"], 1))
        kin.append(r["kinetic_frac"])
        Lstar.append(bhk * ce.Dv_T(t))
        floor.append(math.sqrt(12.0 * r["eps"] / args.Rave))
        bounds["Eq.43a heat"].append(ce.D_heat_of("mean") * bhk)
        bounds["Eq.43c vapour"].append(ce.Dv_T(t) * bhk)
        bounds["Eq.45 kinetic"].append(rfeat)
    eps, nodes, tau = map(np.array, (eps, nodes, tau))
    beta_sub, Lstar, floor, kin = map(np.array, (beta_sub, Lstar, floor, kin))

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.4))
    fig.suptitle(f"How $\\alpha_c(T)$ propagates into the discretisation and the "
                 f"transport regime — {args.label}", fontsize=11.5, color=INK)

    # A: alpha_c(T), Arrhenius over the Libbrecht form -----------------------
    a = axes[0][0]
    a.axhspan(ce.ALPHA_LIT_LO, ce.ALPHA_LIT_HI, color=C[0], alpha=0.10,
              label="literature band (Libbrecht 2017; Braun 2024)")
    for i, s in enumerate(sig_family):
        a.plot(T, np.maximum(a_lib[s], 1e-40), lw=1.2, ls=(0, (4, 2)),
               color=MUTED, alpha=0.35 + 0.15 * i)
        j = int(0.62 * len(T))
        if a_lib[s][j] > 1e-38:
            a.annotate(f"$\\sigma$={s:g}", (T[j], a_lib[s][j]), fontsize=7,
                       color=MUTED, ha="left", va="bottom")
    a.plot(T, a_unclamped, lw=1.0, ls=":", color=C[1], alpha=0.8,
           label="Arrhenius, unclamped")
    a.plot(T, a_arr, lw=2.4, color=C[1], label="Arrhenius $\\alpha_c(T)$ (used)")
    A, Ea = ce.arrhenius_params()
    a.annotate(f"$E_a$ = {Ea/1000:.1f} kJ/mol", (0.03, 0.06), xycoords="axes fraction",
               fontsize=8, color=C[1])
    _ax(a, "T [°C]", "$\\alpha_c$  [-]",
        "A.  $\\alpha_c$: Arrhenius vs the Libbrecht form $e^{-\\sigma_0/\\sigma}$")
    a.set_ylim(1e-12, 1e0)
    a.legend(fontsize=7, frameon=False, loc="upper left")

    # B: the actual Libbrecht data -------------------------------------------
    a = axes[0][1]
    a.plot(T, sig0, lw=2.0, color=C[2], label="$\\sigma_0(T)$ interpolant (as coded)")
    a.plot(ce._SIG0_T[:9], ce._SIG0_S[:9], "o", ms=6, color=C[2],
           markerfacecolor="white", markeredgewidth=1.6, label="Libbrecht table points")
    a.annotate("non-monotonic kink\nat −6/−7 °C (digitisation\nartifact, not physics)",
               (-6.5, 8.0e-3), xytext=(-33, 2.4e-2), fontsize=7.5,
               color=C[1], ha="left", va="center",
               arrowprops=dict(arrowstyle="->", color=C[1], lw=1.0))
    _ax(a, "T [°C]", "$\\sigma_0$  [-]",
        "B.  Libbrecht's data: 10 points, and a kink")
    a.legend(fontsize=7.5, frameon=False, loc="upper right")

    # C: eps and which bound binds -------------------------------------------
    a = axes[0][2]
    for i, (k, v) in enumerate(bounds.items()):
        a.plot(T, v, lw=1.3, ls=(0, (5, 2)), color=MUTED, alpha=0.55)
        j = int((0.06 + 0.26 * i) * len(T))          # stagger in x, not stacked
        a.annotate(k, (T[j], v[j]), fontsize=7, color=MUTED, ha="left", va="bottom")
    a.axhline(args.Rave, lw=1.3, ls=(0, (5, 2)), color=MUTED, alpha=0.55)
    a.annotate("geometric $R_{ave}$", (T[int(0.72 * len(T))], args.Rave), fontsize=7,
               color=MUTED, ha="left", va="bottom")
    a.plot(T, eps, lw=2.4, color=C[0],
           label=f"$\\epsilon$ = {args.safety:g} × min(bounds)  → Eq.45 binds")
    _ax(a, "T [°C]", "$\\epsilon$  [m]",
        "C.  Interface width: with literature $\\alpha_c$, K&P stops binding")
    a.legend(fontsize=7.5, frameon=False, loc="lower left")

    # D: mesh cost ------------------------------------------------------------
    # With alpha_c in the literature band the K&P ceilings go slack, so the
    # mesh is no longer set by them -- it is set by the requirement that the
    # NECK FILLET be resolved, r/R >= sqrt(12*eps/R). Plotting only the K&P
    # mesh would show a flat line and hide the real constraint.
    a = axes[1][0]
    eps_floor = args.Rave * FLOOR_TARGET ** 2 / 12.0        # eps for the target floor
    n_floor = (math.ceil(math.sqrt(2) * args.Lx / eps_floor)
               * math.ceil(math.sqrt(2) * args.Ly / eps_floor))
    a.plot(T, nodes / 1e6, lw=2.4, color=C[3],
           label="from the K&P ceilings ($\\epsilon$ in panel C)")
    a.axhline(n_floor / 1e6, lw=2.0, ls=(0, (4, 2)), color=C[1],
              label=f"to resolve a neck at $r/R$ = {FLOOR_TARGET:g}")
    a.annotate(f"{n_floor/1e6:.1f}M nodes — THIS is what sizes the mesh now",
               (args.Tmin + 1, n_floor / 1e6 * 1.35), fontsize=7.5, color=C[1])
    a.annotate(f"{nodes[0]/1e6:.2f}M", (args.Tmin + 1, nodes[0] / 1e6 * 1.3),
               fontsize=7.5, color=C[3])
    a.set_ylim(nodes.min() / 1e6 / 4, n_floor / 1e6 * 8)
    _ax(a, "T [°C]", "mesh nodes  [millions]",
        "D.  What actually sets the mesh (cost $\\propto \\epsilon^{-2}$)")
    a.legend(fontsize=7.5, frameon=False, loc="lower right")

    # E: timestep and kinetic resistance -------------------------------------
    a = axes[1][1]
    a.plot(T, tau, lw=2.4, color=C[4], label="$\\tau_{sub}$  (sets $dt_{max}$) [s]")
    a.plot(T, beta_sub, lw=2.0, ls=(0, (4, 2)), color=C[5],
           label="$\\beta_{sub}$  (kinetic resistance) [s/m]")
    _ax(a, "T [°C]", "value", "E.  Timestep scale and attachment resistance")
    a.legend(fontsize=7.5, frameon=False, loc="upper right")

    # F: THE REGIME MAP -------------------------------------------------------
    a = axes[1][2]
    a.fill_between(T, 1e-9, Lstar, color=C[2], alpha=0.12)
    a.plot(T, Lstar * 1e6, lw=2.4, color=C[2],
           label="$L^* = \\beta_{HK} D_v$  (crossover)")
    a.axhline(args.rneck * 1e6, lw=1.8, ls=(0, (4, 2)), color=C[1])
    a.annotate(f"neck radius {args.rneck*1e6:.0f} µm", (args.Tmin + 1, args.rneck * 1e6),
               fontsize=7.5, color=C[1], ha="left", va="bottom")
    a.axhline(args.Rave * 1e6, lw=1.8, ls=(0, (1, 3)), color=C[0])
    a.annotate(f"grain radius {args.Rave*1e6:.0f} µm", (args.Tmin + 1, args.Rave * 1e6),
               fontsize=7.5, color=C[0], ha="left", va="bottom")
    a.set_ylim(min(Lstar.min() * 1e6, args.rneck * 1e6) / 4,
               max(Lstar.max() * 1e6, args.Rave * 1e6) * 4)
    _ax(a, "T [°C]", "length  [µm]",
        "F.  Regime: $L^*$ above the neck ⇒ attachment-limited")
    a.legend(fontsize=7.5, frameon=False, loc="lower left")

    fig.tight_layout(rect=(0, 0, 1, 0.965))
    args.out.mkdir(parents=True, exist_ok=True)
    png = args.out / "alpha_kinetics.png"
    fig.savefig(png, dpi=150)
    print(f"plot -> {png}")

    # A table view, so identity/values are never colour-alone.
    csv = args.out / "alpha_kinetics.csv"
    with open(csv, "w") as fh:
        fh.write("T_C,alpha_c,eps_m,Nx_times_Ny,tau_sub_s,beta_sub_s_per_m,"
                 "Lstar_m,rneck_over_Lstar,res_floor_rR,kinetic_frac\n")
        for i in range(0, len(T), 4):
            fh.write(f"{T[i]:.2f},{a_arr[i]:.4e},{eps[i]:.4e},{nodes[i]:.0f},"
                     f"{tau[i]:.4e},{beta_sub[i]:.4e},{Lstar[i]:.4e},"
                     f"{args.rneck/Lstar[i]:.4f},{floor[i]:.4f},{kin[i]:.4f}\n")
    print(f"csv  -> {csv}")

    print(f"\n  Arrhenius: Ea = {Ea/1000:.2f} kJ/mol, A = {A:.4e}, "
          f"clamp [{ce.ALPHA_LIT_LO:g}, {ce.ALPHA_LIT_HI:g}]")
    for t in (-3.0, -20.0):
        i = int(np.argmin(abs(T - t)))
        print(f"  T = {t:6.1f} C : alpha_c = {a_arr[i]:.3e}, eps = {eps[i]:.3e} m, "
              f"{nodes[i]/1e6:.2f}M nodes, tau_sub = {tau[i]:.3e} s, "
              f"r_neck/L* = {args.rneck/Lstar[i]:.2f} "
              f"({'attachment' if args.rneck/Lstar[i] < 1 else 'diffusion'}-limited)")


if __name__ == "__main__":
    main()
