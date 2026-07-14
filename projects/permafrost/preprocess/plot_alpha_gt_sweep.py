#!/usr/bin/env python3
"""Gibbs-Thomson parameters across the physical alpha_c range — no Libbrecht.

Two things on one figure:

1. THE SWEEP (top row): alpha_c in [1e-3, 1e-1] -> beta_sub(alpha; T), the
   M&F SI Eq. 9 phase-field parameters (tau_sub, mob_sub, alph_sub at the
   production eps), and the mesh consequence: the K&P thin-interface
   validity ceiling eps_max ~ 0.1*D_ice*beta_HK/(a1*a2) shrinks ~ 1/alpha,
   so FAST kinetics (high alpha) demand fine meshes. This panel shows the
   feasible (alpha, mesh) band explicitly.

2. THE ARRHENIUS PROPOSAL (bottom row): alpha_c(T) = A*exp(-Q/(R*T)),
   anchored so alpha spans [1e-3, 1e-1] smoothly across the temperature
   range of interest (defaults: 1e-1 at T_warm = -2 C, 1e-3 at
   T_cold = -40 C -> Q ~ 64 kJ/mol, comparable to ice self-diffusion
   activation energies — physically plausible shape). Bottom panels show
   the induced beta_sub(T) and the mesh requirement along the Arrhenius
   path: bounded everywhere, unlike the Libbrecht exp(-sigma0/sigma) model
   whose alpha spans ~30 decades and demands impossible meshes.

Anchors are CLI-tunable: --Twarm/--awarm/--Tcold/--acold.

Usage:
    python plot_alpha_gt_sweep.py [--eps 4.64e-7] [--Lx 3.848e-4]
                                  [--Twarm -2 --Tcold -40]
"""

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import comp_eps as ce  # noqa: E402

import numpy as np                 # noqa: E402
import matplotlib                  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt    # noqa: E402

A1A2 = ce._A1 * ce._A2


def beta_pair(T_C, alpha):
    """(beta_HK scaled, beta_sub unscaled) at T, alpha."""
    bhk = ce.beta_HK(T_C, alpha)
    return bhk, bhk / (ce.rho_vs_sat(T_C) / ce._RHO_ICE)


def pf_params(T_C, alpha, eps):
    """M&F SI Eq. 9 derived parameters at fixed eps."""
    bhk, _ = beta_pair(T_C, alpha)
    rho_rat = ce.rho_vs_sat(T_C) / ce._RHO_ICE
    d0_sub = ce.capillary_length(T_C) * rho_rat   # d0_sub0/(rho_i/rho_vs), K&P scaling
    lam = ce._A1 * eps / d0_sub
    Di = ce._K_I / ce._C_I if hasattr(ce, "_K_I") else 1.2722e-6
    Da = ce._K_A / ce._C_A if hasattr(ce, "_K_A") else 1.4855e-5
    Dstar = 0.5 * (Di + Da)
    Dv = ce.Dv_T(T_C)
    tau = eps * lam * (bhk / ce._A1 + ce._A2 * eps / Dstar + ce._A2 * eps / Dv)
    return lam, tau, eps / (3.0 * tau), lam / tau  # lam, tau, mob, alph


def eps_ceiling(T_C, alpha):
    """K&P thin-interface validity ceiling (heat-in-ice channel, corr>0.9)."""
    bhk, _ = beta_pair(T_C, alpha)
    Di = 1.2722e-6
    return 0.1 * Di * bhk / A1A2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eps", type=float, default=4.64e-7, help="production eps [m]")
    ap.add_argument("--Lx", type=float, default=3.848e-4,
                    help="domain length for the mesh panels [m]")
    ap.add_argument("--Twarm", type=float, default=-2.0)
    ap.add_argument("--awarm", type=float, default=1e-1)
    ap.add_argument("--Tcold", type=float, default=-40.0)
    ap.add_argument("--acold", type=float, default=1e-3)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).parent / "alpha_gt_sweep.png")
    args = ap.parse_args()

    # ----- Arrhenius fit through the two anchors: alpha = A*exp(-Q/(R*T))
    Tw, Tc = args.Twarm + 273.15, args.Tcold + 273.15
    QoverR = math.log(args.awarm / args.acold) / (1.0 / Tc - 1.0 / Tw)
    Apre = args.awarm * math.exp(QoverR / Tw)
    Q_kJ = QoverR * 8.314 / 1e3
    print(f"Arrhenius: alpha_c(T) = {Apre:.3e} * exp(-{Q_kJ:.1f} kJ/mol / RT)")
    print(f"  anchors: alpha({args.Twarm:g} C) = {args.awarm:g}, "
          f"alpha({args.Tcold:g} C) = {args.acold:g}")

    alpha = np.logspace(-3, -1, 200)
    T = np.linspace(args.Tcold, args.Twarm, 300)
    alpha_T = Apre * np.exp(-QoverR / (T + 273.15))

    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.6))
    fig.suptitle("Gibbs–Thomson parameters across the physical $\\alpha_c$ range "
                 "and a smooth Arrhenius $\\alpha_c(T)$ (no Libbrecht data)",
                 fontsize=12)
    C = {"m20": "#3d74d9", "m5": "#e8883a", "arr": "#3a8f8a", "ref": "#3f434a"}

    def mark_runs(ax):
        for a in [2e-3, 1e-2, 3e-2, 1e-1]:
            ax.axvline(a, color=C["ref"], lw=0.6, ls=":", alpha=0.5)

    # --- A: beta_sub(alpha) at two temperatures
    ax = axs[0, 0]
    for Tv, key, lab in [(-20.0, "m20", "−20 °C"), (-5.0, "m5", "−5 °C")]:
        b = [beta_pair(Tv, a)[1] for a in alpha]
        ax.loglog(alpha, b, lw=2, color=C[key], label=lab)
    mark_runs(ax)
    ax.set_xlabel(r"$\alpha_c$"); ax.set_ylabel(r"$\beta_{sub}$ [s/m]")
    ax.set_title(r"Kinetic coefficient $\beta_{sub}\propto 1/\alpha_c$"
                 "\n(dotted: this campaign's runs)")
    ax.legend(fontsize=9); ax.grid(alpha=0.25, which="both", lw=0.4)

    # --- B: derived PF parameters at fixed eps, T=-20
    ax = axs[0, 1]
    lam, tau, mob, alph = zip(*[pf_params(-20.0, a, args.eps) for a in alpha])
    ax.loglog(alpha, tau, lw=2, color=C["m20"], label=r"$\tau_{sub}$ [s]")
    ax.loglog(alpha, np.array(mob) * 1e9, lw=2, color=C["m5"],
              label=r"$mob_{sub}\times10^9$ [m/s]")
    ax.loglog(alpha, np.array(alph) / 1e6, lw=2, color=C["arr"],
              label=r"$\alpha_{sub}/10^6$ [1/s]")
    mark_runs(ax)
    ax.set_xlabel(r"$\alpha_c$")
    ax.set_title(f"M&F Eq. 9 parameters at eps = {args.eps:.2e}, −20 °C"
                 "\n(alph_sub saturates: eps-correction terms take over)")
    ax.legend(fontsize=9); ax.grid(alpha=0.25, which="both", lw=0.4)

    # --- C: mesh feasibility vs alpha
    ax = axs[0, 2]
    for Tv, key, lab in [(-20.0, "m20", "−20 °C"), (-5.0, "m5", "−5 °C")]:
        em = np.array([eps_ceiling(Tv, a) for a in alpha])
        Nx = args.Lx * math.sqrt(2) / em
        ax.loglog(alpha, Nx, lw=2, color=C[key], label=f"Nx required, {lab}")
    ax.axhline(args.Lx * math.sqrt(2) / args.eps, color=C["ref"], lw=1, ls="--",
               label=f"current mesh (eps = {args.eps:.2g})")
    mark_runs(ax)
    ax.set_xlabel(r"$\alpha_c$"); ax.set_ylabel(r"$N_x$ for validity (corr > 0.9)")
    ax.set_title("Mesh demanded by the K&P validity ceiling\n"
                 r"(eps$_{max}\propto 1/\alpha_c$: fast kinetics need fine meshes)")
    ax.legend(fontsize=8); ax.grid(alpha=0.25, which="both", lw=0.4)

    # --- D: the Arrhenius alpha_c(T)
    ax = axs[1, 0]
    ax.semilogy(T, alpha_T, lw=2.5, color=C["arr"])
    ax.axhspan(1e-3, 1e-1, color=C["arr"], alpha=0.08)
    ax.plot([args.Twarm, args.Tcold], [args.awarm, args.acold], "o",
            color=C["ref"], ms=6)
    ax.set_ylim(3e-4, 3e-1)
    ax.set_xlabel("T [°C]"); ax.set_ylabel(r"$\alpha_c$")
    ax.set_title(f"Proposed Arrhenius $\\alpha_c(T)$\n"
                 f"Q = {Q_kJ:.0f} kJ/mol (cf. ice self-diffusion ~60–70)")
    ax.grid(alpha=0.25, lw=0.4)

    # --- E: induced beta_sub(T) along the Arrhenius path
    ax = axs[1, 1]
    b_arr = [beta_pair(t, a)[1] for t, a in zip(T, alpha_T)]
    b_lo = [beta_pair(t, 1e-3)[1] for t in T]
    b_hi = [beta_pair(t, 1e-1)[1] for t in T]
    ax.fill_between(T, b_lo, b_hi, color=C["m20"], alpha=0.10,
                    label=r"envelope $\alpha_c\in[10^{-3},10^{-1}]$")
    ax.semilogy(T, b_arr, lw=2.5, color=C["arr"], label="Arrhenius path")
    ax.set_xlabel("T [°C]"); ax.set_ylabel(r"$\beta_{sub}$ [s/m]")
    ax.set_title(r"Induced $\beta_{sub}(T)$: smooth, bounded"
                 "\n(vs Libbrecht spanning ~30 decades)")
    ax.legend(fontsize=9); ax.grid(alpha=0.25, which="both", lw=0.4)

    # --- F: mesh requirement along the Arrhenius path
    ax = axs[1, 2]
    Nx_arr = [args.Lx * math.sqrt(2) / eps_ceiling(t, a)
              for t, a in zip(T, alpha_T)]
    ax.semilogy(T, Nx_arr, lw=2.5, color=C["arr"], label="Nx required (validity)")
    ax.axhline(args.Lx * math.sqrt(2) / args.eps, color=C["ref"], lw=1, ls="--",
               label=f"current mesh")
    ax.set_xlabel("T [°C]"); ax.set_ylabel(r"$N_x$")
    ax.set_title("Mesh along the Arrhenius path\n(bounded across the whole range)")
    ax.legend(fontsize=9); ax.grid(alpha=0.25, which="both", lw=0.4)

    for ax in axs.flat:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(args.out, dpi=150)
    print(f"figure -> {args.out}")

    # reference values along the path
    print(f"\n{'T [C]':>7} | {'alpha_c':>9} | {'beta_sub [s/m]':>14} | "
          f"{'eps_max [m]':>11} | {'Nx(385um)':>9}")
    for t in [-2, -5, -10, -20, -30, -40]:
        a = Apre * math.exp(-QoverR / (t + 273.15))
        _, b = beta_pair(t, a)
        em = eps_ceiling(t, a)
        print(f"{t:7.0f} | {a:9.3e} | {b:14.3e} | {em:11.3e} | "
              f"{args.Lx*math.sqrt(2)/em:9.0f}")


if __name__ == "__main__":
    main()
