#!/usr/bin/env python3
"""Plot every temperature-dependent Gibbs-Thomson / kinetic parameter.

Companion to comp_eps.py (imports its functions, which are exact ports of
material_properties.c / monitoring.c where noted there). Produces a 2x3
panel figure over T in [-40, -1] C and prints a reference table at the
temperatures we run (-5, -20 C).

The central subtlety this figure exposes: with the code's Libbrecht
attachment model (flag_Tdep branch), alpha = exp(-sigma0(T)/sigma) depends
on the LOCAL supersaturation sigma, not just T. At Gibbs-Thomson-scale
driving (sigma ~ d0*kappa ~ 1e-5..1e-3, the regime of a saturated sintering
experiment), alpha collapses by tens to thousands of orders of magnitude
relative to the constant alpha_c ~ 2e-3 the mesh derivations have assumed.
beta, the eps bound (K&P Eq. 45), and every timescale inherit that swing —
so "the" mesh size at a given T is undefined until a characteristic sigma
is chosen. Sweep before trusting.

Usage:
    python plot_gt_temperature_dependence.py [--out FIG.png] [--vn 1e-9]
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


def beta_unscaled(T_C: float, alpha_c: float) -> float:
    """K&P beta0 (= M&F beta_sub) [s/m] = beta_HK / (rho_vs/rho_ice)."""
    rho_rat = ce.rho_vs_sat(T_C) / ce._RHO_ICE
    return ce.beta_HK(T_C, alpha_c) / rho_rat


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).parent / "gt_temperature_dependence.png")
    ap.add_argument("--vn", type=float, default=1.0e-9,
                    help="Front velocity for the Eq. 45 eps bound [m/s]")
    ap.add_argument("--Rref", type=float, default=36.5e-6,
                    help="Reference grain radius for tau_kin [m]")
    args = ap.parse_args()

    T = np.linspace(-40.0, -1.0, 400)
    sigmas = [1e-4, 1e-3, 1e-2, 1e-1]          # characteristic supersaturations
    alpha_const = 2e-3                          # the constant-alpha assumption
    Tmarks = [-5.0, -20.0]

    rho_vs = np.array([ce.rho_vs_sat(t) for t in T])
    Dv     = np.array([ce.Dv_T(t) for t in T])
    d0_scr = np.array([ce.capillary_length(t) for t in T])
    d0_cod = np.array([ce.capillary_length_code(t) for t in T])
    sig0   = np.array([ce.sigma0(t) for t in T])

    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle("Temperature dependence of Gibbs-Thomson / kinetic parameters "
                 "(ports of material_properties.c + monitoring.c flag_Tdep)",
                 fontsize=12)

    def mark(ax):
        for tm in Tmarks:
            ax.axvline(tm, color="gray", ls=":", lw=0.8)

    ax = axs[0, 0]
    ax.semilogy(T, rho_vs, label=r"$\rho_{vs}(T)$ [kg/m$^3$]")
    ax.semilogy(T, ce._RHO_ICE / rho_vs, label=r"$\rho_i/\rho_{vs}$ [-]")
    ax.set_title("Saturation vapor density (ASHRAE)")
    ax.set_xlabel("T [°C]"); ax.legend(); mark(ax); ax.grid(alpha=0.3)

    ax = axs[0, 1]
    ax.plot(T, Dv * 1e5, label=r"$D_v(T)\times10^5$ [m$^2$/s]")
    ax.plot(T, d0_scr * 1e9, label=r"$d_0$ script: $\gamma V_m/RT$ [nm]")
    ax.plot(T, d0_cod * 1e9, "--", label=r"$d_0$ code: $2.548\!\cdot\!10^{-7}/T_K$ [nm]")
    ax.set_title("Vapor diffusivity & capillary length\n(note ~1% script-vs-code d0 offset)")
    ax.set_xlabel("T [°C]"); ax.legend(); mark(ax); ax.grid(alpha=0.3)

    ax = axs[0, 2]
    ax.semilogy(T, sig0)
    ax.set_title(r"Libbrecht critical supersaturation $\sigma_0(T)$"
                 "\n(Sigma0 lookup; non-monotonic bump at -6/-7 °C is in the table)")
    ax.set_xlabel("T [°C]"); mark(ax); ax.grid(alpha=0.3)

    ax = axs[1, 0]
    for s in sigmas:
        a = np.array([ce.alpha_libbrecht(t, s) for t in T])
        ax.semilogy(T, a, label=rf"$\sigma = {s:g}$")
    ax.axhline(alpha_const, color="k", ls="--", lw=1,
               label=rf"constant $\alpha_c = {alpha_const:g}$")
    ax.set_ylim(1e-32, 3)
    ax.set_title(r"Attachment coefficient $\alpha = e^{-\sigma_0(T)/\sigma}$"
                 "\n(floor 1e-30 = code's cutoff)")
    ax.set_xlabel("T [°C]"); ax.legend(fontsize=8); mark(ax); ax.grid(alpha=0.3)

    ax = axs[1, 1]
    for s in sigmas:
        b = np.array([beta_unscaled(t, ce.alpha_libbrecht(t, s)) for t in T])
        ax.semilogy(T, b, label=rf"$\sigma = {s:g}$")
    b_const = np.array([beta_unscaled(t, alpha_const) for t in T])
    ax.semilogy(T, b_const, "k--", lw=1, label=rf"constant $\alpha_c = {alpha_const:g}$")
    ax.set_title(r"Kinetic coefficient $\beta_{sub}(T;\sigma)$ [s/m] (K&P $\beta_0$)")
    ax.set_xlabel("T [°C]"); ax.legend(fontsize=8); mark(ax); ax.grid(alpha=0.3)

    ax = axs[1, 2]
    for s in sigmas:
        b = np.array([beta_unscaled(t, ce.alpha_libbrecht(t, s)) for t in T])
        eps_bound = d0_scr / (b * args.vn)
        ax.semilogy(T, eps_bound, label=rf"$\sigma = {s:g}$")
    eps_bound_const = d0_scr / (b_const * args.vn)
    ax.semilogy(T, eps_bound_const, "k--", lw=1,
                label=rf"constant $\alpha_c = {alpha_const:g}$")
    ax.set_title(rf"K&P Eq. 45 eps bound $d_0/(\beta_{{sub}} v_n)$ [m], "
                 rf"$v_n = {args.vn:g}$"
                 "\n(mesh: h = 0.5·bound/√2; smaller bound = finer mesh)")
    ax.set_xlabel("T [°C]"); ax.legend(fontsize=8); mark(ax); ax.grid(alpha=0.3)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out, dpi=150)
    print(f"figure -> {args.out}")

    # Reference table at run temperatures
    print(f"\n{'':>10} | {'sigma0':>10} | {'rho_vs':>10} | {'d0 [m]':>10} | "
          f"{'beta(a=2e-3)':>12} | {'beta(s=1e-3)':>12} | {'beta(s=1e-2)':>12}")
    for t in Tmarks:
        b_c  = beta_unscaled(t, alpha_const)
        b_s3 = beta_unscaled(t, ce.alpha_libbrecht(t, 1e-3))
        b_s2 = beta_unscaled(t, ce.alpha_libbrecht(t, 1e-2))
        print(f"T = {t:5.1f} C | {ce.sigma0(t):10.3e} | {ce.rho_vs_sat(t):10.3e} | "
              f"{ce.capillary_length(t):10.3e} | {b_c:12.3e} | {b_s3:12.3e} | {b_s2:12.3e}")
    print("\n(beta(s=...) uses the code's Libbrecht alpha = exp(-sigma0/sigma); "
          "note the orders-of-magnitude spread vs the constant-alpha column — "
          "the characteristic sigma must be chosen before the mesh can be sized.)")


if __name__ == "__main__":
    main()
