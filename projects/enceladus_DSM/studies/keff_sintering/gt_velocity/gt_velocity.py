#!/usr/bin/env python3
"""gt_velocity.py — interface velocity from the Gibbs–Thomson condition against temperature.

    venv_enceladus/bin/python studies/keff_sintering/gt_velocity/gt_velocity.py

A plot to look at, not a manuscript figure. The interface condition the model
obeys is

    sigma = d0 kappa + (beta_sub + dbeta) v_n ,
    v_n   = (sigma - d0 kappa) / (beta_sub + dbeta) ,

with sigma the vapor supersaturation relative to flat ice, kappa the curvature
(positive for a convex grain; 1/R in two dimensions, as in the aggregates),
beta_sub = (rho_ice/rho_vs) sqrt(2 pi m / kB T) / alpha_c and the
thin-interface term dbeta = a1 a2 eps (1/D_T + 1/D_v) rho_ice/rho_vs.

Cases, with numbers taken from the aggregate runs (mean grain radius 50 um,
eps = 1 um, alpha_c = 1.020e-3, |sigma| up to about 3e-5 in the pores):
  neck    a concave neck of radius 5 um in vapor in equilibrium with the mean
          grain (sigma = d0 / 50 um): grows
  small   a convex grain of radius 15 um in the same vapor: sublimates
  flat    a flat surface at sigma = 1e-5: grows
Material laws from preprocess/comp_eps.py, the ones the solver uses.

Writes gt_velocity.png and gt_velocity.txt beside this file.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "preprocess"))
import comp_eps as ce                                               # noqa: E402

RHO_I, EPS, ALPHA = 919.0, 1.0e-6, 1.020e-3
A1, A2, D_T = 5.0, 0.1581, 7.78e-6
R_MEAN, R_NECK, R_SMALL, SIG_FLAT = 50e-6, 5e-6, 15e-6, 1e-5
DAY, R_GAS = 86400.0, 8.314


def coeffs(T, alpha=ALPHA):
    rvs = ce.rho_vs_sat(T)
    d0 = ce.capillary_length(T)
    beta = ce.beta_HK(T, alpha) * RHO_I / rvs
    dbeta = A1 * A2 * EPS * (1 / D_T + 1 / ce.Dv_T(T)) * RHO_I / rvs
    return d0, beta, dbeta, rvs


def main():
    T = np.linspace(-40, -5, 141)
    c = np.array([coeffs(t) for t in T])
    d0, beta, dbeta, rvs = c.T
    sig_bg = d0 / R_MEAN                                            # vapor in equilibrium with the mean grain
    cases = {
        rf"neck, concave radius {R_NECK * 1e6:.0f} $\mu$m (grows)": (sig_bg + d0 / R_NECK) / (beta + dbeta),
        rf"small grain, radius {R_SMALL * 1e6:.0f} $\mu$m (sublimates)": (sig_bg - d0 / R_SMALL) / (beta + dbeta),
        rf"flat surface, $\sigma = 10^{{-5}}$ (grows)": SIG_FLAT / (beta + dbeta),
    }
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for (lab, v), col in zip(cases.items(), ("C0", "C3", "C2")):
        ax[0].semilogy(T, np.abs(v) * 1e6 * DAY, color=col, lw=1.8, label=lab)
    ax[0].set(xlabel=r"$T$ [$^\circ$C]", ylabel=r"$|v_n|$ [$\mu$m per day]",
              title=rf"$\alpha_c = {ALPHA * 1e3:.2f}\times10^{{-3}}$, $\varepsilon = 1\ \mu$m")
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].grid(alpha=0.25, which="both")

    i = int(np.argmin(np.abs(T + 20)))
    for y, lab, col in ((1 / (beta + dbeta), r"$1/(\beta_\mathrm{sub}+\Delta\beta)$: mobility", "k"),
                        (rvs, r"$\rho_{vs}^{I}$: saturation vapor density", "C1"),
                        (d0, r"$d_0$: capillary length", "C4"),
                        (np.sqrt(T + 273.15), r"$\sqrt{T}$", "0.6")):
        ax[1].semilogy(T, y / y[i], color=col, lw=1.8, label=lab)
    ax[1].set(xlabel=r"$T$ [$^\circ$C]", ylabel=r"value relative to $-20\,^\circ$C", title="where the temperature dependence comes from")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].grid(alpha=0.25, which="both")
    fig.tight_layout(); fig.savefig(HERE / "gt_velocity.png", dpi=200)

    v = list(cases.values())[0]
    Tk = T + 273.15
    Ea = -R_GAS * np.polyfit(1 / Tk, np.log(v), 1)[0] / 1e3
    L = [f"alpha_c = {ALPHA:.3e}, eps = {EPS * 1e6:g} um; background vapor in equilibrium with a {R_MEAN * 1e6:.0f} um grain.", "",
         "T [C]   d0 [nm]   sigma_bg     beta_sub [s/m]   dbeta/beta   v_n neck [um/d]   v_n small grain [um/d]   v_n flat, sigma 1e-5 [um/d]"]
    for t in (-40, -30, -20, -10, -5):
        j = int(np.argmin(np.abs(T - t)))
        L.append(f"{t:4d}    {d0[j] * 1e9:.3f}    {sig_bg[j]:.2e}     {beta[j]:.3e}        {dbeta[j] / beta[j]:.3f}        "
                 + "            ".join(f"{vv[j] * 1e6 * DAY:+.3f}" for vv in cases.values()))
    L += ["", f"neck velocity, -5 C over -40 C: {v[-1] / v[0]:.1f}x;  apparent activation energy {Ea:.0f} kJ/mol",
          f"of that ratio: mobility {(1 / (beta + dbeta))[-1] / (1 / (beta + dbeta))[0]:.1f}x (rho_vs {rvs[-1] / rvs[0]:.1f}x), capillary length {d0[-1] / d0[0]:.2f}x"]
    txt = "\n".join(L)
    print(txt)
    (HERE / "gt_velocity.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
