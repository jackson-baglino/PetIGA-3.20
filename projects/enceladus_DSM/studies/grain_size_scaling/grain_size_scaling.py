#!/usr/bin/env python3
"""How the mesh, interface width and cost change with mean grain radius.

    venv_enceladus/bin/python studies/grain_size_scaling/grain_size_scaling.py

Sweeps R_ave from 0.05 um to 500 um through the SAME sizing code the campaign
uses (preprocess/comp_eps.py compute_eps, via generate_study_opts.derived_vn
and tau_sub_of), at -20 C, alpha_c = 1e-3, safety 0.5, periodic L = 40 R_ave.

Two ways to change grain size:

  SCALED (the campaign rule): R_feat = R_ave/25, so eps = R_ave/50 and the
      element count is the same at every size. This is what we would do.
  FIXED eps = 1 um (today's mesh spacing kept): the element count grows as
      R_ave^2, and below R_ave = 20 um eps exceeds the geometric bound
      eps <= R_ave/20 (B-CURV) -- the grains are not resolved at all.

Cost model, per 30-day-equivalent run at -20 C (measured 2026-09-30,
studies/keff_sintering/scaling/README.md): 118 core-seconds per million DoF
per phase-field step at 200k DoF/core, 5.2 s x 121 ranks per k_eff sample
(scaled with DoF), 179 samples, $0.012 per core-hour.
Steps = t_equiv/dtmax + 70 (early CFL-limited steps), dtmax = 1.09 tau_sub.

t_equiv is the time to reach the SAME sintering state as our 30 d at 50 um.
Attachment-limited, every time in the problem scales as R^2; with vapour
transport in series (rate ~ 1/(beta_HK + R/D_v)) the time is
    t_equiv = 30 d * (R/50um)^2 * (1 + R/L*) / (1 + 50um/L*),  L* = D_v beta_HK.
This is an estimate of the model's own scaling, not of real snow -- see the
README for what the model leaves out at sub-micron sizes.

Writes grain_size_scaling.csv and grain_size_scaling.{png,pdf} here.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "preprocess"))
from comp_eps import (compute_eps, Dv_T, beta_HK, capillary_length,  # noqa: E402
                      rho_vs_sat)
from generate_study_opts import derived_vn, tau_sub_of, DTMAX_OVER_TAU  # noqa: E402

T_C, ALPHA, SAFETY, L_OVER_R = -20.0, 1.0e-3, 0.5, 40.0
R0 = 50e-6                     # campaign grain radius
T0_DAYS = 30.0                 # campaign run length at R0
DOF_PER_NODE = 3
CORE_S_PER_MDOF_STEP = 23.4 * 121 / 24.0   # measured, 121 ranks, 24.0M DoF
KEFF_CORE_S_PER_MDOF = 5.2 * 121 / 24.0
N_KEFF, EARLY_STEPS, RATE = 179, 70, 0.012
TARGET_DOFS_PER_CORE = 200_000
FIXED_EPS = 1.0e-6

# reference lengths
L_STAR = Dv_T(T_C) * beta_HK(T_C, ALPHA)          # attachment/diffusion crossover
D0 = capillary_length(T_C)                         # ice capillary length
LAMBDA_AIR = 0.0665e-6 * (T_C + 273.15) / 293.15   # air mean free path, 1 atm (~T)
A_MOL = 0.32e-9                                    # water molecular spacing
N_VAPOR = rho_vs_sat(T_C) / 2.99e-26               # vapour molecules per m^3


def t_equiv(R):
    return T0_DAYS * 86400 * (R / R0) ** 2 * (1 + R / L_STAR) / (1 + R0 / L_STAR)


def run_cost(dof, steps):
    """$ per run: phase-field steps + k_eff samples, both ~ DoF."""
    core_s = dof / 1e6 * (steps * CORE_S_PER_MDOF_STEP + N_KEFF * KEFF_CORE_S_PER_MDOF)
    return core_s / 3600 * RATE


def scaled(R):
    L = L_OVER_R * R
    p = compute_eps(Lx=L, Ly=L, Rave=R, T0_C=T_C, alpha_c=ALPHA, safety=SAFETY,
                    v_n=derived_vn(T_C, ALPHA, R / 25))
    tau = tau_sub_of(p["eps"], p["beta_uns"], p["d0"])
    dt = DTMAX_OVER_TAU * tau
    te = t_equiv(R)
    steps = te / dt + EARLY_STEPS
    dof = DOF_PER_NODE * p["Nx"] * p["Ny"]
    return dict(R_m=R, L_m=L, eps_m=p["eps"], binding=p["binding"], h_m=p["eps"] / math.sqrt(2),
                band_1_99_m=9.2 * p["eps"], Nx=p["Nx"], dof=dof,
                ranks=math.ceil(dof / TARGET_DOFS_PER_CORE), tau_sub_s=tau, dtmax_s=dt,
                t_equiv_s=te, steps=steps, usd=run_cost(dof, steps),
                delta_heat=p["delta_heat"], R_over_Lstar=R / L_STAR)


def fixed(R):
    """Keep eps = 1 um. Valid only while eps <= R/20 (B-CURV)."""
    L = L_OVER_R * R
    eps = FIXED_EPS
    p = compute_eps(Lx=L, Ly=L, Rave=R, T0_C=T_C, alpha_c=ALPHA, safety=SAFETY,
                    v_n=derived_vn(T_C, ALPHA, 2e-6))
    nx = math.ceil(math.sqrt(2) * L / eps)
    dof = DOF_PER_NODE * nx * nx
    tau = tau_sub_of(eps, p["beta_uns"], p["d0"])
    steps = t_equiv(R) / (DTMAX_OVER_TAU * tau) + EARLY_STEPS
    return dict(Nx=nx, dof=dof, steps=steps, usd=run_cost(dof, steps),
                valid=eps <= 0.05 * R)


def main():
    Rs = np.logspace(math.log10(0.05e-6), math.log10(500e-6), 121)
    A = [scaled(R) for R in Rs]
    B = [fixed(R) for R in Rs]
    rows = [{**a, "fixed_eps_Nx": b["Nx"], "fixed_eps_dof": b["dof"],
             "fixed_eps_steps": b["steps"], "fixed_eps_usd": b["usd"],
             "fixed_eps_valid": b["valid"]} for a, b in zip(A, B)]
    with open(HERE / "grain_size_scaling.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    def at(R):
        return min(A, key=lambda r: abs(math.log(r["R_m"] / R)))
    um = Rs * 1e6
    g = lambda k, X=A: np.array([x[k] for x in X])
    valid = np.array([b["valid"] for b in B])

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                         "legend.fontsize": 8})
    INK, MUTED = "#1a1a1a", "#6b6b6b"
    C1, C2, C3, C4 = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.6), constrained_layout=True)

    def mark(a):
        a.axvspan(um[0], 1.0, color="#f2f2f2", zorder=0, lw=0)
        a.axvspan(L_STAR * 1e6, um[-1], color="#eef4fa", zorder=0, lw=0)
        for x, s in ((50, "campaign\n50 µm"), (0.1, "0.1 µm")):
            a.axvline(x, color=MUTED, lw=0.8, ls=":")
        a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(um[0], um[-1])
        a.grid(True, which="major", alpha=0.25)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.set_xlabel(r"mean grain radius $R_\mathrm{ave}$  [µm]")

    # (a) length scales
    a = ax[0, 0]; mark(a)
    a.plot(um, um, color=INK, lw=2, label=r"grain radius $R_\mathrm{ave}$")
    a.plot(um, 0.49 * um, color=INK, lw=1.2, ls="--", label="median pore throat (≈0.49 R)")
    a.plot(um, g("band_1_99_m") * 1e6, color=C2, lw=2, label=r"visible interface, 1–99 % (9.2 $\varepsilon$)")
    a.plot(um, g("eps_m") * 1e6, color=C1, lw=2, label=r"$\varepsilon$ = R/50")
    a.plot(um, g("h_m") * 1e6, color=C3, lw=2, label=r"element size $h=\varepsilon/\sqrt{2}$")
    for y, s in ((LAMBDA_AIR * 1e6, "air mean free path"), (D0 * 1e6, r"capillary length $d_0$"),
                 (A_MOL * 1e6, "water molecule")):
        a.axhline(y, color=MUTED, lw=0.9, ls="-.")
        a.text(um[-1] * 0.8, y * 1.15, s, ha="right", va="bottom", color=MUTED, fontsize=8)
    a.set_ylabel("length  [µm]")
    a.set_title("(a) length scales (mesh scales with the grains)", loc="left")
    a.legend(frameon=False, loc="upper left")

    # (b) mesh size
    a = ax[0, 1]; mark(a)
    a.plot(um, g("dof"), color=C1, lw=2, label=r"scaled: $\varepsilon$ = R/50 (campaign rule)")
    fd = g("dof", B)
    a.plot(um[valid], fd[valid], color=C2, lw=2, label=r"fixed $\varepsilon$ = 1 µm")
    a.axvline(20, color=C2, lw=0.8, ls="--")
    a.text(19, 0.55, r"fixed $\varepsilon$ impossible below 20 µm" "\n" r"($\varepsilon$ > R/20: grains unresolved)",
           ha="right", color=C2, fontsize=8, transform=a.get_xaxis_transform())
    a.set_ylabel("unknowns per run (3 per node)")
    a.set_title(f"(b) mesh: Nx = {A[0]['Nx']} at every size when scaled", loc="left")
    a.legend(frameon=False, loc="upper left")

    # (c) time scales
    a = ax[1, 0]; mark(a)
    a.plot(um, g("t_equiv_s"), color=INK, lw=2, label="simulated time for our 30-d state")
    a.plot(um, g("dtmax_s"), color=C1, lw=2, label=r"time-step cap $1.09\,\tau_\mathrm{sub}$")
    st = g("steps")
    a.text(0.03, 0.80, f"steps per run: {st[0]:.0f} at 0.05 µm, {at(10e-6)['steps']:.0f} at 10 µm,\n"
           f"{at(50e-6)['steps']:.0f} at 50 µm, {st[-1]:.0f} at 500 µm",
           transform=a.transAxes, ha="left", va="top", fontsize=8, color=INK)
    for y, s in ((60, "1 min"), (86400, "1 day")):
        a.axhline(y, color=MUTED, lw=0.8, ls="-.")
        a.text(um[0] * 1.2, y * 1.2, s, color=MUTED, fontsize=8)
    a.set_ylabel("time  [s]")
    a.set_title("(c) all times scale as R², so the step count barely moves", loc="left")
    a.legend(frameon=False, loc="upper left")

    # (d) cost
    a = ax[1, 1]; mark(a)
    a.plot(um, g("usd"), color=C1, lw=2, label=r"scaled: $\varepsilon$ = R/50")
    fu = g("usd", B)
    a.plot(um[valid], fu[valid], color=C2, lw=2, label=r"fixed $\varepsilon$ = 1 µm")
    a.axvline(20, color=C2, lw=0.8, ls="--")
    a.set_ylabel("$ per run at −20 °C (30-d-equivalent)")
    a.set_title("(d) cost per run", loc="left")
    a.legend(frameon=False, loc="upper left")

    for a in (ax[0, 1], ax[1, 1]):
        a.text(0.35, 0.02, "continuum\nquestionable", transform=a.get_xaxis_transform(),
               ha="center", va="bottom", color=MUTED, fontsize=7.5)
        a.text(math.sqrt(L_STAR * 1e6 * um[-1]), 0.02, "vapour-\ndiffusion\nlimited",
               transform=a.get_xaxis_transform(), ha="center", va="bottom", color=MUTED, fontsize=7.5)
    fig.savefig(HERE / "grain_size_scaling.png", dpi=200)
    fig.savefig(HERE / "grain_size_scaling.pdf")

    print(f"L* = {L_STAR*1e6:.0f} um, d0 = {D0*1e9:.2f} nm, air mfp = {LAMBDA_AIR*1e9:.0f} nm")
    for R in (0.05e-6, 0.1e-6, 1e-6, 10e-6, 50e-6, 500e-6):
        r = at(R); b = B[A.index(r)]
        print(f"R={r['R_m']*1e6:8.3f} um  L={r['L_m']*1e6:9.2f} um  eps={r['eps_m']*1e9:9.2f} nm "
              f"({r['binding']})  Nx={r['Nx']}  DoF={r['dof']:.3g}  tau={r['tau_sub_s']:.3g}s  "
              f"t_eq={r['t_equiv_s']:.3g}s  steps={r['steps']:.0f}  ${r['usd']:.2f}  "
              f"delta={r['delta_heat']:.3f} | fixed-eps DoF={b['dof']:.3g} ${b['usd']:.3g} valid={b['valid']}")
    print(f"vapour molecules in a 100 nm cube at -20 C: {N_VAPOR*1e-21:.0f}")
    print(f"wrote {HERE/'grain_size_scaling.png'}, .pdf, .csv")


if __name__ == "__main__":
    main()
