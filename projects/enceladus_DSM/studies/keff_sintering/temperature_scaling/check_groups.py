#!/usr/bin/env python3
"""check_groups.py — the dimensionless groups of the aggregate runs at each temperature.

    venv_enceladus/bin/python studies/keff_sintering/temperature_scaling/check_groups.py

Numbers behind nondimensionalization.tex. Material laws from
preprocess/comp_eps.py (the solver's). R = 50 um, eps = 1 um,
alpha_c = 1.020e-3, xi_v = 1e-3. Writes groups.txt beside this file.
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "preprocess"))
import comp_eps as ce                                               # noqa: E402

RHO_I, K_ICE, L_SUB, M_W, KB = 919.0, 2.29, 2.83e6, 2.9915e-26, 1.380649e-23
R, EPS, ALPHA, XI_V = 50e-6, 1e-6, 1.020e-3, 1e-3
A1, A2, D_T = 5.0, 0.1581, 7.78e-6

L = ["T [C]  rho_vs [kg/m3]  D_v [m2/s]   d0 [nm]  beta_sub [s/m]  dbeta/beta   tau = R^2 beta/d0 [d]  tau/tau(-5)  "
     "V = d0/(R beta) [um/d]  L* [um]  Lambda = L*/R   P/xi_v      G_T"]
rows = []
for T in (-40, -30, -20, -10, -5):
    Tk = T + 273.15
    rvs, Dv, d0 = ce.rho_vs_sat(T), ce.Dv_T(T), ce.capillary_length(T)
    bHK = ce.beta_HK(T, ALPHA)                       # sqrt(2 pi m / kB T) / alpha_c
    assert abs(bHK - np.sqrt(2 * np.pi * M_W / (KB * Tk)) / ALPHA) / bHK < 1e-3
    beta = bHK * RHO_I / rvs
    dbeta = A1 * A2 * EPS * (1 / D_T + 1 / Dv) * RHO_I / rvs
    tau = R ** 2 * beta / d0
    V = d0 / (R * beta)
    Ls = Dv * bHK                                    # = D_v beta rho_vs / rho_ice
    P = d0 / (Dv * beta)                             # = R^2 / (D_v tau): vapor storage
    GT = RHO_I * L_SUB ** 2 * M_W * R / (KB * Tk ** 2 * K_ICE * beta)   # latent-heat feedback on saturation
    rows.append((T, rvs, Dv, d0, beta, dbeta / beta, tau, V, Ls, P, GT))
t5 = rows[-1][6]
for T, rvs, Dv, d0, beta, db, tau, V, Ls, P, GT in rows:
    L.append(f"{T:4d}   {rvs:.3e}      {Dv:.3e}   {d0 * 1e9:.3f}    {beta:.3e}      {db:.3f}       {tau / 86400:8.1f}           "
             f"{tau / t5:6.2f}        {V * 1e6 * 86400:.4f}            {Ls * 1e6:5.1f}     {Ls / R:.2f}        {P / XI_V:.1e}    {GT:.1e}")
r = rows
L += ["", f"-40 C against -5 C:  tau {r[0][6] / r[-1][6]:.1f}x   rho_vs {r[-1][1] / r[0][1]:.1f}x   d0 {r[0][3] / r[-1][3]:.3f}x   "
          f"sqrt(T) {np.sqrt(268.15 / 233.15):.3f}x   L* {r[-1][8] / r[0][8]:.3f}x   D_v {r[-1][2] / r[0][2]:.3f}x",
      f"check: tau(-40)/tau(-5) = (rho_vs ratio) x (d0 ratio)^-1... = {r[-1][1] / r[0][1]:.2f} x {np.sqrt(268.15 / 233.15):.3f} / {r[0][3] / r[-1][3]:.3f}"
      f" = {r[-1][1] / r[0][1] * np.sqrt(268.15 / 233.15) / (r[0][3] / r[-1][3]):.1f}",
      f"check: tau_sub = (eps/R)^2 tau (1 + dbeta/beta) at -20 C = {(EPS / R) ** 2 * r[2][6] * (1 + r[2][5]):.0f} s (Table 2: 7822 s)"]
txt = "\n".join(L)
print(txt)
(HERE / "groups.txt").write_text(txt + "\n")
