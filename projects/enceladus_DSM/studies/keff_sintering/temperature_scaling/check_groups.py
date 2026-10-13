#!/usr/bin/env python3
"""check_groups.py — numbers behind time_temperature_scaling.tex.

    venv_enceladus/bin/python studies/keff_sintering/temperature_scaling/check_groups.py

Evaluates every scale and dimensionless group of the scaling analysis at the
five campaign temperatures, with the solver's material laws
(preprocess/comp_eps.py), and checks them against independent values:
the closed form of tau, the tau_sub of manuscript Table 2, and the
activation energy. R = 50 um, eps = 1 um, alpha_c = 1.020e-3, xi_v = 1e-3.
Writes groups.txt beside this file.
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "preprocess"))
import comp_eps as ce                                               # noqa: E402

RHO_I, K_ICE, L_SUB, M_W, KB, GAMMA = 919.0, 2.29, 2.83e6, 2.9915e-26, 1.380649e-23, 0.109
R, EPS, ALPHA, XI_V = 50e-6, 1e-6, 1.020e-3, 1e-3
A1, A2, D_T = 5.0, 0.1581, 7.78e-6
R_GAS, M_MOL = 8.314462, 0.0180153
TABLE2_TAU_SUB = {-40: 60390, -30: 20840, -20: 7822, -10: 3164, -5: 2063}
T_REF = -20

rows = {}
for T in (-40, -30, -20, -10, -5):
    Tk = T + 273.15
    rvs, Dv, d0 = ce.rho_vs_sat(T), ce.Dv_T(T), ce.capillary_length(T)
    bHK = np.sqrt(2 * np.pi * M_W / (KB * Tk)) / ALPHA
    assert abs(ce.beta_HK(T, ALPHA) / bHK - 1) < 1e-3 and abs(d0 / (GAMMA * M_W / (RHO_I * KB * Tk)) - 1) < 1e-3
    b_sub = bHK * RHO_I / rvs
    delta = A1 * A2 * EPS * (1 / D_T + 1 / Dv) / bHK              # dbeta / beta_sub
    beta = b_sub * (1 + delta)
    tau = R ** 2 * beta / d0
    tau_closed = (1 + delta) * R ** 2 * RHO_I ** 2 * np.sqrt(2 * np.pi * KB * Tk / M_W) / (ALPHA * GAMMA * rvs)
    rows[T] = dict(Tk=Tk, rvs=rvs, Dv=Dv, d0=d0, b_sub=b_sub, delta=delta, beta=beta, tau=tau, tau_closed=tau_closed,
                   V=d0 / (R * beta), Ls=Dv * beta * rvs / RHO_I, P=d0 / (Dv * beta),
                   GT=RHO_I * L_SUB ** 2 * M_W * R / (KB * Tk ** 2 * K_ICE * beta), sig=d0 / R,
                   dT=RHO_I * L_SUB * (d0 / (R * beta)) * R / K_ICE, tau_sub=(EPS / R) ** 2 * tau)
ref = rows[T_REF]
L = ["T [C]   rho_vs      D_v        d0 [nm]  beta_sub   delta    beta       d0/R       V [um/d]  tau [d]   a_T      L* [um]  Lambda  P/xi_v    G_T      dT [K]    tau_sub [s] (Table 2)"]
for T, r in rows.items():
    L.append(f"{T:4d}   {r['rvs']:.3e}  {r['Dv']:.3e}  {r['d0'] * 1e9:.3f}    {r['b_sub']:.3e}  {r['delta']:.4f}   {r['beta']:.3e}  {r['sig']:.3e}  "
             f"{r['V'] * 1e6 * 86400:.4f}    {r['tau'] / 86400:7.1f}   {r['tau'] / ref['tau']:6.3f}   {r['Ls'] * 1e6:6.1f}   {r['Ls'] / R:.3f}   "
             f"{r['P'] / XI_V:.1e}   {r['GT']:.1e}  {r['dT']:.1e}   {r['tau_sub']:7.0f} ({TABLE2_TAU_SUB[T]})")
    assert abs(r["tau_closed"] / r["tau"] - 1) < 2e-3, "closed form of tau"
    assert abs(r["tau_sub"] / TABLE2_TAU_SUB[T] - 1) < 2e-3, "tau_sub against Table 2"
a, b = rows[-40], rows[-5]
Tk = np.array([r["Tk"] for r in rows.values()]); tau = np.array([r["tau"] for r in rows.values()])
Ea = R_GAS * np.polyfit(1 / Tk, np.log(tau), 1)[0] / 1e3
L += ["",
      f"-40 C against -5 C: tau {a['tau'] / b['tau']:.2f}x = rho_vs {b['rvs'] / a['rvs']:.2f}x  *  sqrt(T) {np.sqrt(a['Tk'] / b['Tk']):.4f}x  *  (1+delta) {(1 + a['delta']) / (1 + b['delta']):.4f}x"
      f"  = {b['rvs'] / a['rvs'] * np.sqrt(a['Tk'] / b['Tk']) * (1 + a['delta']) / (1 + b['delta']):.2f}",
      f"                    L* {b['Ls'] / a['Ls']:.3f}x;  D_v {b['Dv'] / a['Dv']:.3f}x;  (T ratio)^1.31 = {(b['Tk'] / a['Tk']) ** 1.31:.3f}",
      f"activation energy of 1/tau, fitted over the five temperatures: {Ea:.1f} kJ/mol;  L_sub M_w - R T/2 at -20 C = {(L_SUB * M_MOL - 0.5 * R_GAS * 253.15) / 1e3:.1f} kJ/mol",
      "all assertions passed: closed form of tau; tau_sub against Table 2 (0.2 %); beta_HK and d0 against comp_eps"]
txt = "\n".join(L)
print(txt)
(HERE / "groups.txt").write_text(txt + "\n")
