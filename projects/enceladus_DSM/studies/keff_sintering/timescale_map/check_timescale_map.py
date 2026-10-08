#!/usr/bin/env python3
"""Every number in timescale_map_derivation.tex, recomputed.

    venv_enceladus/bin/python studies/keff_sintering/timescale_map/check_timescale_map.py [<campaign dir>]

Prints (1) the closed form of tau_sub against the solver's values, (2) the
map's temperature scaling against the solver's over the simulated range,
(3) the effective activation energy, (4) the crossover length D_v beta_HK,
(5) the map's times at a few (T, R), and, if the campaign folder is given,
(6) what the target age theta = 331 means in k_eff and SSA. Writes the same
text to check_timescale_map.txt next to this file.

The map itself is fig_timescales() in ../figures/sample_figures.py; the
formulas here are copied from it and from preprocess/comp_eps.py on purpose,
so a change there shows up as a disagreement here.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
M, KB, RHO_I, GAMMA, ALPHA_C = 2.9915e-26, 1.380649e-23, 919.0, 0.109, 1e-3
DV0, RGAS = 2.178e-5, 8.314
YEAR, DAY = 365.25 * 86400, 86400.0
TAU_REF, T_REF, R_REF, EPS_REF, THETA = 7822.3, 253.15, 50e-6, 1e-6, 331.0

# the solver's own values (outp.txt banners; Table "keff_temperature")
T_SIM = np.array([-40, -30, -20, -10, -5.0]) + 273.15
RHO_SOLVER = np.array([1.055e-4, 3.123e-4, 8.487e-4, 2.139e-3, 3.311e-3])
TAU_SOLVER = np.array([60390, 20840, 7822.3, 3164, 2063.0])


def psat(T):                                  # Murphy & Koop (2005), Pa
    return np.exp(9.550426 - 5723.265 / T + 3.53068 * np.log(T) - 0.00728332 * T)


def tau_closed(T, rho_vs, eps):               # eps^2 beta_sub / d0, written out
    return eps ** 2 * RHO_I ** 2 * np.sqrt(2 * np.pi * M * KB * T) / (ALPHA_C * GAMMA * M * rho_vs)


def rate(T):                                  # 1/tau_sub up to a constant, ideal-gas rho_vs = p m / (k T)
    return psat(T) / T ** 1.5


def t_map(T, R):
    return THETA * TAU_REF * rate(T_REF) / rate(T) * (R / R_REF) ** 2


def fmt(t):
    for u, s in ((YEAR * 1e9, "Gyr"), (YEAR * 1e6, "Myr"), (YEAR * 1e3, "kyr"), (YEAR, "yr"), (DAY, "d"), (3600, "h")):
        if t >= u:
            return f"{t / u:.3g} {s}"
    return f"{t:.3g} s"


out = []
p = lambda s="": (print(s), out.append(s))
p("1. tau_sub = eps^2 rho_ice^2 sqrt(2 pi m k T) / (alpha_c gamma m rho_vs), with the solver's rho_vs:")
for T, a, b in zip(T_SIM, tau_closed(T_SIM, RHO_SOLVER, EPS_REF), TAU_SOLVER):
    p(f"   T = {T - 273.15:5.0f} C: closed form {a:8.0f} s, solver {b:8.0f} s ({100 * (a / b - 1):+.2f} %)")
p()
p("2. The map's temperature scaling (Murphy-Koop p_sat, ideal gas) against the solver's (ASHRAE, fixed air density):")
rho_mk = psat(T_SIM) * M / (KB * T_SIM)
tau_m = TAU_REF * rate(T_REF) / rate(T_SIM)
for T, r, a, b in zip(T_SIM, rho_mk / RHO_SOLVER, tau_m, TAU_SOLVER):
    p(f"   T = {T - 273.15:5.0f} C: rho_vs ratio {r:.3f}; tau_sub map {a:7.0f} s, solver {b:7.0f} s ({100 * (a / b - 1):+.1f} %)")
p()
p("3. Effective activation energy of 1/tau_sub, E = Q_sub - (3/2) R T:")
for T in (T_REF, 180.0):
    q = -(-5723.265 - 3.53068 * T + 0.00728332 * T ** 2)
    p(f"   T = {T:.2f} K: Q_sub = {q * RGAS / 1e3:.1f} kJ/mol, E = {(q - 1.5 * T) * RGAS / 1e3:.1f} kJ/mol")
p()
p("4. Crossover length L* = D_v beta_HK (air at 1 atm, alpha_c = 1e-3):")
for T in (T_REF, 180.0):
    dv = DV0 * (T / 273.15) ** 1.81
    bhk = np.sqrt(2 * np.pi * M / (KB * T)) / ALPHA_C
    p(f"   T = {T:.2f} K: D_v = {dv:.3e} m2/s, beta_HK = {bhk:.3f} s/m, L* = {dv * bhk * 1e6:.0f} um")
p()
p(f"5. Time to theta = {THETA:g} (= 30 d / tau_sub at -20 C = {30 * DAY / TAU_REF:.1f}):")
for T in (253.15, 220, 200, 180, 150, 120, 110, 80, 60):
    p(f"   T = {T:6.1f} K: " + " | ".join(f"R = {R * 1e6:g} um: {fmt(t_map(T, R))}" for R in (1e-6, 5e-6, 6e-6, 50e-6)))

if len(sys.argv) > 1:
    sys.path.insert(0, str(HERE.parents[2] / "postprocess"))
    from plot_keff import load, read_tau_sub
    kk, ss = [], []
    pat = "packing_2D_phi0.[23]*_Rave50um_LR40_seed1[678]0[1-5]_L2mm_eps1000nm_perxy_T-10__snow_T-10_h1.00_30d"
    for d in sorted(Path(sys.argv[1]).glob(pat)):
        r = load(d); th = r["t"] / read_tau_sub(d)
        kk.append(np.interp(THETA, th, r["kiso"]) / np.interp(30, th, r["kiso"]))
        ss.append(np.interp(THETA, th, r["ssa"]) / np.interp(30, th, r["ssa"]))
    p()
    p(f"6. At theta = {THETA:g}, {len(kk)} packings with porosity <= 0.375 (-10 C runs):")
    p(f"   k_eff / k_eff,r = {np.mean(kk):.3f} +- {np.std(kk, ddof=1):.3f};  SSA / SSA_r = {np.mean(ss):.3f} +- {np.std(ss, ddof=1):.3f}")
(HERE / "check_timescale_map.txt").write_text("\n".join(out) + "\n")
