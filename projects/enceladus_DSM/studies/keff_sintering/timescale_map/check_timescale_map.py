#!/usr/bin/env python3
"""Every number in timescale_map_derivation.tex, recomputed.

    venv_enceladus/bin/python studies/keff_sintering/timescale_map/check_timescale_map.py [<campaign dir>]

Prints (1) the closed form of tau_sub against the solver's values, (2) the
map's tau_sub against the solver's at the simulated temperatures (they must
agree: the map uses the model's own rho_vs), (3) the effective activation
energy, (4) the crossover length D_v beta_HK, (5) the map's times at a few
(T, R), (6) the grain-size dependence in Molaro et al. (2019), Table 6, and,
if the campaign folder is given, (7) what the target age theta = 331 means
in k_eff and SSA. Writes the same
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


sys.path.insert(0, str(HERE.parents[2] / "preprocess"))
from comp_eps import rho_vs_sat               # noqa: E402  the model's own rho_vs^I(T_C)


def rho_vs(T):                                # [kg/m3] at T [K]
    return np.vectorize(lambda t: rho_vs_sat(float(t) - 273.15))(T)


def tau_closed(T, rho_vs, eps):               # eps^2 beta_sub / d0, written out
    return eps ** 2 * RHO_I ** 2 * np.sqrt(2 * np.pi * M * KB * T) / (ALPHA_C * GAMMA * M * rho_vs)


def rate(T):                                  # 1/tau_sub up to a constant: rho_vs / sqrt(T)
    return rho_vs(T) / np.sqrt(T)


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
p("2. The map's tau_sub = tau_ref * [sqrt(T)/rho_vs(T)] / [sqrt(T_ref)/rho_vs(T_ref)] against the solver's:")
tau_m = TAU_REF * rate(T_REF) / rate(T_SIM)
for T, a, b in zip(T_SIM, tau_m, TAU_SOLVER):
    p(f"   T = {T - 273.15:5.0f} C: map {a:7.0f} s, solver {b:7.0f} s ({100 * (a / b - 1):+.2f} %)")
p()
p("3. Effective activation energy of 1/tau_sub, E = -R d ln[rho_vs/sqrt(T)] / d(1/T):")
for T in (T_REF, 180.0):
    h = 0.5
    e = -RGAS * (np.log(rate(T + h)) - np.log(rate(T - h))) / (1 / (T + h) - 1 / (T - h))
    p(f"   T = {T:.2f} K: E = {float(e) / 1e3:.1f} kJ/mol")
p()
p("4. Crossover length L* = D_v beta_HK (air at 1 atm, alpha_c = 1e-3):")
for T in (T_REF, 180.0):
    dv = DV0 * (T / 273.15) ** 1.81
    bhk = np.sqrt(2 * np.pi * M / (KB * T)) / ALPHA_C
    p(f"   T = {T:.2f} K: D_v = {dv:.3e} m2/s, beta_HK = {bhk:.3f} s/m, L* = {dv * bhk * 1e6:.0f} um")
p()
p(f"5. Time to theta = {THETA:g} (= 30 d / tau_sub at -20 C = {30 * DAY / TAU_REF:.1f}):")
for T in (253.15, 220, 200, 180, 150, 120, 110, 80, 50):
    p(f"   T = {T:6.1f} K: " + " | ".join(f"R = {R * 1e6:g} um: {fmt(t_map(T, R))}" for R in (1e-6, 5e-6, 6e-6, 50e-6)))

p()
p("6. Grain-size dependence in Molaro et al. (2019), Table 6 (tau = f T^g, tau in years):")
TAB6 = {0.1: (3.189e39, -19.18), 1: (8.945e77, -35.54), 10: (1.204e88, -39.26),
        30: (3.743e89, -39.51), 50: (1.043e90, -39.51), 100: (4.159e90, -39.51)}
mol = lambda T, r: TAB6[r][0] * T ** TAB6[r][1]
for T in (130.0, 150.0, 180.0):
    ex = [np.log(mol(T, b) / mol(T, a)) / np.log(b / a) for a, b in ((1, 10), (10, 30), (30, 50), (50, 100))]
    p(f"   T = {T:.0f} K: d ln tau / d ln R = " + ", ".join(f"{e:.2f}" for e in ex)
      + "  (1-10, 10-30, 30-50, 50-100 um)")
p("   Molaro's stage-1 timescale against this map (theta = 331), at 180 K:")
for r in (1, 10, 50):
    p(f"   R = {r:3d} um: Molaro {fmt(mol(180.0, r) * YEAR)}; this map {fmt(t_map(180.0, r * 1e-6))}")

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
    p(f"7. At theta = {THETA:g}, {len(kk)} packings with porosity <= 0.375 (-10 C runs):")
    p(f"   k_eff / k_eff,r = {np.mean(kk):.3f} +- {np.std(kk, ddof=1):.3f};  SSA / SSA_r = {np.mean(ss):.3f} +- {np.std(ss, ddof=1):.3f}")
(HERE / "check_timescale_map.txt").write_text("\n".join(out) + "\n")
