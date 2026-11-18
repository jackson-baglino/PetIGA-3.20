"""The parameter map between the derived phase-field equation (T07) and the code's A1 row, checked with SymPy.

    venv_workflows/bin/python studies/derived_model/verification/coupling_map.py

T07 (beth-mathematician, run 20261118-124517):
    phi_t = -M [6 gamma/eps phi(1-phi)(1-2phi) - 6 gamma eps lap phi] + 30 M rho_i R_v T phi^2(1-phi)^2 (rho_v - rho_vs)/rho_vs
A1 (src/assembly.c, Residual_A1, phi row):
    phi_t = -[3 mob/eps phi(1-phi)(1-2phi) - 3 mob eps lap phi] + (alpha/rho_i) phi^2(1-phi)^2 (rho_v - rho_vs)
with, in src/lunar_main.c (setup of the kinetic parameters):
    mob = eps/(3 tau), alpha = lambda/tau, lambda = a1 eps/d0_sub, a1 = 5,
    d0_sub = d0 / (rho_i/rho_vs0), d0 = gamma (m_w/rho_i)/(k_B T0)   [D0_CAPILLARY, include/NASA_types.h]
Prints each identity and exits 1 if one fails.
"""
import sys

import sympy as sp

gamma, eps, tau, rho_i, rho_vs0, rho_vs, T0, T, k_B, m_w, M, a1 = sp.symbols(
    "gamma eps tau rho_i rho_vs0 rho_vs T0 T k_B m_w M a1", positive=True)
R_v = k_B / m_w

mob = eps / (3 * tau)
d0 = gamma * (m_w / rho_i) / (k_B * T0)
d0_sub = d0 / (rho_i / rho_vs0)
lam = a1 * eps / d0_sub
alpha = lam / tau

M_map = mob / (2 * gamma)           # the map used by -residual_model derived (src/lunar_main.c)
derived_coupling = 30 * M * rho_i * R_v * T / rho_vs

checks = {
    "gradient term: 6 gamma M eps = 3 mob eps": sp.simplify(6 * gamma * M_map * eps - 3 * mob * eps),
    "well term: 6 gamma M / eps = 3 mob / eps": sp.simplify(6 * gamma * M_map / eps - 3 * mob / eps),
    "coupling at T = T0, rho_vs = rho_vs0, a1 = 5: 30 M rho_i R_v T/rho_vs = alpha/rho_i":
        sp.simplify((derived_coupling.subs({M: M_map, T: T0, rho_vs: rho_vs0}) - alpha / rho_i).subs(a1, 5)),
}
ratio = sp.simplify(derived_coupling.subs(M, M_map) / (alpha / rho_i))
ok = True
for name, diff in checks.items():
    good = diff == 0
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}  {name}   (difference: {diff})")
print(f"general ratio derived/A1 coupling = {ratio}   (1 at T = T0 and rho_vs = rho_vs0 when a1 = 5)")
alpha_sym = sp.Symbol("alpha", positive=True)
print("alpha_sub in terms of the derived M, frozen at T0: alpha_sub =", sp.simplify(sp.solve(sp.Eq(alpha_sym / rho_i, 30 * M * rho_i * R_v * T0 / rho_vs0), alpha_sym)[0]))

# Heat capacity per unit volume at phi = 1/2, code's values (src/lunar_main.c defaults).
vals = {"rho_i": 919.0, "rho_a": 1.341, "c_i": 1.96e3, "c_a": 1.044e3}
phi = sp.Rational(1, 2)
mix = phi * vals["rho_i"] * vals["c_i"] + (1 - phi) * vals["rho_a"] * vals["c_a"]
prod = (phi * vals["rho_i"] + (1 - phi) * vals["rho_a"]) * (phi * vals["c_i"] + (1 - phi) * vals["c_a"])
print(f"heat capacity at phi = 1/2: mixture {float(mix):.6g} J m^-3 K^-1, product {float(prod):.6g}, "
      f"product/mixture = {float(prod / mix):.4f}")
p_, ri, ra, ci, ca = sp.symbols("phi rho_i rho_a c_i c_a")
gap = sp.factor(sp.expand((p_ * ri + (1 - p_) * ra) * (p_ * ci + (1 - p_) * ca) - (p_ * ri * ci + (1 - p_) * ra * ca)))
print(f"product - mixture = {gap}")
sys.exit(0 if ok else 1)
