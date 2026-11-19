"""M1 check (SymPy + NumPy): the symbolic claims behind ## Mechanisms and ## Scales
of spec_problem.md.

  - dimensionless form d theta/d tau = d2 theta/d xi2 + 1 (no parameters)
  - scaling with theta' = T/T0 gives source term Q
  - steady parabola T0 + q x (L-x)/(2k) solves the steady problem with the BCs
  - peak at x = L/2, value q L^2/(8k); theta range [0, 1/8]
  - steady end flux: outward j.n = q L/2 at both ends; sum equals q L
  - slowest mode sin(pi x/L) decays as exp(-pi^2 kappa t/L^2)
  - transient stays in [T0, T0 + qL^2/(8k)] for q >= 0 (maximum-principle range)
  - at tau = t_99 the remaining error relative to the steady peak (informational
    claim "steady state at tau >~ 0.47")

Usage: python m1_nondim.py   (no numeric inputs; all symbolic)
"""
import json
import sys

import numpy as np
import sympy as sp

print(json.dumps({"inputs": "symbolic; parameters L,k,rho,cp,q,T0 positive; "
                  "series truncation N=4001 odd modes for numeric checks"}))
results = []


def emit(name, claim, ok, details):
    results.append(bool(ok))
    print(json.dumps({"check": name, "claim": claim, "pass": bool(ok), "details": details}))


x, t, L, k, rho, cp, q, T0 = sp.symbols("x t L k rho c_p q T_0", positive=True)
xi, tau = sp.symbols("xi tau", positive=True)
kappa = k / (rho * cp)
dTq = q * L**2 / k
Qg = q * L**2 / (k * T0)

# 1. Dimensionless form: T = T0 + dTq*theta(xi,tau), xi=x/L, tau=kappa t/L^2
th = sp.Function("theta")
T_expr = T0 + dTq * th(x / L, kappa * t / L**2)
resid = rho * cp * sp.diff(T_expr, t) - k * sp.diff(T_expr, x, 2) - q
# divide by q (the scale of every term) and substitute back to xi, tau
r = sp.simplify(resid / q)
r = r.subs({x: xi * L, t: tau * L**2 / kappa})
r = sp.simplify(r.doit())
target = sp.Derivative(th(xi, tau), tau) - sp.Derivative(th(xi, tau), xi, 2) - 1
# Evaluate with a test function to avoid Subs/Derivative representation issues
f_test = sp.sin(3 * xi) * sp.exp(-2 * tau) + xi**3 * tau
def apply(expr_fun):
    return expr_fun
T_test = T0 + dTq * f_test.subs({xi: x / L, tau: kappa * t / L**2})
res_test = rho * cp * sp.diff(T_test, t) - k * sp.diff(T_test, x, 2) - q
res_test = sp.simplify((res_test / q).subs({x: xi * L, t: tau * L**2 / kappa}))
tgt_test = sp.diff(f_test, tau) - sp.diff(f_test, xi, 2) - 1
diff1 = sp.simplify(res_test - tgt_test)
emit("dimensionless_form", "rho c_p T_t = k T_xx + q  <=>  theta_tau = theta_xixi + 1 "
     "(divided by q; no parameters remain)", diff1 == 0, {"residual_difference": str(diff1),
     "free_symbols_left": str(sorted(map(str, tgt_test.free_symbols)))})

# BCs / IC map to zeros
emit("dimensionless_bc_ic", "T=T0 at ends and t=0  <=>  theta=0",
     sp.simplify((T0 - T0) / dTq) == 0, {"theta_at_T0": 0})

# 2. theta' = T/T0 form gives +Q
T_test2 = T0 * f_test.subs({xi: x / L, tau: kappa * t / L**2})
res2 = rho * cp * sp.diff(T_test2, t) - k * sp.diff(T_test2, x, 2) - q
# divide by k T0 / L^2 (the conduction scale for theta')
res2 = sp.simplify((res2 / (k * T0 / L**2)).subs({x: xi * L, t: tau * L**2 / kappa}))
tgt2 = sp.diff(f_test, tau) - sp.diff(f_test, xi, 2) - Qg
d2 = sp.simplify(res2 - tgt2)
emit("Q_form", "theta' = T/T0 gives theta'_tau = theta'_xixi + Q with Q = q L^2/(k T0)",
     d2 == 0, {"residual_difference": str(d2)})

# 3. steady parabola
Tss = T0 + q * x * (L - x) / (2 * k)
ode = sp.simplify(k * sp.diff(Tss, x, 2) + q)
bc0 = sp.simplify(Tss.subs(x, 0) - T0)
bcL = sp.simplify(Tss.subs(x, L) - T0)
emit("steady_parabola", "T0 + q x (L-x)/(2k) solves k T'' + q = 0 with T(0)=T(L)=T0",
     ode == 0 and bc0 == 0 and bcL == 0, {"ode": str(ode), "bc0": str(bc0), "bcL": str(bcL)})

xs = sp.solve(sp.diff(Tss, x), x)
peak = sp.simplify(Tss.subs(x, xs[0]) - T0)
emit("peak_location", "peak of steady parabola at x = L/2", xs == [L / 2], {"x_star": str(xs)})
emit("peak_value", "T(L/2) - T0 = q L^2/(8k)", sp.simplify(peak - q * L**2 / (8 * k)) == 0,
     {"peak": str(peak)})
theta_ss = sp.simplify(((Tss - T0) / dTq).subs(x, xi * L))
emit("theta_range", "steady theta in [0, 1/8] on xi in [0,1]",
     sp.simplify(theta_ss.subs(xi, 0)) == 0 and sp.simplify(theta_ss.subs(xi, sp.Rational(1, 2))) == sp.Rational(1, 8),
     {"theta_ss": str(theta_ss)})
emit("Q_over_8", "Delta T_max / T0 = Q/8", sp.simplify(peak / T0 - Qg / 8) == 0, {})

# conduction/source balance: k*DeltaT/L^2 ~ q  => DeltaT ~ qL^2/k
dT = sp.symbols("DeltaT", positive=True)
sol = sp.solve(sp.Eq(k * dT / L**2, q), dT)
emit("balance_scale", "balancing k DeltaT/L^2 with q gives DeltaT ~ q L^2/k",
     sp.simplify(sol[0] - dTq) == 0, {"DeltaT": str(sol)})
# storage vs conduction balance: rho cp DeltaT / t ~ k DeltaT / L^2 => t ~ L^2/kappa
tt = sp.symbols("t_c", positive=True)
solt = sp.solve(sp.Eq(rho * cp * dT / tt, k * dT / L**2), tt)
emit("time_scale", "storage balances conduction on t ~ L^2/kappa",
     sp.simplify(solt[0] - L**2 / kappa) == 0, {"t_c": str(solt)})

# 4. end fluxes: j = -k T', outward j.n with n=-1 at 0, +1 at L
j = -k * sp.diff(Tss, x)
out0 = sp.simplify(j.subs(x, 0) * (-1))
outL = sp.simplify(j.subs(x, L) * (+1))
emit("end_flux", "steady outward flux j.n = q L/2 at each end (sum = q L, the source)",
     sp.simplify(out0 - q * L / 2) == 0 and sp.simplify(outL - q * L / 2) == 0,
     {"out_x0": str(out0), "out_xL": str(outL)})

# 5. slowest mode decay
n = sp.symbols("n", positive=True, integer=True)
lam = sp.symbols("lambda")
mode = sp.exp(-lam * t) * sp.sin(sp.pi * x / L)
eq = sp.simplify((rho * cp * sp.diff(mode, t) - k * sp.diff(mode, x, 2)) / mode)
lam_sol = sp.solve(eq, lam)
emit("slowest_mode_rate", "sin(pi x/L) decays as exp(-pi^2 kappa t/L^2)",
     len(lam_sol) == 1 and sp.simplify(lam_sol[0] - sp.pi**2 * kappa / L**2) == 0,
     {"lambda": str(lam_sol)})

# 6. series solution in theta: theta = theta_ss - sum b_n exp(-n^2 pi^2 tau) sin(n pi xi)
b_n = sp.simplify(2 * sp.integrate(theta_ss * sp.sin(n * sp.pi * xi), (xi, 0, 1)))
b1 = float(b_n.subs(n, 1))
N = 4001
ns = np.arange(1, N + 1)
bn = np.array([float(b_n.subs(n, int(m))) for m in ns[:41]] +
              [0.0] * (N - 41))
# use closed form for speed: b_n = 2(1-(-1)^n)/(n^3 pi^3)
bn_cf = 2 * (1 - (-1.0) ** ns) / (ns**3 * np.pi**3)
cf_ok = np.allclose(bn[:41], bn_cf[:41], rtol=1e-12, atol=1e-15)
emit("series_coeffs", "sine coefficients of theta_ss are nonnegative (b_n = 2(1-(-1)^n)/(n pi)^3)",
     cf_ok and np.all(bn_cf >= 0), {"b_n": str(b_n), "b1": b1})
xg = np.linspace(0, 1, 401)
S = np.sin(np.pi * np.outer(ns, xg))
ok_range = True
worst = None
for tv in [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.3, 1.0]:
    th_t = xg * (1 - xg) / 2 - (bn_cf * np.exp(-(ns * np.pi) ** 2 * tv)) @ S
    lo, hi = th_t.min(), th_t.max()
    if lo < -1e-6 or hi > 0.125 + 1e-12:
        ok_range = False
        worst = (tv, lo, hi)
emit("max_principle_range", "for q >= 0, T(x,t) in [T0, T0 + qL^2/(8k)] (theta in [0,1/8])",
     ok_range, {"violation": worst})

tau99 = np.log(100) / np.pi**2
err = np.max(np.abs((bn_cf * np.exp(-(ns * np.pi) ** 2 * tau99)) @ S))
rel_to_peak = err / 0.125
emit("steady_at_tau99", "at tau = t_99 ~= 0.47 the transient is within ~1 % of the steady peak "
     "(first mode at 1 % of its initial amplitude; higher modes negligible)",
     rel_to_peak <= 0.0105, {"tau99": tau99, "max_abs_theta_err": float(err),
                              "err_over_peak": float(rel_to_peak),
                              "first_mode_amp_ratio": float(np.exp(-np.pi**2 * tau99))})

sys.exit(0 if all(results) else 1)
