"""Manufactured solution for heat1d-fd-be (DEV TOY DRY RUN ONLY; not a result, do not cite).

Prints the source term f(xi, tau) that must be added to every residual row (spec.md, ## Manufactured
solutions) so that theta'_m below is an exact solution of the forced problem

    d_tau theta' = d_xixi theta' + Q + f(xi, tau),  0 < xi < 1,
    theta'(0, tau) = theta'(1, tau) = 1,  theta'(xi, 0) = 1,

which is S29 (ledger) with the forcing f added. It also checks, symbolically, that theta'_m meets the
production boundary and initial data, so the MMS option changes nothing but the forcing, and that the
data are compatible at the corners (no corner singularity, full order expected).

Run:  <python> -I mms_sources.py [A B Q]      (defaults A = 0.5, B = 0.25, Q = 1)
"""
import sys

import sympy as sp

xi, tau = sp.symbols('xi tau', real=True)
A, B, Q = sp.symbols('A B Q', real=True)
x, t = sp.symbols('x t', real=True)
L, k, rho, cp, T0 = sp.symbols('L k rho c_p T_0', positive=True)

theta_m = 1 + A * sp.sin(sp.pi * xi) * (1 - sp.exp(-tau)) + B * xi * (1 - xi) * sp.sin(2 * sp.pi * tau)
f = sp.simplify(sp.diff(theta_m, tau) - sp.diff(theta_m, xi, 2) - Q)

print("theta'_m(xi,tau) =", theta_m)
print("f(xi,tau)        =", sp.expand(f))
print("python: f = %s" % sp.pycode(sp.expand(f)))
print()

# --- checks against the production data ------------------------------------
checks = {
    "theta'_m(0,tau) - 1": sp.simplify(theta_m.subs(xi, 0) - 1),
    "theta'_m(1,tau) - 1": sp.simplify(theta_m.subs(xi, 1) - 1),
    "theta'_m(xi,0) - 1": sp.simplify(theta_m.subs(tau, 0) - 1),
    # corner compatibility: the end condition theta' = 1 forces d_tau theta' = 0 at the corners, so the
    # right side d_xixi theta' + Q + f must vanish there at tau = 0 (it does not for f = 0: that is the
    # corner incompatibility of the production problem, selection.md item 1)
    "corner (0,0): d_xixi theta'_m + Q + f": sp.simplify((sp.diff(theta_m, xi, 2) + Q + f).subs({xi: 0, tau: 0})),
    "corner (1,0): d_xixi theta'_m + Q + f": sp.simplify((sp.diff(theta_m, xi, 2) + Q + f).subs({xi: 1, tau: 0})),
    "PDE residual (all xi,tau)": sp.simplify(sp.diff(theta_m, tau) - sp.diff(theta_m, xi, 2) - Q - f),
    "d_tau theta'_m at ends (should be 0)": sp.simplify(sp.diff(theta_m, tau).subs(xi, 0)),
}
for name, val in checks.items():
    print("%-40s %s" % (name, val))
print("A = B = 0 gives f =", f.subs({A: 0, B: 0}), " (cancels the source: theta' stays 1)")
print()

# --- dimensional form (only for a solver that works in SI) ------------------
# T = T0 theta', x = L xi, t = L^2 tau / kappa. The S28 -> S29 step divides by kappa T0 / L^2 and
# multiplies q by L^2/(k T0); so a dimensionless row forcing f corresponds to a volumetric source
# q_m = (k T0 / L^2) f  [W m^-3].
kappa = k / (rho * cp)
q_m = (k * T0 / L**2) * f.subs({xi: x / L, tau: kappa * t / L**2})
print("q_m(x,t) [W m^-3] = (k T0/L^2) f(x/L, kappa t/L^2) =", sp.simplify(q_m))
print()

# --- numeric samples for cross-checking an implementation ------------------
vals = [0.5, 0.25, 1.0]
if len(sys.argv) == 4:
    vals = [float(s) for s in sys.argv[1:4]]
sub = {A: vals[0], B: vals[1], Q: vals[2]}
print("samples with A=%g B=%g Q=%g:" % tuple(vals))
print("   %6s %6s %22s %22s" % ("xi", "tau", "theta'_m", "f"))
for xv in (0.0, 0.25, 0.5, 0.75, 1.0):
    for tv in (0.0, 0.1, 0.5):
        tm = float(theta_m.subs(sub).subs({xi: xv, tau: tv}))
        fv = float(f.subs(sub).subs({xi: xv, tau: tv}))
        print("   %6.3f %6.3f %22.15e %22.15e" % (xv, tv, tm, fv))
