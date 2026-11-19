"""Scratch copy of derivation/mms/mms_sources.py (same formulas, verbatim body), because the guard runs
only scripts under derivation/scratch/. My own check, not evidence for anyone else."""
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

checks = {
    "theta'_m(0,tau) - 1": sp.simplify(theta_m.subs(xi, 0) - 1),
    "theta'_m(1,tau) - 1": sp.simplify(theta_m.subs(xi, 1) - 1),
    "theta'_m(xi,0) - 1": sp.simplify(theta_m.subs(tau, 0) - 1),
    "corner (0,0): d_xixi theta'_m + Q + f": sp.simplify((sp.diff(theta_m, xi, 2) + Q + f).subs({xi: 0, tau: 0})),
    "corner (1,0): d_xixi theta'_m + Q + f": sp.simplify((sp.diff(theta_m, xi, 2) + Q + f).subs({xi: 1, tau: 0})),
    "PDE residual (all xi,tau)": sp.simplify(sp.diff(theta_m, tau) - sp.diff(theta_m, xi, 2) - Q - f),
    "d_tau theta'_m at ends (should be 0)": sp.simplify(sp.diff(theta_m, tau).subs(xi, 0)),
}
for name, val in checks.items():
    print("%-40s %s" % (name, val))
print("A = B = 0 gives f =", f.subs({A: 0, B: 0}), " (cancels the source: theta' stays 1)")
print()

kappa = k / (rho * cp)
q_m = (k * T0 / L**2) * f.subs({xi: x / L, tau: kappa * t / L**2})
print("q_m(x,t) [W m^-3] = (k T0/L^2) f(x/L, kappa t/L^2) =", sp.simplify(q_m))
print()

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
