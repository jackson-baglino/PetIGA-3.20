import sympy as sp
x = sp.symbols('x', real=True)
g, e, A, k, p = sp.symbols('gamma eps A kappa p', positive=True)
W = p**2*(1-p)**2
# S2: stationary points of W
print("W' roots:", sp.solve(sp.diff(W, p), p), " W'' at 0,1:", [sp.diff(W, p, 2).subs(p, v) for v in (0, 1)])
# S4: profile with delta = sqrt(kappa/(2A))
d = sp.sqrt(k/(2*A))
ph = (1 + sp.tanh(x/(2*d)))/2
EL = k*sp.diff(ph, x, 2) - A*sp.diff(W, p).subs(p, ph)
print("EL residual:", sp.simplify(EL.rewrite(sp.exp)))
FI = k/2*sp.diff(ph, x)**2 - A*W.subs(p, ph)
print("first integral residual:", sp.simplify(FI.rewrite(sp.exp)))
# S5: gamma = int kappa phi'^2 dx = int sqrt(2 A kappa) p(1-p) dp
gam = sp.integrate(sp.sqrt(2*A*k)*p*(1-p), (p, 0, 1))
print("gamma(A,kappa):", gam)
# S6: choose delta = eps, gamma = given
sol = sp.solve([sp.Eq(d, e), sp.Eq(gam, g)], [A, k], dict=True)
print("A, kappa:", sol)
# direct numerical check of energy with final coefficients
Af, kf = 3*g/e, 3*g*e  # kappa/2 = 3 g e -> kappa = 6 g e
phf = (1 + sp.tanh(x/(2*e)))/2
dens = Af*W.subs(p, phf) + kf*sp.diff(phf, x)**2
val = sp.integrate(dens.subs({g: 1, e: sp.Rational(1, 3)}), (x, -sp.oo, sp.oo))
print("energy per area with gamma=1, eps=1/3:", sp.nsimplify(sp.N(val, 30)))
# width: 1/max slope and 10-90 width
s0 = sp.diff(phf, x).subs(x, 0)
print("1/phi'(0):", sp.simplify(1/s0))
x90 = sp.solve(sp.Eq(phf, sp.Rational(9, 10)), x)
print("10-90 width:", [sp.simplify(2*v) for v in x90])
# units: [A]=J m^-3, [kappa]=J m^-1, [kappa |grad phi|^2] = J m^-1 m^-2 = J m^-3
import pint
u = pint.UnitRegistry()
gq, eq = 1*u.J/u.m**2, 1*u.m
print((gq/eq).to_base_units().units, (gq*eq/eq**2).to_base_units().units)
