import sympy as sp
x, g, e, H, k = sp.symbols('x gamma eps H kappa', positive=True)
p = sp.symbols('phi')
W = H*p**2*(1-p)**2
# S2: wells and barrier
print('W(0),W(1):', W.subs(p,0), W.subs(p,1))
print("W'(0),W'(1):", sp.diff(W,p).subs(p,0), sp.diff(W,p).subs(p,1))
print("W''(0),W''(1):", sp.diff(W,p,2).subs(p,0), sp.diff(W,p,2).subs(p,1))
print('W(1/2):', W.subs(p, sp.Rational(1,2)))
# S4: profile solves first integral (k/2) phi'^2 = W and EL k phi'' = W'
d = sp.symbols('delta', positive=True)
prof = (1+sp.tanh(x/(2*d)))/2
dp = sp.diff(prof, x)
fi = sp.simplify((k/2*dp**2 - W.subs(p, prof)).subs(d, sp.sqrt(k/(2*H))).rewrite(sp.exp))
print('first integral residual:', fi)
el = sp.simplify((k*sp.diff(prof,x,2) - sp.diff(W,p).subs(p,prof)).subs(d, sp.sqrt(k/(2*H))).rewrite(sp.exp))
print('EL residual:', el)
# S5: energy
gam = sp.integrate(sp.sqrt(2*k*W), (p, 0, 1))
print('gamma(H,k) =', sp.simplify(gam))
# S6: solve
sol = sp.solve([sp.Eq(sp.sqrt(k/(2*H)), e), sp.Eq(sp.sqrt(2*k*H)/6, g)], [H, k], dict=True)
print('solution:', sol)
Hs, ks = 3*g/e, 6*g*e
# S7: direct integration of f over profile
pr = prof.subs(d, e)
fdens = Hs*pr**2*(1-pr)**2 + ks/2*sp.diff(pr,x)**2
u = sp.symbols('u')
val = sp.integrate(sp.simplify(fdens.rewrite(sp.exp)), (x, -sp.oo, sp.oo))
print('integral of f over x:', sp.simplify(val))
# width
print('max slope at x=0:', sp.diff(pr,x).subs(x,0), ' -> width 1/slope =', 1/sp.diff(pr,x).subs(x,0))
print('10-90 width / eps:', sp.N(2*2*sp.atanh(sp.Rational(4,5))))
print('phi at +-2eps:', sp.N(pr.subs(x,-2*e)), sp.N(pr.subs(x,2*e)))
# units with pint
import pint
U = pint.UnitRegistry()
gq, eq = 1*U.J/U.m**2, 1*U.m
print('units W term:', (3*gq/eq).to_base_units().units, ' grad term:', (3*gq*eq/U.m**2).to_base_units().units, ' target:', (1*U.J/U.m**3).to_base_units().units)
