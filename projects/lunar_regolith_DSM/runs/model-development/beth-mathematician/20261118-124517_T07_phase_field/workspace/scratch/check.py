import sympy as sp
from sympy.calculus.euler import euler_equations
x, y, z, t = sp.symbols('x y z t')
g, e = sp.symbols('gamma eps', positive=True)
phi = sp.Function('phi')(x, y, z, t)
f = 3*g/e*phi**2*(1-phi)**2 + 3*g*e*(phi.diff(x)**2 + phi.diff(y)**2 + phi.diff(z)**2)
# S1: variational derivative (Euler-Lagrange expression)
EL = sp.diff(f, phi) - sum(sp.diff(sp.diff(f, phi.diff(v)), v) for v in (x, y, z))
lap = phi.diff(x, 2) + phi.diff(y, 2) + phi.diff(z, 2)
dF = 6*g/e*phi*(1-phi)*(1-2*phi) - 6*g*e*lap
print('S1 dF/dphi check:', sp.simplify(EL - dF))
# S3: h(phi) properties
p = sp.symbols('p')
h = p**3*(10 - 15*p + 6*p**2)
hp = sp.diff(h, p)
print('S3 h(0),h(1):', h.subs(p, 0), h.subs(p, 1), ' h-30p^2(1-p)^2:', sp.factor(hp - 30*p**2*(1-p)**2))
# S2/S4 equilibrium profile and integral of phi'^2
X = sp.symbols('X', real=True)
phi0 = (1 + sp.tanh(X/(2*e)))/2
res = 6*g/e*phi0*(1-phi0)*(1-2*phi0) - 6*g*e*phi0.diff(X, 2)
print('S2 profile solves dF=0:', sp.simplify(res.rewrite(sp.exp)))
I2 = sp.integrate(sp.simplify(phi0.diff(X)**2).rewrite(sp.exp), (X, -sp.oo, sp.oo))
print('S8 int phi0_x^2 =', sp.simplify(I2), ' expect 1/(6 eps)')
gam = sp.integrate((2*3*g/e*(phi0*(1-phi0))**2).rewrite(sp.exp), (X, -sp.oo, sp.oo))
print('S2 interface energy =', sp.simplify(gam))
w = 2*e*sp.log(9)  # 10-90 width: phi0=0.9 at X = 2 eps ln 9
print('S2 phi0(eps ln9)=', sp.nsimplify(sp.simplify(phi0.subs(X, e*sp.log(9)))))
# units with pint
import pint
u = pint.UnitRegistry()
gamma = 1*u.J/u.m**2; eps = 1*u.m; M = 1*u.m**3/u.J/u.s
rho_i = 1*u.kg/u.m**3; Rv = 1*u.J/u.kg/u.K; T = 1*u.K; rho = 1*u.kg/u.m**3
terms = [1/u.s, M*gamma/eps, M*gamma*eps/u.m**2, M*rho_i*Rv*T*rho/rho]
print('units:', [str(tt.to_base_units().units) for tt in terms])
# Hertz-Knudsen typical M at 200 K
import math
Rv_=461.5; T_=200; psat=0.163; rvs=psat/(Rv_*T_); ri=917; ep=1e-6
Mtyp = rvs*math.sqrt(Rv_*T_/(2*math.pi))/(6*ep*ri**2*Rv_*T_)
print('rho_vs(200K)=', rvs, ' M typical=', Mtyp)
