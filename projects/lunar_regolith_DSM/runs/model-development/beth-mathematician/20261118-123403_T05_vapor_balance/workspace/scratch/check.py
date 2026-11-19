import sympy as sp
x, y, z, t = sp.symbols('x y z t')
X = (x, y, z)
phi = sp.Function('phi')(*X, t)
T = sp.Function('T')(*X, t)
rv = sp.Function('rho_v')(*X, t)
rho_i = sp.Symbol('rho_i', positive=True)
D = sp.Function('D')

# S1-S3: total water mass density m and flux J
m = rho_i*phi + (1 - phi)*rv
J = [-(1 - phi)*D(T)*sp.diff(rv, xi) for xi in X]
balance = sp.diff(m, t) + sum(sp.diff(J[i], X[i]) for i in range(3))

R = (sp.Derivative((1 - phi)*rv, t) + rho_i*sp.Derivative(phi, t)
     - sum(sp.Derivative((1 - phi)*D(T)*sp.Derivative(rv, xi), xi) for xi in X))
print("S4 balance - R =", sp.simplify(balance - R.doit()))

# Expanded form check: (1-phi) d_t rho_v + (rho_i - rho_v) d_t phi - div(...)
alt = (1 - phi)*sp.diff(rv, t) + (rho_i - rv)*sp.diff(phi, t) \
    - sum(sp.diff((1 - phi)*D(T)*sp.diff(rv, xi), xi) for xi in X)
print("S5 expanded - R =", sp.simplify(alt - R.doit()))

# 1D integral check: with zero-flux ends, total mass is conserved
# (divergence integrates to boundary flux) -- verified symbolically in 1D
f = sp.Function('f')
print("1D div integral:", sp.integrate(sp.diff(f(x), x), (x, 0, 1)))

# units
import pint
u = pint.UnitRegistry()
rv_u = 1*u.kg/u.m**3; rhoi_u = 1*u.kg/u.m**3; D_u = 1*u.m**2/u.s
terms = {
    "d_t((1-phi) rho_v)": rv_u/u.s,
    "rho_i d_t phi": rhoi_u/u.s,
    "d_x((1-phi) D d_x rho_v)": D_u*rv_u/u.m**2,
}
for k, v in terms.items():
    print(k, v.to_base_units().units)
