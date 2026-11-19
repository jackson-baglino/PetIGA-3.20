import sympy as sp

xi, s = sp.symbols('xi s')
L, k0, rho, c, Q0, T0, DT = sp.symbols('L k0 rho c Q0 T0 DT', positive=True)
theta = sp.Function('theta')

tau = rho*c*L**2/k0
# T(x,t) = T0 + DT*theta(x/L, t/tau); derivatives via chain rule
dTdt = DT/tau*sp.Derivative(theta(xi, s), s)
d2Tdx2 = DT/L**2*sp.Derivative(theta(xi, s), (xi, 2))
lhs_minus_rhs = rho*c*dTdt - k0*d2Tdx2 - Q0
factor = rho*c*DT/tau
R = sp.simplify(sp.expand(lhs_minus_rhs/factor))
print("factor =", sp.simplify(factor))
print("R =", R)

# direct chain-rule check with a concrete theta
x, t = sp.symbols('x t')
f = sp.sin(3*xi)*sp.exp(-2*s)
T = T0 + DT*f.subs({xi: x/L, s: t/tau})
res = rho*c*sp.diff(T, t) - k0*sp.diff(T, x, 2) - Q0
Rf = R.subs(theta(xi, s), f).doit()
print("chain-rule check:", sp.simplify(res.subs({x: L*xi, t: tau*s})/factor - Rf))

# units: exponents of (kg, m, s, K); W = kg m^2 s^-3
W = sp.Matrix([1, 2, -3, 0])
k_u = W + sp.Matrix([0, -1, 0, -1])
Q_given = W + sp.Matrix([0, -2, 0, 0])
Q_vol = W + sp.Matrix([0, -3, 0, 0])
Lu = sp.Matrix([0, 1, 0, 0]); Ku = sp.Matrix([0, 0, 0, 1])
print("units of Q0*L**2/(k0*DT) with Q0 in W m^-2:", list(Q_given + 2*Lu - k_u - Ku))
print("units with Q0 in W m^-3:", list(Q_vol + 2*Lu - k_u - Ku))
