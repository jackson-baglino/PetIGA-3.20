import sympy as sp
x, t, xi, s = sp.symbols('x t xi s')
L, k0, rho, c, Q0, T0, DT = sp.symbols('L k0 rho c Q0 T0 DT', positive=True)
theta = sp.Function('theta')
tau = rho*c*L**2/k0
T = T0 + DT*theta(x/L, t/tau)
lhs = rho*c*sp.diff(T, t)
rhs = k0*sp.diff(T, x, 2) + Q0
res = sp.simplify(((lhs - rhs)/(k0*DT/L**2)).subs({x: L*xi, t: tau*s}).doit())
R = sp.Derivative(theta(xi, s), s) - sp.Derivative(theta(xi, s), (xi, 2)) - Q0*L**2/(k0*DT)
print(res)
print("difference:", sp.simplify(res.doit() - R.doit()))
# units: Q0 L^2/(k0 DT): (W m^-3)(m^2)/((W m^-1 K^-1)(K)) = 1
