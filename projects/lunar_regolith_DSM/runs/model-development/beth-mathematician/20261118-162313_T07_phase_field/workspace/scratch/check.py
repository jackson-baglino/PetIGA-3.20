import sympy as sp
from sympy.calculus.euler import euler_equations

x, y, z, t = sp.symbols('x y z t', real=True)
g, e = sp.symbols('gamma eps', positive=True)
phi = sp.Function('phi')(x, y, z, t)
X = sp.symbols('X', real=True)


def check_s1():
    dens = 3*g/e*phi**2*(1-phi)**2 + 3*g*e*(sp.diff(phi, x)**2 + sp.diff(phi, y)**2 + sp.diff(phi, z)**2)
    eq = euler_equations(dens, phi, [x, y, z, t])[0].lhs  # = d f/d phi - div(...)
    claim = 6*g/e*phi*(1-phi)*(1-2*phi) - 6*g*e*(sp.diff(phi, x, 2) + sp.diff(phi, y, 2) + sp.diff(phi, z, 2))
    print('S1 dF/dphi residual:', sp.simplify(eq - claim))


def check_s2():
    p = sp.symbols('p')
    h = p**3*(10 - 15*p + 6*p**2)
    print("S2 h' - 30 p^2(1-p)^2:", sp.factor(sp.diff(h, p) - 30*p**2*(1-p)**2),
          ' h(0), h(1):', h.subs(p, 0), h.subs(p, 1), " h'(0),h'(1):", sp.diff(h, p).subs(p, 0), sp.diff(h, p).subs(p, 1))


def check_s5():
    # 1D equilibrium profile and int phi'^2 = 1/(6 eps); energy = gamma; width
    p0 = (1 + sp.tanh(X/(2*e)))/2
    ode = sp.simplify(6*g/e*p0*(1-p0)*(1-2*p0) - 6*g*e*sp.diff(p0, X, 2))
    I = sp.integrate(sp.diff(p0, X)**2, (X, -sp.oo, sp.oo))
    E = sp.integrate(3*g/e*p0**2*(1-p0)**2 + 3*g*e*sp.diff(p0, X)**2, (X, -sp.oo, sp.oo))
    w = sp.solve(sp.Eq(p0, sp.Rational(9, 10)), X)[0] - sp.solve(sp.Eq(p0, sp.Rational(1, 10)), X)[0]
    print('S5 profile ODE residual:', ode, ' int phi0^2 =', sp.simplify(I), ' energy =', sp.simplify(E),
          ' 10-90 width =', sp.simplify(w))


def check_s6():
    # traveling wave solvability: v int phi'^2 = M dg int h'(phi) phi' = M dg
    M, dg = sp.symbols('M Delta_g', positive=True)
    p0 = (1 + sp.tanh(X/(2*e)))/2
    lhs_coef = sp.integrate(sp.diff(p0, X)**2, (X, -sp.oo, sp.oo))
    rhs = sp.integrate(30*p0**2*(1-p0)**2*sp.diff(p0, X), (X, -sp.oo, sp.oo))
    print('S6 int h\'(phi0) phi0\' dX =', sp.simplify(rhs), ' => v =', sp.simplify(M*dg*rhs/lhs_coef))


def check_units():
    import pint
    u = pint.UnitRegistry()
    gam = 1*u.J/u.m**2; ep = 1*u.m; M = 1*u.m**3/u.J/u.s
    rho_i = 1*u.kg/u.m**3; Rv = 1*u.J/u.kg/u.K; T = 1*u.K; rv = 1*u.kg/u.m**3
    terms = {'dphi/dt': 1/u.s, 'M gamma/eps': M*gam/ep, 'M gamma eps /m^2': M*gam*ep/u.m**2,
             'M rho_i R_v T (rho_v-rho_vs)/rho_vs': M*rho_i*Rv*T*rv/rv}
    for k, v in terms.items():
        print('units', k, v.to_base_units().units)


if __name__ == '__main__':
    check_s1(); check_s2(); check_s5(); check_s6(); check_units()
