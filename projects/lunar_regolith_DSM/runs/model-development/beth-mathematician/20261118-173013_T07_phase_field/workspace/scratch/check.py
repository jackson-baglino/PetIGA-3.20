import sympy as sp
from sympy.calculus.euler import euler_equations

x, y, z, t, n = sp.symbols('x y z t n', real=True)
gam, eps, M, alpha = sp.symbols('gamma epsilon M alpha', positive=True)
phi = sp.Function('phi')(x, y, z, t)
p = sp.Symbol('p')


def check_s1():
    L = 3*gam/eps*phi**2*(1-phi)**2 + 3*gam*eps*(phi.diff(x)**2 + phi.diff(y)**2 + phi.diff(z)**2)
    eq = euler_equations(L, phi, [x, y, z, t])[0].lhs  # dL/dphi - div(dL/dgrad)
    target = 6*gam/eps*phi*(1-phi)*(1-2*phi) - 6*gam*eps*(phi.diff(x, 2)+phi.diff(y, 2)+phi.diff(z, 2))
    print('S1 dF/dphi matches:', sp.simplify(eq - target) == 0)


def check_s3():
    q = sp.Symbol('q')
    pol = 30*q**2*(1-q)**2
    print('S3 p(0), p(1):', pol.subs(q, 0), pol.subs(q, 1), ' integral 0..1:', sp.integrate(pol, (q, 0, 1)))


def check_s5():
    phi0 = 1/(1+sp.exp(-n/eps))
    eqn = 6*gam/eps*phi0*(1-phi0)*(1-2*phi0) - 6*gam*eps*phi0.diff(n, 2)
    print('S5 profile solves dF/dphi=0:', sp.simplify(eqn) == 0)
    w = sp.solve(sp.Eq(phi0, sp.Rational(9, 10)), n)[0] - sp.solve(sp.Eq(phi0, sp.Rational(1, 10)), n)[0]
    print('S5 10-90 width:', sp.simplify(w))
    E = sp.integrate(sp.simplify(3*gam/eps*phi0**2*(1-phi0)**2 + 3*gam*eps*phi0.diff(n)**2), (n, -sp.oo, sp.oo))
    print('S5 energy per area:', sp.simplify(E))
    I2 = sp.integrate(sp.simplify(phi0.diff(n)**2), (n, -sp.oo, sp.oo))
    print('S5 int phi0_n^2 dn:', sp.simplify(I2))


def check_s6():
    # traveling flat front phi = phi0(n + V t)  (ice at n>0 grows toward -n at speed V)
    phi0 = lambda s: 1/(1+sp.exp(-s/eps))
    V, D = sp.symbols('V Delta', real=True)
    s = sp.Symbol('s', real=True)
    # solvability: V*int phi0'^2 = alpha*D*int p(phi0) phi0'
    lhs = V*sp.integrate(sp.simplify(phi0(s).diff(s)**2), (s, -sp.oo, sp.oo))
    rhs = alpha*D*sp.integrate(30*sp.Symbol('u')**2*(1-sp.Symbol('u'))**2, (sp.Symbol('u'), 0, 1))
    print('S6 V =', sp.solve(sp.Eq(lhs, rhs), V))
    # curvature: coefficient of laplacian -> V_curv = -6 gamma eps M kappa
    print('S6 diffusion coefficient on laplacian:', 6*gam*eps*M)


def check_units():
    import pint
    u = pint.UnitRegistry()
    gam_ = 1*u.J/u.m**2; eps_ = 1*u.m; L = 1*u.m; tt = 1*u.s; rho = 1*u.kg/u.m**3
    M_ = 1*u.m**3/(u.J*u.s); alpha_ = 1*u.m**3/(u.kg*u.s)
    terms = {'dphi/dt': 1/tt, 'M*gam/eps': M_*gam_/eps_, 'M*gam*eps/L^2': M_*gam_*eps_/L**2, 'alpha*rho': alpha_*rho}
    for k, v in terms.items():
        print('units', k, v.to_base_units().units)


for f in [check_s1, check_s3, check_s5, check_s6, check_units]:
    f()
