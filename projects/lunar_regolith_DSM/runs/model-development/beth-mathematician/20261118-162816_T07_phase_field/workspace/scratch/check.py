import sympy as sp
import pint

x, y, z, t = sp.symbols('x y z t')
gamma, eps, M, alpha, rv, rvs = sp.symbols('gamma eps M alpha rho_v rho_vs', positive=True)
phi = sp.Function('phi')(x, y, z, t)
p = sp.Symbol('p')


def dF(expr_phi):
    """Euler-Lagrange derivative of the density with respect to phi."""
    f = 3*gamma/eps*expr_phi**2*(1-expr_phi)**2 + 3*gamma*eps*sum(sp.diff(expr_phi, v)**2 for v in (x, y, z))
    from sympy.calculus.euler import euler_equations
    eq = euler_equations(f, [expr_phi], [x, y, z])[0]
    return sp.simplify(eq.lhs)


def check_s1():
    fp = sp.diff(3*gamma/eps*p**2*(1-p)**2, p)
    print('S1 df/dphi - 6g/e p(1-p)(1-2p) =', sp.simplify(fp - 6*gamma/eps*p*(1-p)*(1-2*p)))


def check_s2():
    ph = sp.Function('phi')(x, y, z)  # t is a parameter in the variation
    el = dF(ph)
    target = 6*gamma/eps*ph*(1-ph)*(1-2*ph) - 6*gamma*eps*sum(sp.diff(ph, v, 2) for v in (x, y, z))
    print('S2 deltaF/deltaphi - target =', sp.simplify(el - target))


def check_s3():
    # flat profile phi0 = (1+tanh(x/(2 eps)))/2 is stationary: deltaF/deltaphi = 0
    X = sp.Symbol('X')
    p0 = (1 + sp.tanh(X/(2*eps)))/2
    r = 6*gamma/eps*p0*(1-p0)*(1-2*p0) - 6*gamma*eps*sp.diff(p0, X, 2)
    print('S3 deltaF/deltaphi on tanh profile =', sp.simplify(r.rewrite(sp.exp)))
    # source g = p^2(1-p)^2 vanishes with zero slope at 0 and 1
    g = p**2*(1-p)**2
    print('S3 g(0), g(1), g\'(0), g\'(1) =', [g.subs(p, 0), g.subs(p, 1), sp.diff(g, p).subs(p, 0), sp.diff(g, p).subs(p, 1)])


def check_s4():
    # dissipation: dF/dt = -M int (dF/dphi)^2 <= 0 when source off; sign of growth: rho_v>rho_vs -> phi_t>0
    print('S4 source sign at phi=1/2, rho_v-rho_vs=+1:', sp.sign(alpha*sp.Rational(1, 16)*1))


def check_units():
    u = pint.UnitRegistry()
    gam = 1*u.J/u.m**2; e = 1*u.m; Mq = 1*u.m**3/u.J/u.s; al = 1*u.m**3/u.kg/u.s; rho = 1*u.kg/u.m**3
    terms = {'phi_t': 1/u.s, 'M*g/e': Mq*gam/e, 'M*g*e/L^2': Mq*gam*e/u.m**2, 'alpha*drho': al*rho}
    for k, v in terms.items():
        print('units', k, v.to_base_units().units)


if __name__ == '__main__':
    check_s1(); check_s2(); check_s3(); check_s4(); check_units()
