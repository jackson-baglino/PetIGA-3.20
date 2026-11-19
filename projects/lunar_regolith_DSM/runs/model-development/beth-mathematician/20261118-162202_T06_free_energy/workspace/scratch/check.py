import sympy as sp

X = sp.symbols('X', real=True)
p = sp.symbols('p', real=True)
W, kap, g, e = sp.symbols('W kappa gamma epsilon', positive=True)

g_well = p**2*(1-p)**2


def check_s1():
    # double well: minima at 0 and 1, equal depth, barrier between
    d = sp.diff(g_well, p)
    print("roots of g':", sp.solve(d, p))
    print("g(0), g(1):", g_well.subs(p, 0), g_well.subs(p, 1))
    print("g''(0), g''(1):", sp.diff(g_well, p, 2).subs(p, 0), sp.diff(g_well, p, 2).subs(p, 1))
    print("g(1/2):", g_well.subs(p, sp.Rational(1, 2)))


def check_s3():
    # first integral: kappa/2 phi'^2 = W g(phi); verify profile phi=1/(1+exp(-X/l)), l=sqrt(kappa/(2W))
    l = sp.sqrt(kap/(2*W))
    ph = 1/(1+sp.exp(-X/l))
    lhs = kap/2*sp.diff(ph, X)**2
    rhs = W*g_well.subs(p, ph)
    print("first integral residual:", sp.simplify(lhs-rhs))
    # Euler-Lagrange: kappa phi'' = W g'(phi)
    el = kap*sp.diff(ph, X, 2) - W*sp.diff(g_well, p).subs(p, ph)
    print("EL residual:", sp.simplify(el))


def check_s4():
    # excess energy sigma = int (W g + kappa/2 phi'^2) dX = int sqrt(2 W kappa) p(1-p) dp
    sig = sp.integrate(sp.sqrt(2*W*kap)*p*(1-p), (p, 0, 1))
    print("sigma =", sp.simplify(sig))
    # direct integration in X with numbers
    Wn, kn = 3, 2
    l = sp.sqrt(sp.Rational(kn, 2*Wn))
    ph = 1/(1+sp.exp(-X/l))
    dens = Wn*g_well.subs(p, ph) + sp.Rational(kn, 2)*sp.diff(ph, X)**2
    val = sp.Integral(dens, (X, -sp.oo, sp.oo)).evalf()
    print("numeric sigma:", val, "formula:", sp.sqrt(2*Wn*kn).evalf()/6)


def check_s5():
    # choose l = eps/2 and sigma = gamma
    sol = sp.solve([sp.sqrt(2*W*kap)/6 - g, kap/(2*W) - e**2/4], [W, kap], dict=True)
    print("solution:", sol)


def check_s6():
    # profile with W=6g/e, kappa=3 g e: phi=(1+tanh(X/e))/2; width by max slope and 10-90
    Wv, kv = 6*g/e, 3*g*e
    ph = (1+sp.tanh(X/e))/2
    el = kv*sp.diff(ph, X, 2) - Wv*sp.diff(g_well, p).subs(p, ph)
    print("EL residual tanh profile:", sp.simplify(el.rewrite(sp.exp)))
    slope0 = sp.diff(ph, X).subs(X, 0)
    print("1/max slope =", sp.simplify(1/slope0))
    x90 = sp.solve(sp.Eq(ph, sp.Rational(9, 10)), X)
    print("10-90 width =", [sp.simplify(2*s) for s in x90], "=", [sp.N(2*s/e) for s in x90], "eps")
    dens = Wv*g_well.subs(p, ph) + kv/2*sp.diff(ph, X)**2
    val = sp.integrate(sp.simplify(dens.rewrite(sp.exp)), (X, -sp.oo, sp.oo))
    print("sigma direct =", sp.simplify(val))


def check_units():
    import pint
    u = pint.UnitRegistry()
    gam = 1*u.J/u.m**2
    eps = 1*u.m
    grad = 1/u.m
    t1 = 6*gam/eps
    t2 = 1.5*gam*eps*grad**2
    print("bulk term:", t1.to_base_units().units, " gradient term:", t2.to_base_units().units)
    print("same:", t1.dimensionality == t2.dimensionality == (u.J/u.m**3).dimensionality)


if __name__ == '__main__':
    for f in [check_s1, check_s3, check_s4, check_s5, check_s6, check_units]:
        print('---', f.__name__)
        try:
            f()
        except Exception as ex:
            print('ERROR', ex)
