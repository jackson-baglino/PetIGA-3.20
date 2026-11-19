import sympy as sp

p, u, g, e, A, K = sp.symbols('phi u gamma epsilon A K', positive=True)
W = p**2*(1-p)**2

def check_s1():
    # wells: W(0)=W(1)=0, W''>0 there; barrier max at 1/2
    d = sp.diff(W, p)
    print("S1 roots W':", sp.solve(d, p), " W''(0),W''(1):", sp.diff(W, p, 2).subs(p, 0), sp.diff(W, p, 2).subs(p, 1), " W(1/2):", W.subs(p, sp.Rational(1, 2)))

def check_s3():
    # first integral: d/dx[K/2 phi'^2 - A W] = phi'(K phi'' - A W') = 0 on EL solutions
    x = sp.symbols('x')
    f = sp.Function('f')(x)
    lag = K/2*f.diff(x)**2 + A*W.subs(p, f)
    el = sp.euler_equations(lag, f, x)[0].lhs
    H = K/2*f.diff(x)**2 - A*W.subs(p, f)
    print("S3 dH/dx + phi'*EL :", sp.simplify(H.diff(x) + f.diff(x)*el))

def check_s4():
    x = sp.symbols('x', real=True)
    k = sp.sqrt(2*A/K)
    prof = 1/(1+sp.exp(-k*x))
    res = sp.simplify(K*prof.diff(x, 2) - A*sp.diff(W, p).subs(p, prof))
    print("S4 EL residual of logistic profile:", res)

def check_s5():
    gam = sp.integrate(sp.sqrt(2*A*K)*p*(1-p), (p, 0, 1))
    print("S5 gamma =", sp.simplify(gam))

def check_s6():
    sol = sp.solve([sp.sqrt(2*A*K)/6 - g, sp.sqrt(2*A/K) - 1/e], [A, K], dict=True)
    print("S6 A,K =", sol)

def check_s7():
    x = sp.symbols('x', real=True)
    prof = 1/(1+sp.exp(-x/e))
    print("S7 profile == (1+tanh(x/2eps))/2:", sp.simplify((prof - (1+sp.tanh(x/(2*e)))/2).rewrite(sp.exp)))
    slope0 = prof.diff(x).subs(x, 0)
    print("S7 width 1/phi'(0) =", sp.simplify(1/slope0))
    x9 = sp.solve(sp.Eq(prof, sp.Rational(9, 10)), x)
    print("S7 10-90 width =", sp.simplify(2*x9[0]), "=", sp.N(2*x9[0]/e), "eps")
    # energy check by integrating f over x with final coefficients
    f = 3*g/e*W.subs(p, prof) + 3*g*e*prof.diff(x)**2
    print("S7 int f dx =", sp.simplify(sp.integrate(sp.simplify(f.rewrite(sp.exp)), (x, -sp.oo, sp.oo))))

def check_units():
    import pint
    ur = pint.UnitRegistry()
    gq = 1*ur.J/ur.m**2; eq = 1*ur.m
    t1 = (3*gq/eq)
    t2 = 3*gq*eq*(1/ur.m)**2
    print("units:", t1.to_base_units().units, "|", t2.to_base_units().units, "| J/m^3 =", (1*ur.J/ur.m**3).to_base_units().units)

for fn in [check_s1, check_s3, check_s4, check_s5, check_s6, check_s7, check_units]:
    fn()
