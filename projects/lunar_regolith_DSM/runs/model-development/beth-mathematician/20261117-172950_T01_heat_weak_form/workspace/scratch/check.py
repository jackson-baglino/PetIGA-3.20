import sympy as sp

x, t, L, rho, c = sp.symbols('x t L rho c', positive=True)
T = sp.Function('T')(x, t)
w = sp.Function('w')(x)
k = sp.Function('k')(x)
Q = sp.Function('Q')(x)
flux = k*sp.Derivative(T, x)

R = (sp.Integral(rho*c*sp.Derivative(T, t)*w, (x, 0, L))
     + sp.Integral(k*sp.Derivative(T, x)*sp.Derivative(w, x), (x, 0, L))
     - sp.Integral(Q*w, (x, 0, L))
     - sp.Subs(flux*w, x, L)
     + sp.Subs(flux*w, x, 0))

# integral of (strong residual) * w; must equal R for any smooth fields
strong = sp.Integral((rho*c*sp.Derivative(T, t) - sp.Derivative(flux, x) - Q)*w, (x, 0, L))

Tf = sp.sin(2*x)*sp.exp(-t) + x**3
wf = sp.cos(x) + x**2
kf = 1 + x**2
Qf = sp.exp(x)


def concrete(e):
    return e.xreplace({T: Tf, w: wf, k: kf, Q: Qf}).doit()


diff = sp.simplify(concrete(R) - concrete(strong))
print("symbolic difference:", diff)
vals = {L: sp.Rational(13, 10), rho: 2, c: 3, t: sp.Rational(1, 2)}
print("numeric difference:", sp.N((concrete(R) - concrete(strong)).subs(vals)))

# second set of fields
Tf = x**2*t + sp.cos(x)
wf = sp.exp(-x)
kf = 2 + sp.sin(x)
Qf = x
print("second check:", sp.simplify(concrete(R) - concrete(strong)))
