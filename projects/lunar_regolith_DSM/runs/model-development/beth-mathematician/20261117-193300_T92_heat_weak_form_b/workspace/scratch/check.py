import sympy as sp

x, t, L = sp.symbols('x t L', positive=True)
rho, c = sp.symbols('rho c', positive=True)
T = sp.Function('T')(x, t)
w = sp.Function('w')(x)
k = sp.Function('k')(x)
Q = sp.Function('Q')(x)

# Strong residual weighted by w
strong = rho*c*sp.diff(T, t) - sp.diff(k*sp.diff(T, x), x) - Q

# Integration by parts identity (pointwise): w*(k T_x)_x = (w k T_x)_x - w_x k T_x
lhs = w*sp.diff(k*sp.diff(T, x), x)
rhs = sp.diff(w*k*sp.diff(T, x), x) - sp.diff(w, x)*k*sp.diff(T, x)
print("IBP pointwise identity holds:", sp.simplify(lhs - rhs) == 0)

# Weak-form integrand + total derivative should equal w*strong
weak_integrand = rho*c*w*sp.diff(T, t) + k*sp.diff(w, x)*sp.diff(T, x) - w*Q
total_deriv = -sp.diff(w*k*sp.diff(T, x), x)
print("w*strong == weak integrand + boundary flux derivative:",
      sp.simplify(w*strong - (weak_integrand + total_deriv)) == 0)

# Concrete test with polynomials
Tt = x**2*t + x**3
wt = 1 + x
kt = 2 + x
Qt = x**2
s = (rho*c*sp.diff(Tt, t) - sp.diff(kt*sp.diff(Tt, x), x) - Qt)*wt
I_strong = sp.integrate(s, (x, 0, L))
flux = wt*kt*sp.diff(Tt, x)
R = (sp.integrate(rho*c*wt*sp.diff(Tt, t) + kt*sp.diff(wt, x)*sp.diff(Tt, x) - wt*Qt, (x, 0, L))
     - flux.subs(x, L) + flux.subs(x, 0))
print("Concrete check R == integral of w*strong:", sp.simplify(R - I_strong) == 0)
