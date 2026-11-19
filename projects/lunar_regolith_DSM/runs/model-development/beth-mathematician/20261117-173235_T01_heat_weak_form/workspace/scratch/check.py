import sympy as sp

x, t, L, rho, c = sp.symbols('x t L rho c', positive=True)
T = sp.Function('T')(x, t)
w = sp.Function('w')(x)
k = sp.Function('k')(x)
Q = sp.Function('Q')(x)

# Strong residual times w, integrated
strong = w*(rho*c*sp.diff(T, t) - sp.diff(k*sp.diff(T, x), x) - Q)

# Weak residual after integration by parts
R = (sp.Integral(rho*c*w*sp.diff(T, t), (x, 0, L))
     + sp.Integral(k*sp.diff(w, x)*sp.diff(T, x), (x, 0, L))
     - sp.Integral(w*Q, (x, 0, L))
     - sp.Subs(w*k*sp.diff(T, x), x, L)
     + sp.Subs(w*k*sp.diff(T, x), x, 0))

# Check pointwise identity: integrand difference is an exact x-derivative of w k T_x
integrand_weak = rho*c*w*sp.diff(T, t) + k*sp.diff(w, x)*sp.diff(T, x) - w*Q
flux = w*k*sp.diff(T, x)
print("product rule check:", sp.simplify(strong - (integrand_weak - sp.diff(flux, x))))

# Concrete check with polynomial choices
Tc = x**2*t + sp.sin(x)*t**2
wc = 1 + x + x**3
kc = 2 + x**2
Qc = sp.cos(x)
sub = {T: Tc, w: wc, k: kc, Q: Qc}
lhs = sp.integrate((wc*(rho*c*sp.diff(Tc, t) - sp.diff(kc*sp.diff(Tc, x), x) - Qc)), (x, 0, L))
rhs = (sp.integrate(rho*c*wc*sp.diff(Tc, t) + kc*sp.diff(wc, x)*sp.diff(Tc, x) - wc*Qc, (x, 0, L))
       - (wc*kc*sp.diff(Tc, x)).subs(x, L) + (wc*kc*sp.diff(Tc, x)).subs(x, 0))
print("concrete check:", sp.simplify(lhs - rhs))
print(R)
