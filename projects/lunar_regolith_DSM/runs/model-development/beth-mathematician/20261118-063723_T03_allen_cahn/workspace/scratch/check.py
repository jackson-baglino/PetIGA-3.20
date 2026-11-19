from sympy import *
from sympy.calculus.euler import euler_equations
x, y, t = symbols('x y t')
kappa, W, M = symbols('kappa W M', positive=True)
phi = Function('phi')
p = phi(x, y, t)
L = kappa/2*(Derivative(p, x)**2 + Derivative(p, y)**2) + W*p**2*(1-p)**2
# Euler-Lagrange: dL/dphi - d/dx dL/dphi_x - d/dy dL/dphi_y = 0  -> lhs is mu
eq = euler_equations(L, [p], [x, y, t])[0]
mu_EL = eq.lhs
mu = 2*W*p*(1-p)*(1-2*p) - kappa*(Derivative(p, (x, 2)) + Derivative(p, (y, 2)))
print('mu check:', simplify(mu_EL.doit() - mu.doit()))
R = Derivative(p, t) + M*(2*W*p*(1-p)*(1-2*p) - kappa*(Derivative(p, (x, 2)) + Derivative(p, (y, 2))))
print('R check:', simplify(R - (Derivative(p, t) + M*mu)))
print('df check:', simplify(diff(W*symbols('u')**2*(1-symbols('u'))**2, symbols('u')) - 2*W*symbols('u')*(1-symbols('u'))*(1-2*symbols('u'))))
