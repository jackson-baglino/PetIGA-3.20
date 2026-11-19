import sympy as sp
x, y, z, t = sp.symbols('x y z t')
X = (x, y, z)
T = sp.Function('T')(x, y, z, t)
phi = sp.Function('phi')(x, y, z, t)
rho_i = sp.Symbol('rho_i', positive=True)
rho_v = sp.Function('rho_v')(x, y, z, t)
rho_a = sp.Function('rho_a')(x, y, z, t)
h_i, h_v, h_a = [sp.Function(n) for n in ('h_i', 'h_v', 'h_a')]
k = sp.Function('k')(phi, T)
j = [sp.Function('j' + c)(x, y, z, t) for c in 'xyz']
div = lambda F: sum(sp.diff(F[n], X[n]) for n in range(3))
grad = lambda f: [sp.diff(f, v) for v in X]

# S1: enthalpy density, S2: energy flux
H = phi*rho_i*h_i(T) + (1-phi)*rho_v*h_v(T) + (1-phi)*rho_a*h_a(T)
q = [-k*g + h_v(T)*jj for g, jj in zip(grad(T), j)]
energy = sp.diff(H, t) + div(q)                       # S3: = 0
water = sp.diff(phi*rho_i + (1-phi)*rho_v, t) + div(j)  # S4: = 0
air = sp.diff((1-phi)*rho_a, t)                         # S5: = 0
# S6: subtract h_v*water + h_a*air (both zero) from energy
S6 = sp.expand(energy - h_v(T)*water - h_a(T)*air)
C = phi*rho_i*sp.diff(h_i(T), T) + (1-phi)*rho_v*sp.diff(h_v(T), T) + (1-phi)*rho_a*sp.diff(h_a(T), T)
L = h_v(T) - h_i(T)
claimed = (C*sp.diff(T, t) - div([k*g for g in grad(T)]) - rho_i*L*sp.diff(phi, t)
           + sum(jj*sp.diff(h_v(T), v) for jj, v in zip(j, X)))
print('S6 identity residual:', sp.simplify(sp.expand(S6 - claimed)))

# units check with pint
import pint
u = pint.UnitRegistry()
C_ = 2e6*u.J/u.m**3/u.K; Tt = 1*u.K/u.s; k_ = 2*u.W/u.m/u.K; Txx = 1*u.K/u.m**2
rho_i_ = 917*u.kg/u.m**3; L_ = 2.8e6*u.J/u.kg; phit = 1/u.s
terms = [C_*Tt, k_*Txx, rho_i_*L_*phit]
print([str(tt.to_base_units().units) for tt in terms])
# size of dropped vapor-enthalpy transport term relative to conduction (per unit grad T)
c_v = 1850*u.J/u.kg/u.K; D = 1e-3*u.m**2/u.s; rho_v_ = 1.7e-6*u.kg/u.m**3
print('c_v D rho_v / k =', (c_v*D*rho_v_/k_).to_base_units())
print('rho_v c_v / C =', (rho_v_*c_v/C_).to_base_units())
