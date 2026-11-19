import sympy as sp

x, y, z, t = sp.symbols('x y z t')
X = (x, y, z)
T = sp.Function('T')(x, y, z, t)
phi = sp.Function('phi')(x, y, z, t)
rho_i = sp.Symbol('rho_i', positive=True)
h_i = sp.Function('h_i'); h_v = sp.Function('h_v'); h_a = sp.Function('h_a')
k = sp.Function('k')
rho_v = sp.Function('rho_v')(x, y, z, t)
m_a = sp.Function('m_a')(x, y, z)  # air mass per unit volume, (1-phi) rho_a, time independent (A4)
jv = [sp.Function('j_%s' % c)(x, y, z, t) for c in 'xyz']


def E():
    return rho_i*phi*h_i(T) + m_a*h_a(T) + (1-phi)*rho_v*h_v(T)


def check_s4():
    # energy balance expanded, then water-mass balance substituted
    q = [-k(phi, T)*sp.diff(T, xi) + jv[n]*h_v(T) for n, xi in enumerate(X)]
    energy = sp.diff(E(), t) + sum(sp.diff(q[n], X[n]) for n in range(3))
    # water mass: d/dt(rho_i phi + (1-phi) rho_v) + div j_v = 0 -> solve for div j_v
    divj = -sp.diff(rho_i*phi + (1-phi)*rho_v, t)
    energy = energy.subs(sp.Derivative(jv[0], x), divj - sp.diff(jv[1], y) - sp.diff(jv[2], z))
    energy = sp.expand(energy.doit())
    c_i = sp.diff(h_i(T), t)/sp.diff(T, t)
    C = rho_i*phi*sp.diff(h_i(T), T) + m_a*sp.diff(h_a(T), T) + (1-phi)*rho_v*sp.diff(h_v(T), T)
    L = h_v(T) - h_i(T)
    target = (C*sp.diff(T, t) - sum(sp.diff(k(phi, T)*sp.diff(T, xi), xi) for xi in X)
              - rho_i*L*sp.diff(phi, t)
              + sum(jv[n]*sp.diff(h_v(T), X[n]) for n in range(3)))
    print('S4 residual (should be 0):', sp.simplify(sp.expand(energy - target.doit())))


def check_units():
    import pint
    u = pint.UnitRegistry()
    C = 1.4e6*u('J/m^3/K'); k = 2.0*u('W/m/K'); rho_i = 917*u('kg/m^3'); L = 2.8e6*u('J/kg')
    T = 200*u.K; s = 1*u.s; m = 1*u.m
    terms = [C*T/s, k*T/m**2, rho_i*L/s]
    for tm in terms:
        print(tm.to('W/m^3'))


def check_size():
    # vapor enthalpy-flux term j_v c_v grad T vs conduction k lap T; j_v ~ D grad rho_v
    D = 2e-5; rho_vs = 1e-7; dTdx = 1e3; k = 2.0; c_v = 1.4e3
    # rho_v gradient ~ rho_vs * (L/(R_v T^2)) dT/dx
    Lm = 2.8e6; Rv = 461.5; T = 200.
    jv = D*rho_vs*Lm/(Rv*T**2)*dTdx
    print('j_v c_v |gradT| =', jv*c_v*dTdx, 'W/m^3;  k |gradT|^2/T-scale', k*dTdx**2/1, 'per K  -> ratio', jv*c_v/k)


check_s4()
check_units()
check_size()
