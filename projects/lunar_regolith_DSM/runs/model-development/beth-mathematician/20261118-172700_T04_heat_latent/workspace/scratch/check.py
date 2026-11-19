import sympy as sp

x, y, z, t = sp.symbols('x y z t')
X = (x, y, z, t)
T = sp.Function('T')(*X)
phi = sp.Function('phi')(*X)
rho_i, rho_a = sp.symbols('rho_i rho_a', positive=True)
rho_v = sp.Function('rho_v')(*X)
h_i, h_a, h_v = [sp.Function(n) for n in ('h_i', 'h_a', 'h_v')]
k = sp.Function('k')
jv = [sp.Function('jv%s' % c)(*X) for c in 'xyz']
ja = [sp.Function('ja%s' % c)(*X) for c in 'xyz']
crd = (x, y, z)
div = lambda F: sum(sp.diff(F[i], crd[i]) for i in range(3))
grad = lambda f: [sp.diff(f, c) for c in crd]

H = phi*rho_i*h_i(T) + (1-phi)*(rho_a*h_a(T) + rho_v*h_v(T))
q = [-k(phi, T)*g for g in grad(T)]


def energy_residual():
    return sp.diff(H, t) + div([q[i] + h_v(T)*jv[i] + h_a(T)*ja[i] for i in range(3)])


def check_s4():
    # use the mass balances (S2, S3) to eliminate div j_v and div j_a
    E = energy_residual()
    divjv = -rho_i*sp.diff(phi, t) - sp.diff((1-phi)*rho_v, t)
    divja = -sp.diff((1-phi)*rho_a, t)
    E = E.subs({sp.Derivative(jv[0], x): divjv - sp.Derivative(jv[1], y) - sp.Derivative(jv[2], z),
                sp.Derivative(ja[0], x): divja - sp.Derivative(ja[1], y) - sp.Derivative(ja[2], z)})
    ci, ca, cv = [sp.diff(f(T), T) for f in (h_i, h_a, h_v)]
    C = phi*rho_i*ci + (1-phi)*(rho_a*ca + rho_v*cv)
    L = h_v(T) - h_i(T)
    adv = sum((jv[i]*cv + ja[i]*ca)*grad(T)[i] for i in range(3))
    expected = C*sp.diff(T, t) + adv - div([k(phi, T)*g for g in grad(T)]) - rho_i*L*sp.diff(phi, t)
    d = sp.simplify(sp.expand((E - expected).doit()))
    print('S4 energy - expected =', d)


def check_units():
    import pint
    u = pint.UnitRegistry()
    C = 1.0*u.J/u.m**3/u.K; k_ = 1.0*u.W/u.m/u.K; Tq = 1.0*u.K; tq = 1.0*u.s; L = 1.0*u.J/u.kg
    r = 1.0*u.kg/u.m**3; m = 1.0*u.m; j = 1.0*u.kg/u.m**2/u.s; c = 1.0*u.J/u.kg/u.K
    terms = {'C dT/dt': C*Tq/tq, 'div(k grad T)': k_*Tq/m**2, 'rho_i L dphi/dt': r*L/tq,
             'j c grad T': j*c*Tq/m}
    for n, v in terms.items():
        print(n, v.to(u.W/u.m**3))


if __name__ == '__main__':
    check_s4()
    check_units()
