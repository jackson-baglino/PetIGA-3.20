import sympy as sp

x, y, z, t = sp.symbols('x y z t')
phi = sp.Function('phi')(x, y, z, t)
T = sp.Function('T')(x, y, z, t)
rho_i = sp.Symbol('rho_i', positive=True)
hi = sp.Function('h_i')
hv = sp.Function('h_v')
rv = sp.Function('rho_v')(x, y, z, t)
jx, jy, jz = [sp.Function(n)(x, y, z, t) for n in ('j_x', 'j_y', 'j_z')]
mdot = sp.Function('mdot')(x, y, z, t)
X = (x, y, z)
J = (jx, jy, jz)


def div(F):
    return sum(sp.diff(F[i], X[i]) for i in range(3))


def check_s4():
    # energy density e and vapor enthalpy flux; rate of e + div(h_v j_v) minus sensible part
    e = phi*rho_i*hi(T) + (1 - phi)*rv*hv(T)
    lhs = sp.diff(e, t) + div([hv(T)*J[i] for i in range(3)])
    # use mass balances: d(phi rho_i)/dt = -mdot, d((1-phi)rho_v)/dt = mdot - div j
    sens = phi*rho_i*sp.diff(hi(T), t) + (1 - phi)*rv*sp.diff(hv(T), t)
    rest = (hi(T)*(-mdot) + hv(T)*(mdot - div(J)) + div([hv(T)*J[i] for i in range(3)]))
    # lhs with substitution: replace time derivative of masses
    lhs_sub = (sens + hi(T)*sp.diff(phi*rho_i, t) + hv(T)*sp.diff((1 - phi)*rv, t)
               + div([hv(T)*J[i] for i in range(3)]))
    print('product rule ok:', sp.simplify(sp.expand(lhs - lhs_sub)) == 0)
    target = sens + mdot*(hv(T) - hi(T)) + sum(J[i]*sp.diff(hv(T), X[i]) for i in range(3))
    print('mass-balance substitution ok:', sp.simplify(sp.expand(rest + sens - target)) == 0)


def check_s7():
    import pint
    u = pint.UnitRegistry()
    C = 1.4e6*u('J/m^3/K'); k = 3.5*u('W/m/K'); rho = 917*u('kg/m^3'); L = 2.8e6*u('J/kg')
    Tt = 1*u('K/s'); phit = 1*u('1/s'); lap = 1*u('K/m^2')
    for name, q in [('C T_t', C*Tt), ('div k grad T', k*lap), ('rho L phi_t', rho*L*phit)]:
        print(name, q.to('W/m^3'))


if __name__ == '__main__':
    check_s4()
    check_s7()
