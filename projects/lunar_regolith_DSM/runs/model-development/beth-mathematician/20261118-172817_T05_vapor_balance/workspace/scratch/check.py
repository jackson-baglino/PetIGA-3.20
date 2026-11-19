import sympy as sp
import pint

x, y, z, t = sp.symbols('x y z t')
phi = sp.Function('phi')(x, y, z, t)
T = sp.Function('T')(x, y, z, t)
rv = sp.Function('rho_v')(x, y, z, t)
D = sp.Function('D_v')(T)
ri = sp.Function('rho_i')(T)
X = (x, y, z)


def div(F):
    return sum(sp.diff(F[i], X[i]) for i in range(3))


def grad(f):
    return [sp.diff(f, c) for c in X]


def R_final():
    J = [-(1 - phi) * D * g for g in grad(rv)]
    return sp.diff((1 - phi) * rv, t) + sp.diff(ri * phi, t) + div(J)


def check_s2():
    # water mass density per unit total volume
    m = ri * phi + (1 - phi) * rv
    print('S2 m at phi=1:', sp.simplify(m.subs(phi, 1)), ' at phi=0:', sp.simplify(m.subs(phi, 0)))


def check_s4():
    # local form = d/dt m + div J, with J Fick flux in gas fraction
    m = ri * phi + (1 - phi) * rv
    J = [-(1 - phi) * D * g for g in grad(rv)]
    diff = sp.simplify(sp.diff(m, t) + div(J) - R_final())
    print('S4 local form minus R:', diff)


def check_s5():
    # expanded form with constant rho_i: (1-phi) d rho_v/dt + (rho_i - rho_v) dphi/dt - div(...)
    ric = sp.Symbol('rho_i')
    Rc = R_final().subs(ri, ric).doit()
    alt = (1 - phi) * sp.diff(rv, t) + (ric - rv) * sp.diff(phi, t) \
        - div([(1 - phi) * D * g for g in grad(rv)])
    print('S5 expanded minus R (rho_i const):', sp.simplify(Rc - alt))


def check_units():
    u = pint.UnitRegistry()
    rho = 1 * u.kg / u.m**3
    Dv = 1 * u.m**2 / u.s
    L = 1 * u.m
    s = 1 * u.s
    terms = {'d_t((1-phi) rho_v)': rho / s, 'd_t(rho_i phi)': rho / s,
             'div((1-phi) D grad rho_v)': Dv * rho / L**2}
    for k, v in terms.items():
        print('units', k, v.to_base_units().units)


if __name__ == '__main__':
    check_s2(); check_s4(); check_s5(); check_units()
    print('R =', R_final())
