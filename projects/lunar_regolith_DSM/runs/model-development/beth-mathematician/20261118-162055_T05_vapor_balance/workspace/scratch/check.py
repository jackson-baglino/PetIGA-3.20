import sympy as sp

x, y, z, t = sp.symbols('x y z t')
phi = sp.Function('phi')(x, y, z, t)
T = sp.Function('T')(x, y, z, t)
rho_v = sp.Function('rho_v')(x, y, z, t)
rho_i = sp.Function('rho_i')
D = sp.Function('D')
X = (x, y, z)


def check_s1():
    # water mass per unit total volume: ice part + vapor part
    m = phi * rho_i(T) + (1 - phi) * rho_v
    # limits: pure ice gives rho_i, pure gas gives rho_v
    print('S1 phi=1:', sp.simplify(m.subs(phi, 1)), ' phi=0:', sp.simplify(m.subs(phi, 0)))


def check_s4():
    # divergence theorem + arbitrary volume -> local form; check that the local form
    # integrated over a box with zero flux at the walls gives d/dt of total mass = 0 (1D test)
    s = sp.symbols('s')
    f = sp.Function('f')
    J = -D(1, 1) * sp.diff(f(s), s)
    # int_0^L dJ/ds ds = J(L)-J(0)
    L = sp.symbols('L', positive=True)
    lhs = sp.integrate(sp.diff(J, s), (s, 0, L))
    print('S4 int of div J =', sp.simplify(lhs - (J.subs(s, L) - J.subs(s, 0))))


def R_expr():
    m = phi * rho_i(T) + (1 - phi) * rho_v
    flux = [-D(phi, T) * sp.diff(rho_v, xi) for xi in X]
    return sp.diff(m, t) + sum(sp.diff(flux[k], X[k]) for k in range(3))


def check_s5():
    R = R_expr()
    # expanded form with rho_i constant must reduce to the familiar one
    ri = sp.symbols('rho_i0', positive=True)
    Rc = R.replace(rho_i, lambda a: ri).doit()
    familiar = ((1 - phi) * sp.diff(rho_v, t) + (ri - rho_v) * sp.diff(phi, t)
                - sum(sp.diff(D(phi, T) * sp.diff(rho_v, xi), xi) for xi in X))
    print('S5 R(const rho_i) - familiar =', sp.simplify(Rc - familiar))
    # no exchange when phi is constant and T constant: pure diffusion
    print('S5 R expanded:', R)


def check_units():
    import pint
    u = pint.UnitRegistry()
    rho = 1 * u.kg / u.m**3
    Dq = 1 * u.m**2 / u.s
    terms = {
        'd_t((1-phi) rho_v)': rho / u.s,
        'd_t(phi rho_i)': rho / u.s,
        'd_x(D d_x rho_v)': Dq * rho / u.m**2,
        'sigma * d_t phi': rho / u.s,
    }
    for k, v in terms.items():
        print('units', k, v.to_base_units().units)


if __name__ == '__main__':
    check_s1(); check_s4(); check_s5(); check_units()
