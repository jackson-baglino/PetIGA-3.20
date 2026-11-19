"""M3 check: conservation and energy decay of the strong form (S15/S17 from the ledger).

The toy has no ice or vapor (S8), so there is no total water to conserve; the
conserved quantity is internal energy int e dx with e from S3.
  * closed box (insulated ends, j = 0): d/dt int e = int q  (= 0 for q = 0);
  * Dirichlet box: d/dt int e = int q - (F_0 + F_L), with F from S39;
  * the same balance on the ledger's series solution (S54) in theta form;
  * without the source, E = 1/2 int (T - T0)^2 does not increase: symbolic identity
    dE/dt = -kappa int T_x^2 on test fields, and a backward-Euler solve of S17;
  * with the source, the discrete steady state of S17 equals S44.
"""
import os
import random
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S15', 'S17', 'S39', 'S44', 'S54']

x, t, xi, tau = (C.S(n) for n in ('x', 't', 'xi', 'tau'))
Lc, k, rho, cp, kappa, q, T0 = (C.S(n) for n in ('L', 'k', 'rho', 'c_p', 'kappa', 'q', 'T_0'))
F0, FL = C.S('F_0'), C.S('F_L')
PARAMS = {Lc: 1.0, k: 2.0, rho: 1.5, cp: 0.8, T0: 1.0}


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_conservation_energy.py', L,
                     {'numeric_params': {str(a): b for a, b in PARAMS.items()}})
    Tf = C.field('T')
    Tt, Txx = sp.Derivative(Tf, t), sp.Derivative(Tf, (x, 2))

    def rate(sid, Texpr):
        r = C.piece(L.stmt(sid), 0).residual()
        sol = sp.solve(r, Tt)
        if len(sol) != 1:
            raise ValueError(f'{sid} not solvable for dT/dt')
        return C.plug(sol[0], 'T', sp.Lambda((x, t), Texpr))

    def e_law():
        return C.find_piece(L.stmt('S3'), C.lhs_is('e'), 'e = ...').rhs

    def d_dt_energy(Texpr):
        # d/dt int e dx = int (d e / d T) T_t dx, e from S3, T_t from S15
        e = e_law()
        dedT = sp.diff(e.subs(Tf, sp.Symbol('Tv')), sp.Symbol('Tv'))
        return sp.integrate(dedT * rate('S15', Texpr), (x, 0, Lc))

    def closed_box():
        Tc = T0 + x**2 * (Lc - x)**2 * (1 + t) + sp.cos(sp.pi * x / Lc) * sp.exp(-t)
        ends = [sp.simplify(sp.diff(Tc, x).subs(x, 0)), sp.simplify(sp.diff(Tc, x).subs(x, Lc))]
        if any(e_ != 0 for e_ in ends):
            return False, 'test field is not insulated'
        dEdt = sp.simplify(d_dt_energy(Tc))
        src = sp.integrate(q, (x, 0, Lc))
        ok = sp.simplify(dEdt - src) == 0 and sp.simplify(dEdt.subs(q, 0)) == 0
        return ok, {'d/dt int e (S3, S15)': str(dEdt), 'int q': str(src),
                    'conserved for q = 0': sp.simplify(dEdt.subs(q, 0)) == 0}
    rep.check('closed_box_energy', 'closed (insulated) box: total energy changes only by int q; conserved for q = 0',
              closed_box)

    def dirichlet_balance():
        Tc = T0 + sp.sin(2 * x + t) * x + x**3 * sp.exp(-t)
        dEdt = d_dt_energy(Tc)
        jlaw = C.find_piece(L.stmt('S4'), C.lhs_is('j'), 'j = ...').rhs
        jl = sp.Lambda((x, t), C.plug(jlaw, 'T', sp.Lambda((x, t), Tc)))
        F = {}
        for p in C.pieces(L.stmt('S39')):
            F[p.exprs[0]] = C.plug(C.plug(p.exprs[-1], 'T', sp.Lambda((x, t), Tc)), 'j', jl)
        bal = sp.simplify(dEdt - (sp.integrate(q, (x, 0, Lc)) - F[F0] - F[FL]))
        return bal == 0, {'d/dt int e - (int q - F0 - FL)': str(bal)}
    rep.check('dirichlet_energy_balance', 'Dirichlet box: d/dt int e = int q - F_0 - F_L (S3, S15, S39)',
              dirichlet_balance)

    def series_balance():
        # theta form of the same balance: d/dtau int theta = 1 + theta_xi(1) - theta_xi(0), on the S54 series
        ps = [p for p in C.pieces(L.stmt('S54')) if p.ok]
        pv = [p for p in ps if p.exprs[0] == C.field('v')][0]
        th43 = [p for p in C.pieces(L.stmt('S43')) if p.ok and p.exprs[0] == C.field('theta_s')][0].rhs
        sm = [a for a in sp.preorder_traversal(pv.rhs) if isinstance(a, sp.Sum)][0]
        r = C.S('r')
        f = sm.function
        I_term = sp.integrate(f, (xi, 0, 1))
        rows = []
        for tv in (0.02, 0.1, 0.5):
            nn = np.arange(4000, dtype=float)
            g = sp.lambdify((r, tau), I_term, 'numpy')
            gt = sp.lambdify((r, tau), sp.diff(I_term, tau), 'numpy')
            dx0 = sp.lambdify((r, tau), sp.diff(f, xi).subs(xi, 0), 'numpy')
            dx1 = sp.lambdify((r, tau), sp.diff(f, xi).subs(xi, 1), 'numpy')
            d_int_theta = -float(np.sum(gt(nn, tv)))           # theta = theta_s - v
            th_x0 = float(sp.diff(th43, xi).subs(xi, 0)) - float(np.sum(dx0(nn, tv)))
            th_x1 = float(sp.diff(th43, xi).subs(xi, 1)) - float(np.sum(dx1(nn, tv)))
            c25 = -C.piece(L.stmt('S25'), 0).residual().subs(
                {sp.Derivative(C.field('theta'), tau): 0, sp.Derivative(C.field('theta'), (xi, 2)): 0})
            rhs = float(c25) + th_x1 - th_x0
            rows.append((tv, d_int_theta, rhs))
        err = max(abs(a_ - b_) for _, a_, b_ in rows)
        return err < 1e-6, {'(tau, d/dtau int theta, source + theta_xi(1) - theta_xi(0))': rows, 'max diff': err}
    rep.check('series_energy_balance', 'the S54 series obeys the global balance d/dtau int theta = 1 - outward fluxes',
              series_balance)

    def lyapunov_identity():
        # E = 1/2 int (T - T0)^2, q = 0, T = T0 at both ends: dE/dt = int (T - T0) T_t = -kappa int T_x^2
        kdef = C.piece(L.stmt('S16'), 0).exprs[1]
        out, ok = [], True
        for Tc in (T0 + x * (Lc - x) * sp.exp(-t) * (1 + x), T0 + sp.sin(sp.pi * x / Lc) * sp.cos(t)):
            Tt_ = rate('S17', Tc).subs(q, 0).subs(kappa, kdef)
            dE = sp.integrate((Tc - T0) * Tt_, (x, 0, Lc))
            target = -kdef * sp.integrate(sp.diff(Tc, x)**2, (x, 0, Lc))
            diff_ = sp.simplify(dE - target)   # exact, symbolic L, k, rho, c_p
            val = float(dE.subs(PARAMS).subs(t, 0.3))
            ok &= diff_ == 0 and val <= 0
            out.append((str(diff_), val))
        return ok, {'(dE/dt + kappa int T_x^2, dE/dt at t=0.3) on test fields': out}
    rep.check('energy_identity_no_source',
              'without the source, dE/dt = -kappa int T_x^2 <= 0 for E = 1/2 int (T-T0)^2 (S17, Dirichlet)',
              lyapunov_identity)

    def coeffs():
        r17 = C.piece(L.stmt('S17'), 0).residual()
        kdef = C.piece(L.stmt('S16'), 0).exprs[1]
        rhs = sp.solve(r17, Tt)[0].subs(kappa, kdef)
        X = sp.Symbol('Txx_')
        rhs = sp.expand(rhs.subs(Txx, X))
        if rhs.has(Tf):
            raise ValueError(f'S17 right side is not linear in d2T/dx2 only: {rhs}')
        return rhs.coeff(X), rhs.subs(X, 0)

    def fd_matrix(n, h):
        A = np.diag(-2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
        return A / h**2

    def energy_decay():
        alpha, beta = coeffs()
        a_num = float(alpha.subs(PARAMS))
        b_num = float(beta.subs(PARAMS).subs(q, 0))
        n, Lv = 199, PARAMS[Lc]
        h = Lv / (n + 1)
        xs = np.linspace(h, Lv - h, n)
        rng = random.Random(11)
        u = sum(rng.uniform(-1, 1) * np.sin((m + 1) * np.pi * xs / Lv) for m in range(8))
        A = fd_matrix(n, h)
        dt = 1e-3
        M = np.eye(n) - dt * a_num * A
        E = [0.5 * h * float(u @ u)]
        for _ in range(300):
            u = np.linalg.solve(M, u + dt * b_num)
            E.append(0.5 * h * float(u @ u))
        incr = max(np.diff(E))
        ok = incr <= 1e-14 * E[0] and E[-1] < E[0]
        return ok, {'coef of d2T/dx2 (S17)': str(alpha), 'source at q=0': b_num,
                    'E0': E[0], 'E_end': E[-1], 'max step increase': incr}
    rep.check('energy_decay_no_source', 'without the source, E = 1/2 int (T-T0)^2 does not increase '
              '(backward-Euler solve of S17, Dirichlet)', energy_decay)

    def steady_vs_S44():
        alpha, beta = coeffs()
        qv = 3.0
        a_num = float(alpha.subs(PARAMS))
        b_num = float(beta.subs(PARAMS).subs(q, qv))
        n, Lv = 99, PARAMS[Lc]
        h = Lv / (n + 1)
        xs = np.linspace(h, Lv - h, n)
        A = fd_matrix(n, h)
        bc = np.zeros(n)
        bc[0] = bc[-1] = PARAMS[T0] / h**2
        u = np.linalg.solve(a_num * A, -b_num - a_num * bc)
        Ts = C.piece(L.stmt('S44'), 0).rhs.subs(PARAMS).subs(q, qv)
        ref = np.array([float(Ts.subs(x, xv)) for xv in xs])
        rel = float(np.max(np.abs(u - ref)) / np.max(np.abs(ref - PARAMS[T0])))
        return rel <= 0.01, {'max |FD - S44| / max |S44 - T0|': rel, 'q': qv}
    rep.check('steady_state_vs_S44', 'the steady state of S17 (finite differences) equals S44 within 1 %',
              steady_vs_S44)

    rep.finish()


if __name__ == '__main__':
    main()
