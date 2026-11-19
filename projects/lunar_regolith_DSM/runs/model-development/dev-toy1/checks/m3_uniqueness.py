"""M3 check: uniqueness, the q = 0 case, the minimum principle and the sign of T - T0
(S33-S37, S48-S50, S38).

S33: linearity on two arbitrary test fields T1, T2 pushed through S17.
S34: dE/dt from the ledger's own E definition, with w_t replaced from the ledger's S33.
S35: integration by parts on test fields that vanish at the ends (the S33 data).
S36-S37: each premise the logic rests on is recomputed.
S48-S50: the inequalities claimed at an interior/top minimum are tested on functions
that have such a minimum; the conclusion theta >= 0 is tested against an
independent finite-difference solve of the ledger's own S25-S26.
S38: parameter-free theta equation, sign of DeltaT_q, and theta >= 0 numerically.
"""
import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S33', 'S34', 'S35', 'S36', 'S37', 'S48', 'S49', 'S50', 'S38']

x, t, xi, tau = (C.S(n) for n in ('x', 't', 'xi', 'tau'))
Lc, k, rho, cp, kappa, q, T0 = (C.S(n) for n in ('L', 'k', 'rho', 'c_p', 'kappa', 'q', 'T_0'))
LN, KN = sp.Rational(2), sp.Rational(3, 10)  # numeric L, kappa where integrals need numbers
T1_TEST = sp.sin(x + t) * x**2 + sp.exp(-t) * x
T2_TEST = sp.cos(2 * x) * t + x**3


def lam(expr):
    return sp.Lambda((x, t), expr)


def s25_source(L):
    r25 = C.piece(L.stmt('S25'), 0).residual()
    th = C.field('theta')
    rest = sp.expand(r25 - r25.coeff(sp.Derivative(th, tau)) * sp.Derivative(th, tau)
                     - r25.coeff(sp.Derivative(th, (xi, 2))) * sp.Derivative(th, (xi, 2)))
    a = r25.coeff(sp.Derivative(th, tau))
    b = r25.coeff(sp.Derivative(th, (xi, 2)))
    # normalise to  theta_tau - theta_xixi - c = 0
    return sp.simplify(-rest / a), sp.simplify(-b / a)


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_uniqueness.py', L, {'T1_test': str(T1_TEST), 'T2_test': str(T2_TEST),
                                             'L_numeric': str(LN), 'kappa_numeric': str(KN)})
    w = C.field('w')

    def s33():
        st = L.stmt('S33')
        r17 = C.piece(L.stmt('S17'), 0).residual()
        pde = C.find_piece(st, lambda p: p.ok and p.exprs[0] == sp.Derivative(w, t), 'w_t = ...')
        a1 = C.plug(r17, 'T', lam(T1_TEST))
        a2 = C.plug(r17, 'T', lam(T2_TEST))
        aw = C.plug(pde.residual(), 'w', lam(T1_TEST - T2_TEST))
        ok_pde, d = C.same_up_to_sign(sp.expand(a1 - a2), aw)
        # data: every S7 condition, subtracted, gives w = 0 at the same points
        exp_pts = set()
        for p in C.pieces(L.stmt('S7')):
            if p.ok and all(r == '=' for r in p.rels):
                for e in p.exprs[:-1]:
                    if isinstance(e, sp.core.function.AppliedUndef) and e.func.__name__ == 'T':
                        exp_pts.add(e.args)
        got = {}
        for p in C.pieces(st):
            if p.ok and all(r == '=' for r in p.rels):
                for e in p.exprs[:-1]:
                    if isinstance(e, sp.core.function.AppliedUndef) and e.func.__name__ == 'w' \
                            and e.args != (x, t):
                        got[e.args] = p.exprs[-1]
        data_ok = set(got) == exp_pts and all(v == 0 for v in got.values())
        wdef = C.find_piece(st, C.lhs_is('w'), 'w = ...').rhs
        def_ok = sp.simplify(wdef - (C.field('T_1') - C.field('T_2'))) == 0
        return ok_pde and data_ok and def_ok, {
            'S17(T1) - S17(T2) vs S33 pde on w = T1 - T2': d, 'q cancels': not aw.has(q),
            'data points from S7': sorted(map(str, exp_pts)),
            'stated w data': {str(kk): str(v) for kk, v in got.items()}, 'w = T1 - T2': def_ok}
    rep.check('S33', L.stmt('S33'), s33)

    def s34():
        st = L.stmt('S34')
        pE = C.find_piece(st, C.lhs_is('E'), 'E(t) = ...')
        pd = C.find_piece(st, lambda p: p.exprs[0] == sp.Derivative(C.field('E'), t), 'dE/dt = ...')
        pde = C.find_piece(L.stmt('S33'), lambda p: p.exprs[0] == sp.Derivative(w, t), 'w_t = ...')
        wt = sp.solve(pde.residual(), sp.Derivative(w, t))[0]
        # dE/dt computed from the E definition; then w_t replaced by S33
        Edef = pE.rhs
        # differentiate under the integral sign (fixed limits 0, L)
        if not Edef.has(sp.Integral):
            return False, f'E is not defined by an integral: {Edef}'
        dE = Edef.replace(lambda e: isinstance(e, sp.Integral),
                          lambda e: sp.Integral(sp.diff(e.function, t), *e.limits))
        dE = dE.xreplace({sp.Derivative(w, t): wt})  # xreplace: x is bound inside the integral
        if dE.has(sp.Derivative(w, t)):
            return False, 'could not replace w_t by S33'
        tests = [x * (Lc - x) * sp.exp(-t) * (1 + x), sp.sin(sp.pi * x / Lc) * sp.cos(t) + x**2 * t]
        res = []
        for wt_ in tests:
            env = {Lc: LN, kappa: KN}
            a = C.plug(dE, 'w', lam(wt_)).subs(env)
            b = C.plug(pd.rhs, 'w', lam(wt_)).subs(env)
            res.append(sp.simplify(a - b))
        return all(r == 0 for r in res), {'dE/dt from E def with S33': C.short(dE), 'stated rhs': C.short(pd.rhs),
                                          'difference on tests': [str(r) for r in res]}
    rep.check('S34', L.stmt('S34'), s34)

    def s35():
        p34 = C.find_piece(L.stmt('S34'), lambda p: p.exprs[0] == sp.Derivative(C.field('E'), t), 'dE/dt')
        p35 = C.piece(L.stmt('S35'), 0)
        if len(p35.exprs) != 3 or p35.rels != ['=', '<='] or p35.exprs[2] != 0:
            return False, f'S35 is not "dE/dt = ... <= 0": {p35.text}'
        tests = [x * (Lc - x) * (1 + x * t), sp.sin(sp.pi * x / Lc) * sp.exp(-t) + sp.sin(2 * sp.pi * x / Lc) * t]
        res, signs = [], []
        for wt_ in tests:
            ends = [wt_.subs(x, 0), sp.simplify(wt_.subs(x, Lc))]
            if any(e != 0 for e in ends):
                return False, 'test field does not satisfy the S33 end data'
            env = {Lc: LN, kappa: KN}
            a = C.plug(p34.rhs, 'w', lam(wt_)).subs(env)
            b = C.plug(p35.exprs[1], 'w', lam(wt_)).subs(env)
            res.append(sp.simplify(a - b))
            signs.append(float(b.subs(t, sp.Rational(1, 2))))
        # control: the identity must FAIL for a field that does not vanish at the ends
        ctrl = x + 1
        env = {Lc: LN, kappa: KN}
        ctrl_diff = sp.simplify(C.plug(p34.rhs, 'w', lam(ctrl)).subs(env) - C.plug(p35.exprs[1], 'w', lam(ctrl)).subs(env))
        ok = all(r == 0 for r in res) and all(s <= 0 for s in signs)
        return ok, {'S34 rhs - S35 rhs on fields with w = 0 at the ends': [str(r) for r in res],
                    'S35 rhs at t = 1/2 (must be <= 0)': signs,
                    'control: difference when w != 0 at the ends (boundary term)': str(ctrl_diff)}
    rep.check('S35', L.stmt('S35'), s35)

    def s36():
        pE = C.find_piece(L.stmt('S34'), C.lhs_is('E'), 'E(t) = ...')
        w0 = [p for p in C.pieces(L.stmt('S33')) if p.ok and p.exprs[0] == C.Fn('w')(x, 0)]
        w0_zero = bool(w0) and w0[0].exprs[-1] == 0
        E0 = sp.simplify(C.plug(pE.rhs, 'w', lam(sp.Integer(0) * x)))
        Epos = [float(C.plug(pE.rhs, 'w', lam(f)).subs(Lc, LN).doit()) for f in (x, sp.sin(x) + 2, -x**2)]
        p35 = C.piece(L.stmt('S35'), 0)
        dec = p35.rels[-1] == '<=' and p35.exprs[-1] == 0
        concl = C.piece(L.stmt('S36'), 0)
        concl_ok = concl.exprs == [C.field('T_1'), C.field('T_2')] and concl.rels == ['=']
        ok = w0_zero and E0 == 0 and all(v > 0 for v in Epos) and dec and concl_ok
        return ok, {'S33 gives w(x,0) = 0': w0_zero, 'E(0) from S34 definition': str(E0),
                    'E > 0 for nonzero test w (E is a positive multiple of int w^2)': Epos,
                    'S35 states dE/dt <= 0': dec, 'concludes T_1 = T_2': concl_ok}
    rep.check('S36', L.stmt('S36'), s36)

    def s37():
        ps = C.pieces(L.stmt('S37'))
        qp = [p for p in ps if p.ok and p.exprs[0] == q]
        Tp = [p for p in ps if p.ok and p.exprs[0] == C.field('T')]
        if not qp or not Tp:
            return False, 'S37 does not parse as q = 0 => T = ...'
        qv, Tv = qp[0].rhs, Tp[0].rhs
        r17 = C.piece(L.stmt('S17'), 0).residual()
        res = sp.simplify(C.plug(r17, 'T', lam(Tv + 0 * x)).subs(q, qv))
        data_ok = True
        for p in C.pieces(L.stmt('S7')):
            if p.ok and all(r == '=' for r in p.rels):
                for e in p.exprs[:-1]:
                    if isinstance(e, sp.core.function.AppliedUndef) and e.func.__name__ == 'T':
                        data_ok &= sp.simplify(Tv - p.exprs[-1]) == 0
        return res == 0 and data_ok, {'q': str(qv), 'T': str(Tv), 'S17 residual': str(res), 'S7 satisfied': data_ok}
    rep.check('S37', L.stmt('S37'), s37)

    # ---------------------------------------------------------- minimum principle
    def s48():
        st = L.stmt('S48')
        p = [pp for pp in C.pieces(st) if pp.ok and len(pp.exprs) == 3]
        attained = bool(p) and p[0].exprs[2].func.__name__ == 'theta'
        s6 = L.stmt('S6')
        closed = 'C\\big([0,L]\\times[0,\\infty)\\big)' in s6.replace(' ', '')
        # numerical: theta from an independent FD solve is continuous and attains its minimum on the
        # closed set (the minimum over a grid of [0,1]x[0,tau1] is finite and attained at a grid point)
        c, b = s25_source(L)
        xg, tg, th = C.solve_theta_fd(float(c), [0.0, 0.05, 0.2], n_int=200)
        mn = float(th.min())
        return attained and closed, {'S48 names a point where the min is attained': attained,
                                     'S6 gives continuity on the closed set [0,L]x[0,inf)': closed,
                                     'FD min of theta on [0,1]x[0,0.2]': mn}
    rep.check('S48', L.stmt('S48'), s48)

    def s49():
        st = L.stmt('S49')
        ps = [p for p in C.pieces(st) if p.ok]
        th = C.field('theta')
        xs, ts, t1 = C.S('xi_star'), C.S('tau_star'), C.S('tau_1')
        pxx = [p for p in ps if p.exprs[0].has(sp.Derivative(C.Fn('theta')(xi, ts), (xi, 2)))]
        pt = [p for p in ps if p.exprs[0].has(sp.Derivative(C.Fn('theta')(xs, tau), tau))]
        loc_x = [p for p in ps if p.exprs[1:2] == [xs]]
        loc_t = [p for p in ps if p.exprs[1:2] == [ts]]
        if not (pxx and pt and loc_x and loc_t):
            return False, 'S49 pieces not found'
        # functions with a negative minimum at an interior xi and at tau* = tau1 (top) or inside
        ex = [((xi - sp.Rational(2, 5))**2 + sp.Rational(3, 10) * (1 - tau) - 1, sp.Rational(2, 5), 1),
              ((xi - sp.Rational(1, 2))**2 * 3 + (tau - sp.Rational(1, 2))**2 - 2, sp.Rational(1, 2), sp.Rational(1, 2)),
              (sp.sin(sp.pi * xi) * (-tau) * (2 - tau), sp.Rational(1, 2), 1)]
        verdicts = []
        for f, xv, tv in ex:
            vxx = float(sp.diff(f, xi, 2).subs({xi: xv, tau: tv}))
            vt = float(sp.diff(f, tau).subs({xi: xv, tau: tv}))
            okxx = sp.sympify(f'{vxx} {pxx[0].rels[0]} {float(pxx[0].rhs)}')
            okt = sp.sympify(f'{vt} {pt[0].rels[0]} {float(pt[0].rhs)}')
            verdicts.append((vxx, vt, bool(okxx), bool(okt)))
        # location: theta = 0 on xi = 0, 1 and tau = 0 (S26): a negative min is not there,
        # but tau* = tau1 is allowed
        loc_ok = loc_x[0].rels == ['<', '<'] and loc_x[0].exprs[0] == 0 and loc_x[0].exprs[2] == 1 \
            and loc_t[0].rels == ['<', '<='] and loc_t[0].exprs[0] == 0 and loc_t[0].exprs[2] == t1
        s26_zero = all(p.exprs[-1] == 0 for p in C.pieces(L.stmt('S26')) if p.ok)
        ok = all(v[2] and v[3] for v in verdicts) and loc_ok and s26_zero
        return ok, {'(theta_xixi, theta_tau, claim_xx holds, claim_tau holds) at known minima': verdicts,
                    'location claims (0<xi*<1, 0<tau*<=tau1)': loc_ok, 'S26 data are all 0': s26_zero}
    rep.check('S49', L.stmt('S49'), s49)

    def s50():
        p = C.piece(L.stmt('S50'), 0)
        c, b = s25_source(L)
        # contradiction needs: at an interior min theta_tau - b*theta_xixi <= 0 (b > 0), but S25 says = c > 0
        contra = bool(sp.sympify(c).is_positive) and bool(sp.sympify(b).is_positive)
        xg, tg, th = C.solve_theta_fd(float(c), [0.001, 0.01, 0.1, 0.5, 2.0], n_int=300)
        mn, mx = float(th.min()), float(th.max())
        rel = p.rels[0]
        holds = {'>=': mn >= -1e-12, '<=': mx <= 1e-12, '>': mn > 0, '<': mx < 0}.get(rel, False) \
            and p.exprs[0] == C.field('theta') and p.exprs[1] == 0
        return contra and holds, {'S25 source c (needs > 0)': str(c), 'S25 diffusion coefficient (needs > 0)': str(b),
                                  'claimed': p.text, 'FD min, max of theta on tau in [0.001, 2]': [mn, mx]}
    rep.check('S50', L.stmt('S50'), s50)

    def s38():
        r25 = C.piece(L.stmt('S25'), 0).residual()
        params25 = sorted(str(s) for s in r25.free_symbols - {xi, tau})
        params26 = sorted({str(s) for p in C.pieces(L.stmt('S26')) if p.ok for e in p.exprs
                           for s in e.free_symbols} - {'xi', 'tau'})
        DTdef = C.find_piece(L.stmt('S20'), C.lhs_is('DeltaT_q'), 'DeltaT_q').rhs
        sign_ok = bool(sp.simplify(DTdef / q).is_positive)
        stated = C.find_piece(L.stmt('S38'), lambda p: p.exprs[0] == C.field('T') - T0, 'T - T_0 = ...')
        rel_ok = sp.simplify(stated.rhs - C.S('DeltaT_q') * C.field('theta')) == 0
        c, _ = s25_source(L)
        xg, tg, th = C.solve_theta_fd(float(c), [0.01, 0.1, 1.0], n_int=200)
        text_ok = 'sign of } q' in L.stmt('S38') or 'sign of q' in L.stmt('S38')
        ok = not params25 and not params26 and sign_ok and rel_ok and float(th.min()) >= -1e-12 and text_ok
        return ok, {'parameters in S25': params25, 'parameters in S26': params26,
                    'DeltaT_q/q > 0 (k, L > 0)': sign_ok, 'T - T0 = DeltaT_q theta': rel_ok,
                    'FD min theta': float(th.min()), 'claims "sign of q"': text_ok}
    rep.check('S38', L.stmt('S38'), s38)

    rep.finish()


if __name__ == '__main__':
    main()
