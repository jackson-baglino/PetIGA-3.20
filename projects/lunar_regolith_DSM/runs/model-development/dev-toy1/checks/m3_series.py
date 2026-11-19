"""M3 check: transient, series solution, bounds and dominant terms (S51-S57, S32).

Independent reference: my own sine series of the transient, with coefficients
computed here by SymPy from my own steady solution of theta'' = -1, theta(0)=theta(1)=0,
cross-checked against a finite-difference solve of the ledger's S25-S26
(m3_common.solve_theta_fd). Every ledger statement (series, coefficients, bound
formulas, numbers and inequality directions) is parsed from --ledger and evaluated
against that reference.
"""
import os
import re
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S51', 'S52', 'S53', 'S54', 'S55', 'S56', 'S57', 'S32']

xi, tau, ell, r = (C.S(n) for n in ('xi', 'tau', 'ell', 'r'))
R_TERMS = 3000          # odd and even modes up to ell = R_TERMS
XI = np.linspace(0.0, 1.0, 2001)


# ---------------------------------------------------------------- my reference
def my_reference():
    ths = C.field('theta_s')
    mine = sp.dsolve(sp.Derivative(ths, (xi, 2)) + 1, ths,
                     ics={C.Fn('theta_s')(0): 0, C.Fn('theta_s')(1): 0}).rhs
    beta = sp.simplify(2 * sp.integrate(mine * sp.sin(ell * sp.pi * xi), (xi, 0, 1)))
    fb = sp.lambdify(ell, beta, 'numpy')
    ells = np.arange(1, R_TERMS + 1, dtype=float)
    b = np.array([float(fb(int(e))) for e in ells])
    return mine, beta, ells, b


class Ref:
    def __init__(self):
        self.ths_expr, self.beta_expr, self.ells, self.b = my_reference()
        self.ths = sp.lambdify(xi, self.ths_expr, 'numpy')
        self.S = np.sin(np.outer(self.ells * np.pi, XI))          # [modes, xi]

    def decay(self, tv):
        return np.exp(-(self.ells * np.pi) ** 2 * tv)

    def v(self, tv):
        return (self.b * self.decay(tv)) @ self.S

    def theta(self, tv):
        return self.ths(XI) - self.v(tv)

    def theta_t(self, tv):
        return (self.b * (self.ells * np.pi) ** 2 * self.decay(tv)) @ self.S

    def theta_xx(self, tv):
        return -1.0 + (self.b * (self.ells * np.pi) ** 2 * self.decay(tv)) @ self.S

    def at_half(self, tv):
        s = np.sin(self.ells * np.pi / 2)
        return float(np.sum(self.b * (self.ells * np.pi) ** 2 * self.decay(tv) * s))


def numeric_quantity(expr, ref, tv):
    """Evaluate a parsed expression built from theta, theta_s and their derivatives on XI."""
    th, ths = C.field('theta'), C.field('theta_s')
    V, B, Cc, D = sp.symbols('V_v B_tht C_thxx D_ths')
    e = expr.xreplace({sp.Derivative(th, tau): B, sp.Derivative(th, (xi, 2)): Cc})
    # theta = theta_s - v, written out so that theta - theta_s is evaluated as -v (no cancellation)
    e = e.xreplace({th: D - V}).xreplace({ths: D})
    if e.has(sp.Derivative) or e.has(sp.core.function.AppliedUndef):
        raise ValueError(f'cannot evaluate {expr}')
    f = sp.lambdify((V, B, Cc, D), e, 'numpy')
    out = f(ref.v(tv), ref.theta_t(tv), ref.theta_xx(tv), ref.ths(XI))
    return np.broadcast_to(out, XI.shape)


def check_max_claims(L, sid, ref):
    """For every piece 'Max_xi(expr) REL number' under a preceding tau-range condition,
    test it on a tau grid with my reference."""
    ps = C.pieces(L.stmt(sid))
    cond, out, ok, n = None, [], True, 0
    for p in ps:
        if not p.ok:
            continue
        if tau in p.exprs and all(not e.has(sp.Function) for e in p.exprs):
            cond = p
            continue
        lhs = p.exprs[0]
        if not (isinstance(lhs, sp.core.function.AppliedUndef) and lhs.func.__name__ == 'Max_xi'):
            continue
        if len(p.exprs) != 2 or not p.exprs[1].is_number:
            continue
        lo, hi = tau_range(cond)
        grid = np.unique(np.concatenate([np.geomspace(max(lo, 1e-4), hi, 60), [hi]])) if np.isfinite(hi) \
            else np.geomspace(lo, lo * 5, 40)
        vals = np.array([float(np.max(numeric_quantity(lhs.args[0], ref, tv))) for tv in grid])
        num = float(p.exprs[1])
        rel = p.rels[0]
        if rel in ('>=', '>'):
            worst = float(vals.min())
            good = worst >= num if rel == '>=' else worst > num
        else:
            worst = float(vals.max())
            good = worst <= num if rel == '<=' else worst < num
        ok &= good
        n += 1
        out.append({'claim': p.text, 'tau range': [lo, hi], 'worst over range': worst, 'holds': bool(good)})
    return ok and n > 0, out


def tau_range(cond):
    if cond is None:
        raise ValueError('no tau range before the claim')
    ex, rels = cond.exprs, cond.rels
    if len(ex) == 3:            # a < tau <= b
        return float(ex[0]), float(ex[2])
    if ex[0] == tau and rels[0] in ('>', '>='):
        return float(ex[1]), np.inf
    raise ValueError(f'unsupported condition {cond.text}')


def series_summand(p):
    s = [a for a in sp.preorder_traversal(p.rhs) if isinstance(a, sp.Sum)]
    if len(s) != 1:
        raise ValueError(f'no single sum in {p.text}')
    return s[0], sp.simplify(p.rhs - s[0])


def eval_sum_numeric(summand, tv, n_terms=R_TERMS // 2):
    f = sp.lambdify((r, xi, tau), summand, 'numpy')
    rr = np.arange(n_terms, dtype=float)[:, None]
    return np.sum(f(rr, XI[None, :], tv), axis=0)


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_series.py', L, {'R_TERMS': R_TERMS, 'xi_grid_points': len(XI)})
    ref = Ref()
    th, ths, v = C.field('theta'), C.field('theta_s'), C.field('v')

    def s51():
        ps = [p for p in C.pieces(L.stmt('S51')) if p.ok]
        vdef = [p for p in ps if p.exprs[0] == v][0].rhs
        pde = [p for p in ps if p.exprs[0] == sp.Derivative(v, tau)][0].residual()
        th43 = [p for p in C.pieces(L.stmt('S43')) if p.ok and p.exprs[0] == ths][0].rhs
        r25 = C.piece(L.stmt('S25'), 0).residual()
        tests = [sp.sin(2 * xi) * sp.exp(-3 * tau) + xi**3 * tau**2, sp.cos(xi * tau)]
        res = []
        for tt in tests:
            vv = C.plug(vdef, 'theta_s', sp.Lambda(xi, th43))
            vv = C.plug(vv, 'theta', sp.Lambda((xi, tau), tt))
            a51 = C.plug(pde, 'v', sp.Lambda((xi, tau), vv))
            a25 = C.plug(r25, 'theta', sp.Lambda((xi, tau), tt))
            res.append(C.same_up_to_sign(a51, a25))
        # data: v = theta_s - theta at the S26 points, theta_s from S43
        exp_ = {}
        for p in C.pieces(L.stmt('S26')):
            for e in p.exprs[:-1]:
                exp_[e.args] = sp.simplify(th43.subs(xi, e.args[0]) - p.exprs[-1])
        got = {}
        for p in ps:
            for e in p.exprs[:-1]:
                if isinstance(e, sp.core.function.AppliedUndef) and e.func.__name__ == 'v' and e.args != (xi, tau):
                    got[e.args] = p.exprs[-1]
        data_ok = set(exp_) == set(got) and all(sp.simplify(exp_[kk] - got[kk]) == 0 for kk in exp_)
        def_ok = sp.simplify(vdef - (ths - th)) == 0
        return all(x_[0] for x_ in res) and data_ok and def_ok, {
            'S51 pde on v = theta_s - theta_test vs S25 on theta_test': res, 'v = theta_s - theta': def_ok,
            'data expected (S26, S43)': {str(kk): str(vv) for kk, vv in exp_.items()},
            'data stated': {str(kk): str(vv) for kk, vv in got.items()}}
    rep.check('S51', L.stmt('S51'), s51)

    def s52():
        ps = [p for p in C.pieces(L.stmt('S52')) if p.ok]
        if len(ps) != 2:
            return False, 'S52 should have two integrals'
        out, ok = [], True
        for p in ps:
            for lv in range(1, 11):
                exact = sp.simplify(p.lhs.subs(ell, lv).doit())
                stated = sp.simplify(p.rhs.subs(ell, lv))
                good = sp.simplify(exact - stated) == 0
                ok &= good
                if not good:
                    out.append((p.text, lv, str(exact), str(stated)))
        sym = [sp.simplify(sp.simplify(p.lhs.doit()) - p.rhs) for p in ps]
        return ok, {'mismatches for ell = 1..10': out, 'symbolic difference (integer ell)': [str(s) for s in sym]}
    rep.check('S52', L.stmt('S52'), s52)

    def s53():
        p = [pp for pp in C.pieces(L.stmt('S53')) if pp.ok][0]
        raw = L.stmt('S53')
        m = re.search(r'\\begin\{cases\}(.*?)\\end\{cases\}', raw, re.S)
        cases = []
        if m:
            for row in m.group(1).split('\\\\'):
                val, _, cond = row.partition('&')
                cases.append((C.parse_expr_str(val.strip()), 'odd' in cond, 'even' in cond))
        out, ok = [], True
        for lv in range(1, 13):
            mine = sp.simplify(ref.beta_expr.subs(ell, lv))
            for e in p.exprs[1:]:
                val = sp.simplify(C.plug(e.subs(ell, lv), 'theta_s', sp.Lambda(xi, ref.ths_expr)))
                good = sp.simplify(val - mine) == 0
                ok &= good
                if not good:
                    out.append((lv, C.short(e), str(val), str(mine)))
            for val, odd, even in cases:
                if (lv % 2 == 1 and odd) or (lv % 2 == 0 and even):
                    good = sp.simplify(val.subs(ell, lv) - mine) == 0
                    ok &= good
                    if not good:
                        out.append((lv, 'case ' + str(val), str(val.subs(ell, lv)), str(mine)))
        ok &= len(cases) == 2
        return ok, {'my beta_ell': str(ref.beta_expr), 'cases parsed': [str(c_[0]) for c_ in cases],
                    'mismatches for ell = 1..12': out}
    rep.check('S53', L.stmt('S53'), s53)

    def s54():
        ps = [p for p in C.pieces(L.stmt('S54')) if p.ok]
        pv = [p for p in ps if p.exprs[0] == v][0]
        pth = [p for p in ps if p.exprs[0] == th][0]
        sm, rest = series_summand(pv)
        summand = sm.function
        # (i) each term solves the S51 pde and end data
        pde = [p for p in C.pieces(L.stmt('S51')) if p.ok and p.exprs[0] == sp.Derivative(v, tau)][0].residual()
        term_pde = sp.simplify(C.plug(pde, 'v', sp.Lambda((xi, tau), summand)))
        ends = [sp.simplify(summand.subs(xi, 0)), sp.simplify(summand.subs(xi, 1))]
        # (ii) tau = 0: partial sum vs theta_s (uniform convergence)
        v0 = eval_sum_numeric(summand, 0.0) + float(rest) if rest.is_number else None
        ic_err = float(np.max(np.abs(v0 - ref.ths(XI)))) if v0 is not None else None
        # (iii) theta = theta_s - v vs FD of the ledger's S25-S26 and vs my series
        r25 = C.piece(L.stmt('S25'), 0).residual()
        c = -r25.subs({sp.Derivative(th, tau): 0, sp.Derivative(th, (xi, 2)): 0})
        taus = [0.01, 0.05, 0.2, 1.0]
        xg, tg, thfd = C.solve_theta_fd(float(c), taus, n_int=799)
        errs_fd, errs_ref = [], []
        th_rel = pth.rhs
        for i, tv in enumerate(taus):
            vs = eval_sum_numeric(summand, tv) + float(rest)
            val = sp.lambdify((C.S('Vv'), C.S('Ss')), th_rel.xreplace({v: C.S('Vv'), ths: C.S('Ss')}),
                              'numpy')(vs, ref.ths(XI))
            errs_fd.append(float(np.max(np.abs(np.interp(xg, XI, val) - thfd[i]))))
            errs_ref.append(float(np.max(np.abs(val - ref.theta(tv)))))
        ok = term_pde == 0 and all(e == 0 for e in ends) and ic_err is not None and ic_err < 1e-6 \
            and max(errs_fd) < 2e-5 and max(errs_ref) < 1e-9
        return ok, {'term residual in S51 pde': str(term_pde), 'term end values': [str(e) for e in ends],
                    'max |v(xi,0) - theta_s| (partial sum)': ic_err,
                    'max |theta_series - theta_FD| at tau=0.01,0.05,0.2,1': errs_fd,
                    'max |theta_series - my series|': errs_ref, 'tolerances': '1e-6 (IC), 2e-5 (FD h=1/800), 1e-9'}
    rep.check('S54', L.stmt('S54'), s54)

    def s55():
        p54 = [p for p in C.pieces(L.stmt('S54')) if p.ok and p.exprs[0] == v][0]
        p55 = C.piece(L.stmt('S55'), 0)
        s54, _ = series_summand(p54)
        s55, rest55 = series_summand(p55)
        # theta = theta_s - v  =>  d theta / d tau = - d v / d tau, term by term
        term = sp.simplify(-sp.diff(s54.function, tau) - s55.function)
        lhs_ok = p55.lhs == sp.Derivative(th, tau)
        # numeric: compare with my reference derivative and with a time difference of the FD solve
        errs = [float(np.max(np.abs(eval_sum_numeric(s55.function, tv) + float(rest55) - ref.theta_t(tv))))
                for tv in (0.01, 0.1, 0.5)]
        h = 1e-4
        xg, tg, thfd = C.solve_theta_fd(1.0, [0.1 - h, 0.1 + h], n_int=799)
        fd_dt = (thfd[1] - thfd[0]) / (2 * h)
        err_fd = float(np.max(np.abs(np.interp(xg, XI, eval_sum_numeric(s55.function, 0.1)) - fd_dt)))
        ok = term == 0 and lhs_ok and rest55 == 0 and max(errs) < 1e-8 and err_fd < 1e-3
        return ok, {'-(d/dtau S54 term) - S55 term': str(term), 'lhs is d theta/d tau': lhs_ok,
                    'max |S55 - my theta_tau| at tau=0.01,0.1,0.5': errs,
                    'max |S55 - FD time difference| at tau=0.1': err_fd}
    rep.check('S55', L.stmt('S55'), s55)

    def s56():
        p = [pp for pp in C.pieces(L.stmt('S56')) if pp.ok and len(pp.exprs) == 4][0]
        cond = [pp for pp in C.pieces(L.stmt('S56')) if pp.ok and tau in pp.exprs][0]
        lo, hi = tau_range(cond)
        # my g: first two terms of the alternating series at xi = 1/2
        ab = ref.b * (ref.ells * np.pi) ** 2 * np.sin(ref.ells * np.pi / 2)
        g_mine = lambda tv: float(ab[0] * np.exp(-np.pi**2 * tv) + ab[2] * np.exp(-9 * np.pi**2 * tv))  # noqa: E731
        # g as defined in the symbols table, parsed
        gm = L.symbol_meaning('g(\\tau)') or ''
        mm = re.search(r'lower bound (.*?) for', gm)
        g_tab = C.parse_expr_str(mm.group(1).replace('exp(', '\\exp(')) if mm else None
        grid = np.concatenate([np.geomspace(max(lo, 1e-4), hi, 200), [hi]])
        rows, ok = [], True
        for tv in grid:
            vals = [float(np.max(ref.theta_t(tv))), ref.at_half(tv), g_mine(tv), float(p.exprs[3])]
            # member 2 must be d theta/d tau at xi = 1/2 and member 1 a max over xi
            for a_, b_, rel in zip(vals[:-1], vals[1:], p.rels):
                good = (a_ >= b_ - 1e-12) if rel == '>=' else (a_ > b_) if rel == '>' else False
                ok &= good
            if g_tab is not None:
                ok &= abs(float(g_tab.subs(tau, tv)) - g_mine(tv)) < 1e-12
            rows.append(vals)
        member_ok = p.exprs[1] == sp.Derivative(C.Fn('theta')(sp.Rational(1, 2), tau), tau) and \
            p.exprs[0].func.__name__ == 'Max_xi'
        # tau_g: the only critical point of g
        cex = [sp.simplify(ref.beta_expr.subs(ell, lv) * (lv * sp.pi)**2 * sp.sin(lv * sp.pi / 2)) for lv in (1, 3)]
        g_exact = cex[0] * sp.exp(-sp.pi**2 * tau) + cex[1] * sp.exp(-9 * sp.pi**2 * tau)
        crit = [s_ for s_ in sp.solve(sp.diff(g_exact, tau), tau) if s_.is_real]
        if len(crit) != 1:
            return False, f'g has {len(crit)} real critical points: {crit}'
        tg_mine = float(crit[0])
        tgm = L.symbol_meaning('\\tau_g') or ''
        tg_ok = abs(np.log(3) / (8 * np.pi**2) - tg_mine) < 1e-10 and 'ln 3/(8\\pi^2)' in tgm
        worst = min(rr[2] for rr in rows)
        return ok and member_ok and g_tab is not None and tg_ok, {
            'tau range': [lo, hi], 'min over range of g (mine)': worst, 'stated bound': float(p.exprs[3]),
            'g in symbols table equals mine': g_tab is not None, 'tau_g mine': tg_mine, 'tau_g as stated': tg_ok,
            'first/last row [max_xi theta_t, theta_t(1/2), g, bound]': [rows[0], rows[-1]]}
    rep.check('S56', L.stmt('S56'), s56)

    def s57():
        ps = [p for p in C.pieces(L.stmt('S57')) if p.ok]
        formulas = [p for p in ps if p.exprs[0].func.__name__ == 'Max_xi' and not p.exprs[1].is_number]
        if len(formulas) != 2:
            return False, 'expected two bound formulas'
        grid = np.geomspace(0.002, 3.0, 120)
        rows, ok = [], True
        for p in formulas:
            fb = sp.lambdify(tau, p.exprs[1], 'numpy')
            worst = -np.inf
            for tv in grid:
                actual = float(np.max(numeric_quantity(p.exprs[0].args[0], ref, tv)))
                bound = float(fb(tv))
                worst = max(worst, actual / bound)
            dec = bool(np.all(np.diff([fb(tv) for tv in np.linspace(1, 10, 200)]) < 0))
            ok &= worst <= 1.0 + 1e-12 and dec   # 1e-12: floating-point roundoff only
            rows.append({'bound': C.short(p.exprs[1]), 'max actual/bound on [0.002, 3]': worst,
                         'bound decreasing on [1, 10]': dec, 'bound at tau=1': float(fb(1.0))})
        # the tau >= 1 numbers: true value and 1 % agreement with the bound at tau = 1
        num_ok, nums = check_max_claims(L, 'S57', ref)
        for p in ps:
            if p.exprs[0].func.__name__ == 'Max_xi' and p.exprs[1].is_number:
                f = [ff for ff in formulas if ff.exprs[0] == p.exprs[0]][0]
                b1 = float(f.exprs[1].subs(tau, 1))
                within = abs(float(p.exprs[1]) / b1 - 1) <= 0.01 and float(p.exprs[1]) >= b1
                ok &= within
                nums.append({'number': float(p.exprs[1]), 'bound at tau=1': b1, 'within 1% and >= bound': within})
        return ok and num_ok, {'bounds': rows, 'tau>=1 numbers': nums}
    rep.check('S57', L.stmt('S57'), s57)

    def s32():
        ok, rows = check_max_claims(L, 'S32', ref)
        c, b = None, None
        r25 = C.piece(L.stmt('S25'), 0).residual()
        coefs = [r25.coeff(sp.Derivative(th, tau)), -r25.coeff(sp.Derivative(th, (xi, 2))),
                 -r25.subs({sp.Derivative(th, tau): 0, sp.Derivative(th, (xi, 2)): 0})]
        r29 = [p for p in C.pieces(L.stmt('S29')) if p.ok][0].residual()
        Qv = [p for p in C.pieces(L.stmt('S31')) if p.ok and p.exprs[0] == C.S('Q')][0].rhs
        thp = C.field('thetap')
        coefs29 = [r29.coeff(sp.Derivative(thp, tau)), -r29.coeff(sp.Derivative(thp, (xi, 2))),
                   -r29.subs({sp.Derivative(thp, tau): 0, sp.Derivative(thp, (xi, 2)): 0}).subs(C.S('Q'), Qv)]
        coef_ok = all(cc == 1 for cc in coefs) and all(cc == 1 for cc in coefs29)
        # conduction at the ends: theta_xixi(0, tau) = -1 for every tau > 0
        ends = [float(ref.theta_xx(tv)[0]) for tv in (1e-3, 0.05, 1.0)]
        return ok and coef_ok and all(abs(e + 1) < 1e-9 for e in ends), {
            'claims': rows, 'S25 coefficients (storage, conduction, source)': [str(cc) for cc in coefs],
            'S29 coefficients at Q from S31': [str(cc) for cc in coefs29], 'theta_xixi(0, tau)': ends}
    rep.check('S32', L.stmt('S32'), s32)

    rep.finish()


if __name__ == '__main__':
    main()
