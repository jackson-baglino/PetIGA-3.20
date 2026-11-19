"""M3 check: nondimensional forms, data mapping and reference numbers (S20-S31).

Chain-rule and substitution steps are tested on concrete test functions
theta(xi, tau) pushed through the ledger's own change of variables (S20 or S27):
both sides of each claim are evaluated on the same test function and compared
in the physical variables x, t. Division steps are tested on undefined functions.
Reference numbers (S31) are recomputed with Pint from the values of choice S19
and compared within 1 %.
"""
import os
import re
import sys

import pint
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26', 'S27', 'S28', 'S29', 'S30', 'S31']

x, t, xi, tau = (C.S(n) for n in ('x', 't', 'xi', 'tau'))
Lc, k, rho, cp, kappa, q, T0 = (C.S(n) for n in ('L', 'k', 'rho', 'c_p', 'kappa', 'q', 'T_0'))
DT, Q, xc, tc = (C.S(n) for n in ('DeltaT_q', 'Q', 'x_c', 't_c'))
TH_TESTS = [sp.sin(2 * xi) * sp.exp(-3 * tau) + xi**3 * tau**2, sp.cos(xi * tau) + xi * tau]
ureg = pint.UnitRegistry()


def s20_defs(L):
    st = L.stmt('S20')
    d = {}
    for name in ('DeltaT_q', 'x_c', 't_c', 'xi', 'tau'):
        d[name] = C.find_piece(st, C.lhs_is(name), name + ' = ...').rhs
    d['T'] = C.find_piece(st, C.lhs_is('T'), 'T = ...').rhs
    return d


def to_xt(expr, defs):
    """Replace xi, tau by their definitions in x, t (S20)."""
    return expr.subs({xi: defs['xi'], tau: defs['tau']}, simultaneous=True)


def T_of(defs, field_name, th_test, T_rhs):
    """T(x,t) built from the ledger's T = ... with a test theta (or theta')."""
    return to_xt(C.plug(T_rhs, field_name, sp.Lambda((xi, tau), th_test)), defs)


def eval_claim(expr, defs, field_name, th_test, T_rhs):
    """Evaluate an expression containing T and theta-fields on the same test function, in x, t."""
    Txt = T_of(defs, field_name, th_test, T_rhs)
    e = C.plug(expr, 'T', sp.Lambda((x, t), Txt))
    e = C.plug(e, field_name, sp.Lambda((xi, tau), th_test))
    return sp.simplify(to_xt(e, defs))


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_scaling.py', L, {'test_theta': [str(f) for f in TH_TESTS]})

    def s20():
        d = s20_defs(L)
        res = {
            'xi == x/x_c': sp.simplify(d['xi'] - x / d['x_c']) == 0,
            'tau == t/t_c': sp.simplify(d['tau'] - t / d['t_c']) == 0,
            'xi in [0,1] for x in [0,L]': (d['xi'].subs(x, 0), d['xi'].subs(x, Lc)) == (0, 1),
            'DeltaT_q = 0 only if q = 0': sp.solve(sp.Eq(d['DeltaT_q'], 0), q) == [0],
            'T = T_0 + DeltaT_q*theta': sp.simplify(d['T'] - (T0 + DT * C.field('theta'))) == 0,
            'q != 0 stated': '\\ne 0' in L.stmt('S20'),
        }
        return all(res.values()), {**res, 'defs': {kk: str(v) for kk, v in d.items()}}
    rep.check('S20', L.stmt('S20'), s20)

    def chain(sid):
        d = s20_defs(L)
        out = []
        for p in C.pieces(L.stmt(sid)):
            if not p.ok:
                return False, f'unparsed piece {p.text}: {p.error}'
            for th in TH_TESTS:
                r = eval_claim(p.residual(), d, 'theta', th, d['T'])
                out.append((p.text, C.short(r)))
        return bool(out) and all(r == '0' for _, r in out), {'lhs-rhs on test theta': out}
    rep.check('S21', L.stmt('S21'), lambda: chain('S21'))
    rep.check('S22', L.stmt('S22'), lambda: chain('S22'))

    def s23():
        d = s20_defs(L)
        r17 = C.piece(L.stmt('S17'), 0).residual()
        r23 = C.piece(L.stmt('S23'), 0).residual()
        out = []
        for th in TH_TESTS:
            a17 = eval_claim(r17, d, 'theta', th, d['T'])
            a23 = eval_claim(r23, d, 'theta', th, d['T'])
            out.append(C.same_up_to_sign(a17, a23))
        return all(o[0] for o in out), {'S17 vs S23 on test theta': out}
    rep.check('S23', L.stmt('S23'), s23)

    def division(sid_from, sid_to, field_name):
        r_a = C.piece(L.stmt(sid_from), 0).residual()
        r_b = C.piece(L.stmt(sid_to), 0).residual()
        dtau = sp.Derivative(C.field(field_name), tau)
        c = sp.simplify(r_a.coeff(dtau) / r_b.coeff(dtau))
        rest = sp.simplify(r_a - c * r_b)
        return c, rest

    def s24():
        c, rest = division('S23', 'S24', 'theta')
        nonzero_for_q_ne_0 = sp.simplify(c.subs(DT, s20_defs(L)['DeltaT_q'])).has(q)
        return rest == 0 and c != 0, {'divisor (ratio of tau-coefficients)': str(c),
                                      'S23 - divisor*S24': str(rest),
                                      'divisor proportional to q (nonzero iff q != 0)': nonzero_for_q_ne_0}
    rep.check('S24', L.stmt('S24'), s24)

    def s25():
        r24 = C.piece(L.stmt('S24'), 0).residual()
        kdef = C.piece(L.stmt('S16'), 0).exprs[1]
        sub = sp.simplify(r24.subs(kappa, kdef).subs(DT, s20_defs(L)['DeltaT_q']))
        r25 = C.piece(L.stmt('S25'), 0).residual()
        ok, dd = C.same_up_to_sign(sub, r25)
        params = sorted(str(s) for s in r25.free_symbols - {xi, tau})
        return ok and not params, {'S24 with S16, S20': C.short(sub), 'S25': C.short(r25), 'compare': dd,
                                   'parameters left in S25': params}
    rep.check('S25', L.stmt('S25'), s25)

    def map_data(sid_data, map_value):
        """Map every S7 data equation through S20 and return {(xi_arg, tau_arg): value}."""
        d = s20_defs(L)
        out = {}
        for p in C.pieces(L.stmt('S7')):
            if not p.ok or any(r != '=' for r in p.rels):
                continue
            val = p.exprs[-1]
            for e in p.exprs[:-1]:
                if isinstance(e, sp.core.function.AppliedUndef) and e.func.__name__ == 'T':
                    xa, ta = e.args
                    xia = xi if xa == x else sp.simplify(d['xi'].subs(x, xa))
                    taua = tau if ta == t else sp.simplify(d['tau'].subs(t, ta))
                    out[(xia, taua)] = sp.simplify(map_value(val))
        return out

    def stated_data(sid, fname):
        out = {}
        for p in C.pieces(L.stmt(sid)):
            if not p.ok or any(r != '=' for r in p.rels):
                continue
            for e in p.exprs[:-1]:
                if isinstance(e, sp.core.function.AppliedUndef) and e.func.__name__ == fname \
                        and e.args != (xi, tau):
                    out[e.args] = p.exprs[-1]
        return out

    def s26():
        d = s20_defs(L)
        theta_of_T = sp.solve(sp.Eq(sp.Symbol('Tv'), d['T'].subs(C.field('theta'), sp.Symbol('th'))),
                              sp.Symbol('th'))[0]
        exp_ = map_data('S7', lambda v: theta_of_T.subs(sp.Symbol('Tv'), v))
        got = stated_data('S26', 'theta')
        ok = set(exp_) == set(got) and all(sp.simplify(exp_[kk] - got[kk]) == 0 for kk in exp_)
        return ok, {'expected from S7 via S20': {str(kk): str(v) for kk, v in exp_.items()},
                    'stated in S26': {str(kk): str(v) for kk, v in got.items()}}
    rep.check('S26', L.stmt('S26'), s26)

    def s27():
        st = L.stmt('S27')
        T_rhs = C.find_piece(st, C.lhs_is('T'), 'T = ...').rhs
        Qdef = C.find_piece(st, C.lhs_is('Q'), 'Q = ...').rhs
        res = {
            'T = T_0*theta\'': sp.simplify(T_rhs - T0 * C.field('thetap')) == 0,
            'Q finite at q = 0': Qdef.subs(q, 0) == 0,
            'S7 gives T_0 > 0': any(p.ok and p.exprs == [T0, 0] and p.rels == ['>'] for p in C.pieces(L.stmt('S7'))),
            'Q * k T_0 / (q L^2) = 1': sp.simplify(Qdef * k * T0 / (q * Lc**2)) == 1,
        }
        return all(res.values()), {kk: v for kk, v in res.items()}
    rep.check('S27', L.stmt('S27'), s27)

    def s28():
        d = s20_defs(L)
        T_rhs = C.find_piece(L.stmt('S27'), C.lhs_is('T'), 'T = ...').rhs
        r17 = C.piece(L.stmt('S17'), 0).residual()
        r28 = C.piece(L.stmt('S28'), 0).residual()
        out = []
        for th in TH_TESTS:
            a17 = eval_claim(r17, d, 'thetap', th, T_rhs)
            a28 = eval_claim(r28, d, 'thetap', th, T_rhs)
            out.append(C.same_up_to_sign(a17, a28))
        return all(o[0] for o in out), {'S17 vs S28 on test theta\'': out}
    rep.check('S28', L.stmt('S28'), s28)

    def s29():
        c, _ = division('S28', 'S29', 'thetap')
        r28 = C.piece(L.stmt('S28'), 0).residual()
        r29 = C.piece(L.stmt('S29'), 0).residual()
        kdef = C.piece(L.stmt('S16'), 0).exprs[1]
        Qdef = C.find_piece(L.stmt('S27'), C.lhs_is('Q'), 'Q = ...').rhs
        rest = sp.simplify((r28 - c * r29).subs(Q, Qdef).subs(kappa, kdef))
        T_rhs = C.find_piece(L.stmt('S27'), C.lhs_is('T'), 'T = ...').rhs
        thp_of_T = sp.solve(sp.Eq(sp.Symbol('Tv'), T_rhs.subs(C.field('thetap'), sp.Symbol('th'))),
                            sp.Symbol('th'))[0]
        exp_ = map_data('S7', lambda v: thp_of_T.subs(sp.Symbol('Tv'), v))
        got = stated_data('S29', 'thetap')
        data_ok = set(exp_) == set(got) and all(sp.simplify(exp_[kk] - got[kk]) == 0 for kk in exp_)
        return rest == 0 and data_ok, {'divisor': str(c), 'S28 - divisor*S29 (Q, kappa substituted)': str(rest),
                                       'data expected': {str(kk): str(v) for kk, v in exp_.items()},
                                       'data stated': {str(kk): str(v) for kk, v in got.items()}}
    rep.check('S29', L.stmt('S29'), s29)

    def s30():
        d = s20_defs(L)
        Qdef = C.find_piece(L.stmt('S27'), C.lhs_is('Q'), 'Q = ...').rhs
        pQ = C.find_piece(L.stmt('S30'), C.lhs_is('Q'), 'Q = ...')
        pth = C.find_piece(L.stmt('S30'), C.lhs_is('thetap'), "theta' = ...")
        r1 = sp.simplify(Qdef - pQ.rhs.subs(DT, d['DeltaT_q']))
        # T0*theta' = T0 + DT*theta  =>  theta' = 1 + (DT/T0) theta ; compare with stated, Q from S27
        T27 = C.find_piece(L.stmt('S27'), C.lhs_is('T'), 'T = ...').rhs
        thp = sp.solve(sp.Eq(T27, d['T']), C.field('thetap'))[0]
        r2 = sp.simplify((thp - pth.rhs.subs(Q, Qdef)).subs(DT, d['DeltaT_q']))
        return r1 == 0 and r2 == 0, {'Q relation residual': str(r1), "theta' relation residual": str(r2)}
    rep.check('S30', L.stmt('S30'), s30)

    def s31():
        st19 = L.stmt('S19')
        vals = {}
        for m in re.finditer(r'([A-Za-z\\]+(?:_\{?[A-Za-z0-9]+\}?)?)\s*=\s*([0-9.]+)\s*\\,\\mathrm\{', st19):
            content, _ = C._balanced(st19, m.end() - 1)
            unit = content.replace('\\,', ' ').replace('{', '').replace('}', '')
            unit = re.sub(r'\^(-?\d+)', r'**\1', unit)
            unit = re.sub(r'([A-Za-z])\s+([A-Za-z])', r'\1*\2', unit)
            unit = re.sub(r'(\d)\s*([A-Za-z])', r'\1*\2', unit)
            name = C._sanitize(m.group(1)).replace('cp', 'c_p')
            vals[name] = float(m.group(2)) * ureg(unit)
        need = {'L', 'k', 'rho', 'c_p', 'q', 'T_0'}
        if not need <= set(vals):
            return False, f'S19 values found: {sorted(vals)}; need {sorted(need)}'
        kd = C.piece(L.stmt('S16'), 0).exprs[1]
        d = s20_defs(L)
        Qd = C.find_piece(L.stmt('S27'), C.lhs_is('Q'), 'Q = ...').rhs

        def ev(expr, extra=None):
            env = dict(vals)
            if extra:
                env.update(extra)
            f = sp.lambdify([C.S(n) for n in env], expr, modules='math')
            return f(*env.values())
        kap = ev(kd)
        computed = {'kappa': kap, 't_c': ev(d['t_c'], {'kappa': kap}), 'DeltaT_q': ev(d['DeltaT_q']),
                    'Q': ev(Qd)}
        st31 = L.stmt('S31')
        stated = {}
        for name in computed:
            p = C.find_piece(st31, C.lhs_is(name), name)
            m = re.search(re.escape({'kappa': '\\kappa', 't_c': 't_c', 'DeltaT_q': '\\Delta T_q', 'Q': 'Q'}[name])
                          + r'\s*=\s*[0-9.]+\s*(\\,\\mathrm\{)?', st31)
            unit = '1'
            if m and m.group(1):
                content, _ = C._balanced(st31, m.end() - 1)
                unit = content.replace('\\,', ' ')
                unit = re.sub(r'\^\{?(-?\d+)\}?', r'**\1', unit)
                unit = re.sub(r'([A-Za-z0-9])\s+([A-Za-z])', r'\1*\2', unit)
            stated[name] = float(p.rhs) * ureg(unit)
        res = {}
        for name, cval in computed.items():
            cval = ureg.Quantity(cval)
            sval = stated[name]
            same_dim = cval.dimensionality == sval.dimensionality
            ratio = float((cval / sval).to_base_units().magnitude) if same_dim else None
            res[name] = {'computed': str(cval.to_base_units()), 'stated': str(sval.to_base_units()),
                         'ratio': ratio, 'within_1pct': bool(same_dim and abs(ratio - 1) <= 0.01)}
        return all(r['within_1pct'] for r in res.values()), res
    rep.check('S31', L.stmt('S31'), s31)

    rep.finish()


if __name__ == '__main__':
    main()
