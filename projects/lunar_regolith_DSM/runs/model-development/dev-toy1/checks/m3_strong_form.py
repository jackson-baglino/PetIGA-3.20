"""M3 check: integral balance -> strong form (S9-S17).

Each derived step is tested by redoing its named operation with SymPy on the
statements parsed from the ledger (--ledger). Integral identities (S9, S10, S11)
are evaluated on concrete smooth test fields with symbolic end points a, b, so a
wrong sign, factor or exponent in the statement leaves a nonzero difference.
Pointwise steps (S12-S17) are compared on undefined functions T(x,t), e, j.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17']

x, t, a, b, q = (C.S(n) for n in ('x', 't', 'a', 'b', 'q'))
k, rho, cp, kappa = (C.S(n) for n in ('k', 'rho', 'c_p', 'kappa'))
E_TESTS = [sp.sin(2 * x + t) * t + x**3 * t**2 + sp.exp(x - 2 * t), sp.cos(x * t) + x * t**3]
J_TESTS = [sp.cos(2 * x + t) + x**2 * t, sp.exp(-x) * sp.sin(t) + x**4]


def lam(expr):
    return sp.Lambda((x, t), expr)


def zero_on_tests(expr, tests):
    """expr contains fields e and/or j; return list of simplified residuals over test pairs."""
    out = []
    for et, jt in tests:
        r = C.plug(C.plug(expr, 'e', lam(et)), 'j', lam(jt))
        out.append(sp.simplify(r))
    return out


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_strong_form.py', L, {'test_e': [str(e) for e in E_TESTS],
                                              'test_j': [str(j) for j in J_TESTS]})
    pairs = list(zip(E_TESTS, J_TESTS))
    Tf, jf, ef = C.field('T'), C.field('j'), C.field('e')

    def s9():
        p = C.piece(L.stmt('S9'), 0)
        res = zero_on_tests(p.residual(), pairs)
        return all(r == 0 for r in res), {'parsed': p.text, 'lhs-rhs on tests': [C.short(r) for r in res]}
    rep.check('S9', L.stmt('S9'), s9)

    def s10():
        p = C.piece(L.stmt('S10'), 0)
        res = zero_on_tests(p.residual(), pairs)
        return all(r == 0 for r in res), {'parsed': p.text, 'lhs-rhs on tests': [C.short(r) for r in res]}
    rep.check('S10', L.stmt('S10'), s10)

    def s11():
        # S2 <=> S11 for every e, j: (S2 lhs - rhs) - (S11 lhs - rhs) must vanish identically
        r2 = C.piece(L.stmt('S2'), 0).residual()
        r11 = C.piece(L.stmt('S11'), 0).residual()
        d1 = zero_on_tests(r2 - r11, pairs)
        d2 = zero_on_tests(r2 + r11, pairs)
        ok = all(r == 0 for r in d1) or all(r == 0 for r in d2)
        return ok, {'(S2 residual) -/+ (S11 residual) on tests': [C.short(r) for r in d1],
                    'sum form': [C.short(r) for r in d2]}
    rep.check('S11', L.stmt('S11'), s11)

    def s12():
        p11 = C.piece(L.stmt('S11'), 0)
        integ = p11.lhs
        if not isinstance(integ, sp.Integral) or p11.rhs != 0:
            return False, f'S11 is not "one integral = 0": {p11.text}'
        r12 = C.piece(L.stmt('S12'), 0).residual()
        if r12.has(sp.Integral):
            return False, 'S12 still contains an integral'
        ok, d = C.same_up_to_sign(integ.function, r12)
        return ok, {'S11 integrand': C.short(integ.function), 'S12 residual': C.short(r12), 'compare': d}
    rep.check('S12', L.stmt('S12'), s12)

    def s13():
        p3 = C.find_piece(L.stmt('S3'), C.lhs_is('e'), 'e = ...')
        e_of_T = p3.rhs
        r12 = C.piece(L.stmt('S12'), 0).residual()
        sub = C.plug(r12, 'e', lam(e_of_T))
        r13 = C.piece(L.stmt('S13'), 0).residual()
        ok, d = C.same_up_to_sign(sub, r13)
        return ok, {'e from S3': C.short(e_of_T), 'S12 with S3': C.short(sub), 'S13': C.short(r13),
                    'compare': d}
    rep.check('S13', L.stmt('S13'), s13)

    def s14():
        p4 = C.find_piece(L.stmt('S4'), C.lhs_is('j'), 'j = ...')
        jx = sp.diff(p4.rhs, x)
        p14 = C.piece(L.stmt('S14'), 0)
        if p14.lhs != sp.Derivative(jf, x):
            return False, f'S14 left side is not dj/dx: {p14.text}'
        ok, d = C.same_up_to_sign(jx, p14.rhs)
        ok = ok and d == 'equal'
        return ok, {'d/dx of S4 right side': C.short(jx), 'S14 right side': C.short(p14.rhs), 'compare': d}
    rep.check('S14', L.stmt('S14'), s14)

    def s15():
        r13 = C.piece(L.stmt('S13'), 0).residual()
        p14 = C.piece(L.stmt('S14'), 0)
        sub = r13.subs(p14.lhs, p14.rhs)
        if sub.has(jf):
            return False, f'j remains after substituting S14 into S13: {sub}'
        r15 = C.piece(L.stmt('S15'), 0).residual()
        ok, d = C.same_up_to_sign(sub, r15)
        return ok, {'S13 with S14': C.short(sub), 'S15 residual': C.short(r15), 'compare': d}
    rep.check('S15', L.stmt('S15'), s15)

    def s16():
        p = C.piece(L.stmt('S16'), 0)
        if p.exprs[0] != kappa or p.rels[0] != '=':
            return False, f'S16 does not define kappa: {p.text}'
        rhs = p.exprs[1]
        pos = sp.simplify(rhs).is_positive
        claims_pos = len(p.rels) > 1 and p.rels[1] == '>' and p.exprs[2] == 0
        # the definition must make S17 the division of S15 by rho c_p (kappa rho c_p = k)
        ok_def = sp.simplify(rhs * rho * cp - k) == 0
        return bool(pos) and claims_pos and ok_def, {
            'kappa definition': str(rhs), 'positive with k, rho, c_p > 0 (S1)': pos,
            'states > 0': claims_pos, 'kappa*rho*c_p == k (makes S17 = S15/(rho c_p))': ok_def}
    rep.check('S16', L.stmt('S16'), s16)

    def s17():
        r15 = C.piece(L.stmt('S15'), 0).residual()
        kdef = C.piece(L.stmt('S16'), 0).exprs[1]
        r17 = C.piece(L.stmt('S17'), 0).residual()
        divided = sp.expand(r15 / (rho * cp))
        r17k = sp.expand(r17.subs(kappa, kdef))
        ok, d = C.same_up_to_sign(divided, r17k)
        return ok, {'S15/(rho c_p)': C.short(divided), 'S17 with kappa from S16': C.short(r17k), 'compare': d}
    rep.check('S17', L.stmt('S17'), s17)

    rep.finish()


if __name__ == '__main__':
    main()
