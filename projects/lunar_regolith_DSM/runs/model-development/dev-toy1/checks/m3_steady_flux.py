"""M3 check: boundary fluxes and the steady state (S39, S41-S45).

Outward fluxes are recomputed from Fourier's law (S4, read from the ledger) with
my own outward normals n(0) = -1, n(L) = +1 on a test field. The steady problem is
re-solved with SymPy (dsolve) from the ledger's own S25-S26 and compared with each
stated step; the dimensional steady state is also checked directly against S15 with
the S7 end values, independent of the scaling.
"""
import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S39', 'S41', 'S42', 'S43', 'S44', 'S45']

x, t, xi, tau = (C.S(n) for n in ('x', 't', 'xi', 'tau'))
Lc, k, q, T0, c0, c1 = (C.S(n) for n in ('L', 'k', 'q', 'T_0', 'c_0', 'c_1'))
F0, FL = C.S('F_0'), C.S('F_L')
T_TEST = T0 + sp.sin(3 * x + t) * x + x**3 * sp.exp(-t)


def j_law(L):
    return C.find_piece(L.stmt('S4'), C.lhs_is('j'), 'j = ...').rhs


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_steady_flux.py', L, {'test_T': str(T_TEST)})
    ths = C.field('theta_s')

    def s39():
        jT = C.plug(j_law(L), 'T', sp.Lambda((x, t), T_TEST))
        mine = {F0: sp.simplify(jT.subs(x, 0) * (-1)), FL: sp.simplify(jT.subs(x, Lc) * (+1))}
        out, ok = [], True
        for p in C.pieces(L.stmt('S39')):
            if not p.ok:
                return False, f'unparsed: {p.text} ({p.error})'
            name = p.exprs[0]
            if name not in mine:
                return False, f'unexpected left side {name}'
            for e in p.exprs[1:]:
                v = sp.simplify(C.plug(C.plug(e, 'T', sp.Lambda((x, t), T_TEST)), 'j',
                                       sp.Lambda((x, t), C.plug(j_law(L), 'T', sp.Lambda((x, t), T_TEST)))))
                good = sp.simplify(v - mine[name]) == 0
                ok &= good
                out.append((str(name), C.short(e), good))
        ok &= {r[0] for r in out} == {'F_0', 'F_L'}
        return ok, {'j n with n(0)=-1, n(L)=+1 on test T': {str(kk): str(v) for kk, v in mine.items()},
                    'each stated member matches': out}
    rep.check('S39', L.stmt('S39'), s39)

    def s41():
        r25 = C.piece(L.stmt('S25'), 0).residual()
        steady = r25.subs(sp.Derivative(C.field('theta'), tau), 0)
        steady = steady.subs(sp.Derivative(C.field('theta'), (xi, 2)), sp.Derivative(ths, (xi, 2)))
        p41 = C.find_piece(L.stmt('S41'), lambda p: p.exprs[0].has(sp.Derivative), "theta_s'' = ...")
        ok, d = C.same_up_to_sign(steady, p41.residual())
        ends26 = {}
        for p in C.pieces(L.stmt('S26')):
            for e in p.exprs[:-1]:
                if e.args[1] == tau:
                    ends26[e.args[0]] = p.exprs[-1]
        ends41 = {}
        for p in C.pieces(L.stmt('S41')):
            if p.ok and not p.exprs[0].has(sp.Derivative):
                for e in p.exprs[:-1]:
                    ends41[e.args[0]] = p.exprs[-1]
        ends_ok = ends26 == ends41 and set(ends41) == {0, 1}
        return ok and ends_ok, {'S25 with theta_tau = 0': C.short(steady), 'S41': C.short(p41.residual()),
                                'compare': d, 'S26 end values': {str(kk): str(v) for kk, v in ends26.items()},
                                'S41 end values': {str(kk): str(v) for kk, v in ends41.items()}}
    rep.check('S41', L.stmt('S41'), s41)

    def s42():
        p41 = C.find_piece(L.stmt('S41'), lambda p: p.exprs[0].has(sp.Derivative), "theta_s''")
        ode = p41.residual()
        sol = C.piece(L.stmt('S42'), 0).rhs
        sat = sp.simplify(ode.subs(ths, sol).doit()) == 0
        consts = sorted(str(s) for s in sol.free_symbols - {xi})
        gen = sp.dsolve(ode, ths).rhs
        # the stated family must equal the general solution for some choice of constants
        C1, C2 = sorted(gen.free_symbols - {xi}, key=str)
        same_family = sp.solve(sp.Poly(sp.expand(gen - sol), xi).coeffs(), [C1, C2]) != []
        return sat and len(consts) == 2 and same_family, {
            'S42 satisfies S41 ODE': sat, 'free constants': consts, 'dsolve': str(gen),
            'same two-parameter family': same_family}
    rep.check('S42', L.stmt('S42'), s42)

    def s43():
        ps = [p for p in C.pieces(L.stmt('S43')) if p.ok]
        sol = [p for p in ps if p.exprs[0] == ths][0].rhs
        mine = sp.dsolve(C.find_piece(L.stmt('S41'), lambda p: p.exprs[0].has(sp.Derivative), "ode").residual(),
                         ths, ics={C.Fn('theta_s')(0): 0, C.Fn('theta_s')(1): 0}).rhs
        same = sp.simplify(sol - mine) == 0
        crit = sp.solve(sp.diff(mine, xi), xi)
        mx = [sp.simplify(mine.subs(xi, c_)) for c_ in crit]
        pmax = [p for p in ps if p.exprs[0].func.__name__.startswith('Max_')]
        # evaluate theta_s(1/2) with my solution, compare every member of the chain with my max
        chain_ok = bool(pmax)
        details_chain = []
        for e in (pmax[0].exprs[1:] if pmax else []):
            v = sp.simplify(C.plug(e, 'theta_s', sp.Lambda(xi, mine)))
            details_chain.append(str(v))
            chain_ok &= len(mx) == 1 and sp.simplify(v - mx[0]) == 0
        argmax_ok = any(isinstance(e, sp.core.function.AppliedUndef) and e.args == (crit[0],)
                        for e in (pmax[0].exprs[1:] if pmax else []))
        return same and chain_ok and argmax_ok, {'my steady solution': str(mine), 'stated': str(sol),
                                                 'argmax': [str(c_) for c_ in crit], 'max': [str(m) for m in mx],
                                                 'stated chain evaluated': details_chain,
                                                 'stated argmax matches': argmax_ok}
    rep.check('S43', L.stmt('S43'), s43)

    def s44():
        Ts = C.piece(L.stmt('S44'), 0).rhs
        # (i) back-substitution of S20 into S43
        d = {n: C.find_piece(L.stmt('S20'), C.lhs_is(n), n).rhs for n in ('DeltaT_q', 'xi')}
        T_map = C.find_piece(L.stmt('S20'), C.lhs_is('T'), 'T = ...').rhs
        th43 = [p for p in C.pieces(L.stmt('S43')) if p.ok and p.exprs[0] == ths][0].rhs
        back = T_map.subs(C.field('theta'), th43).subs(xi, d['xi']).subs(C.S('DeltaT_q'), d['DeltaT_q'])
        r1 = sp.simplify(back - Ts)
        # (ii) direct: steady S15 and S7 end values
        r15 = C.piece(L.stmt('S15'), 0).residual()
        res15 = sp.simplify(C.plug(r15, 'T', sp.Lambda((x, t), Ts)))
        ends = [sp.simplify(Ts.subs(x, 0) - T0), sp.simplify(Ts.subs(x, Lc) - T0)]
        q0 = sp.simplify(Ts.subs(q, 0) - T0) == 0
        ok = r1 == 0 and res15 == 0 and all(e == 0 for e in ends) and q0
        return ok, {'back-substitution - stated': str(r1), 'steady S15 residual': str(res15),
                    'end values - T0': [str(e) for e in ends], 'q = 0 gives T0 (S37)': q0}
    rep.check('S44', L.stmt('S44'), s44)

    def s45():
        Ts = C.piece(L.stmt('S44'), 0).rhs
        j_on_Ts = C.plug(j_law(L), 'T', sp.Lambda((x, t), Ts + 0 * t))
        vals = {}
        for p in C.pieces(L.stmt('S39')):
            name = p.exprs[0]
            vals[name] = sp.simplify(C.plug(C.plug(p.exprs[-1], 'T', sp.Lambda((x, t), Ts + 0 * t)), 'j',
                                            sp.Lambda((x, t), j_on_Ts)))
        out, ok = [], True
        for p in C.pieces(L.stmt('S45')):
            if not p.ok:
                return False, f'unparsed {p.text}'
            ev = [sp.simplify(e.subs(vals).doit()) for e in p.exprs]
            good = all(sp.simplify(e - ev[0]) == 0 for e in ev)
            ok &= good
            out.append((p.text, [str(e) for e in ev], good))
        balance = sp.simplify(vals[F0] + vals[FL] - sp.integrate(q, (x, 0, Lc)))
        return ok and balance == 0, {'F0, FL from S39 on S44': [str(vals[F0]), str(vals[FL])], 'relations': out,
                                     'F0 + FL - int q': str(balance)}
    rep.check('S45', L.stmt('S45'), s45)

    rep.finish()


if __name__ == '__main__':
    main()
