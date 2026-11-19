"""M3 check: reduction to today's model and the 'inspection' steps (S18, S40, S46, S47).

Reads the heat block R[a][1] and the boundary block of Residual_A1 from
src/assembly.c (--src, main), and the mixture laws from src/material_properties.c.
The weak heat residual is turned into its 1D strong form (integration by parts,
checked on test functions with N = 0 at Dirichlet ends; phase field uniform and
constant in time; Cartesian weight rw = 1) and compared term by term with S15.
"""
import os
import re
import sys

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S18', 'S40', 'S46', 'S47']

x, t, Lc, k, rho, cp, q = (C.S(n) for n in ('x', 't', 'L', 'k', 'rho', 'c_p', 'q'))
GRAIN_ICE = ['psi', 'phi', 'eps', 'eps_s', 'epsilon', 'varepsilon']


def code_heat_weak(src_text):
    m = re.search(r'R\[a\]\[1\]\s*=\s*rw\s*\*\s*\((.*?)\);', src_text, re.S)
    if not m:
        raise ValueError('R[a][1] = rw * ( ... ); not found in assembly.c')
    line = src_text[:m.start()].count('\n') + 1
    body = m.group(1).replace('user->', '').replace('N0[a]', 'N')
    body = re.sub(r'PetscRealPart\(([^()]*)\)', r'\1', body)
    body = ' '.join(body.split())
    names = ['rho', 'cp', 'N', 'tem_t', 'xi_T', 'thcond', 'gN_gtem', 'pc', 'rho_ice', 'lat_sub', 'phi_t']
    loc = {nm: sp.Symbol(nm) for nm in names}
    W = parse_expr(body, local_dict=loc, transformations=standard_transformations)
    return W, loc, line, body


def main():
    args = C.cli(__doc__, with_src=True)
    L = C.Ledger(args.ledger)
    asm = os.path.join(args.src, 'assembly.c')
    mat = os.path.join(args.src, 'material_properties.c')
    rep = C.Reporter('m3_reduction.py', L, {'assembly_c': asm, 'material_properties_c': mat})
    src = open(asm).read()
    mat_src = open(mat).read()
    Tf = C.field('T')

    def r15():
        return C.piece(L.stmt('S15'), 0).residual()

    def s18():
        r = r15()
        names = {str(s) for s in r.free_symbols}
        present = [nm for nm in GRAIN_ICE if nm in names]
        psi = C.S('psi')
        unchanged = sp.simplify(r.subs(psi, 1) - r) == 0
        s8 = L.stmt('S8')
        s8_absent = all(tok in s8 for tok in ('\\psi', '\\phi')) and 'psi\\equiv1' in s8.replace(' ', '')
        return not present and unchanged and s8_absent, {
            'grain/ice symbols in S15': present, 'S15|psi=1 == S15': unchanged,
            'S8 states psi, phi absent (psi == 1)': s8_absent}
    rep.check('S18', L.stmt('S18'), s18)

    def reduction_vs_code():
        W, loc, line, body = code_heat_weak(src)
        N, g = loc['N'], loc['gN_gtem']
        # weak form must be A*N + B*(dN . dT): split it
        B = sp.diff(W, g)
        A = sp.expand((W - B * g) / N)
        if A.has(N) or A.has(g) or B.has(N) or B.has(g):
            return False, f'weak heat residual is not A*N + B*dN.dT: {W}'
        # uniform, steady phase field: phi_t = 0; rho, cp, k are linear mixtures in phi (constant)
        mixtures = all(re.search(pat, mat_src) for pat in (
            r'\(\*cond\)\s*=\s*ice\s*\*\s*cond_ice\s*\+\s*air\s*\*\s*cond_air',
            r'\(\*cp\)\s*=\s*ice\s*\*\s*cp_ice\s*\+\s*air\s*\*\s*cp_air',
            r'\(\*rho\)\s*=\s*ice\s*\*\s*rho_ice\s*\+\s*air\s*\*\s*rho_air'))
        sub = {loc['phi_t']: 0, loc['rho']: rho, loc['cp']: cp, loc['thcond']: k, loc['xi_T']: 1}
        A1, B1 = A.subs(sub), B.subs(sub)
        # integration by parts on test functions with N(0) = N(L) = 0 (Dirichlet ends): the weak
        # integral int (A1 N + B1 N' T') must equal int N * strong, strong = A1 - d/dx(B1 T')
        Tt = sp.Symbol('Tt_')
        strong = A1.subs(loc['tem_t'], sp.Derivative(Tf, t)) - B1 * sp.Derivative(Tf, (x, 2))
        Ttest = sp.sin(2 * x) * sp.exp(-t) + x**3 * t
        Ntest = x * (Lc - x) * (1 + x)
        weak_int = sp.integrate(A1.subs(loc['tem_t'], sp.diff(Ttest, t)) * Ntest
                                + B1 * sp.diff(Ntest, x) * sp.diff(Ttest, x), (x, 0, Lc))
        strong_int = sp.integrate(C.plug(strong, 'T', sp.Lambda((x, t), Ttest)) * Ntest, (x, 0, Lc))
        ibp_ok = sp.simplify(weak_int - strong_int) == 0
        r = sp.expand(r15())
        dTt, dTxx = sp.Derivative(Tf, t), sp.Derivative(Tf, (x, 2))
        storage = sp.simplify(r.coeff(dTt) - strong.coeff(dTt)) == 0
        conduction = sp.simplify(r.coeff(dTxx) - sp.expand(strong).coeff(dTxx)) == 0
        rest = sp.simplify(r - sp.expand(strong))
        only_source = sp.simplify(rest + q) == 0
        ok = ibp_ok and storage and conduction and only_source and mixtures
        return ok, {'code line': line, 'code weak heat (R[a][1])': body,
                    'code strong 1D (phi uniform, xi_T=1, rw=1)': str(strong),
                    'integration by parts checked on test N, T': ibp_ok,
                    'k, c_p, rho are phi-mixtures (constant for uniform phi)': mixtures,
                    'toy S15 residual': str(r), 'storage matches': storage, 'conduction matches': conduction,
                    'toy - code': str(rest), 'difference is exactly -q (the source)': only_source,
                    'assumed': 'xi_T = 1, phi uniform and constant, axisym off, decouple flag irrelevant (phi_t = 0)'}
    rep.check('reduction_vs_Residual_A1',
              "with psi == 1 (no grains, uniform phi) S15 equals today's Residual_A1 heat block term by term, "
              "plus only the source q", reduction_vs_code)

    def s40():
        mb = re.search(r'Rb\[a\]\[1\]\s*=\s*([^;]+);', src)
        heat_bdry = mb.group(1).strip() if mb else None
        heat_zero = heat_bdry is not None and float(sp.sympify(heat_bdry)) == 0.0
        mw = re.search(r'if \(kind == FACE_WALL\) \{(.*?)\} else', src, re.S)
        start = src[:mw.start()].count('\n') + 1 if mw else None
        end = start + mw.group(0).count('\n') - 1 if mw else None
        wall_uses_phi = bool(mw and 'phi' in mw.group(1))
        cited = re.search(r'lines\s*(\d+)-(\d+)', L.op('S40'))
        lines_ok = bool(cited) and (int(cited.group(1)), int(cited.group(2))) == (start, end)
        says_none = 'none' in L.stmt('S40')
        s8_no_phi = '\\phi' in L.stmt('S8') and re.search(r'[Nn]o ', L.stmt('S8')) is not None
        ok = heat_zero and wall_uses_phi and says_none and s8_no_phi and lines_ok
        return ok, {'Rb[a][1] in code': heat_bdry, 'wall term acts only through phi': wall_uses_phi,
                    'FACE_WALL block lines': (start, end), 'operation cites lines': cited.groups() if cited else None,
                    'S40 states none': says_none, 'S8 excludes phi': s8_no_phi}
    rep.check('S40', L.stmt('S40'), s40)

    def s46():
        r = r15()
        names = {str(s) for s in r.free_symbols}
        present = [nm for nm in GRAIN_ICE if nm in names]
        return not present and 'does not apply' in L.stmt('S46'), {'interface symbols in S15': present}
    rep.check('S46', L.stmt('S46'), s46)

    def s47():
        r = r15()
        eps, epss = C.S('eps'), C.S('eps_s')
        lim = sp.simplify(r.subs({eps: 0, epss: 0}) - r) == 0
        first = [sp.simplify(sp.diff(r, vv).subs({eps: 0, epss: 0})) for vv in (eps, epss)]
        ok = lim and all(f == 0 for f in first) and '\\equiv 0' in L.stmt('S47')
        return ok, {'limit leaves S15 unchanged': lim, 'dS15/d(eps, eps_s) at 0': [str(f) for f in first]}
    rep.check('S47', L.stmt('S47'), s47)

    rep.finish()


if __name__ == '__main__':
    main()
