"""M3 check: does each derived step need an assumption that the ledger does not list?

For each operation I state, independently, the hypotheses it needs (table REQUIRED
below, written from the operation names, not from the ledger's from lists), then test
whether the ledger supplies each one: in the step's from list or in the closure of
its premises, and, where the hypothesis is a property, in the text of an assumption.
Where useful, a counterexample shows the hypothesis is really needed.
"""
import os
import re
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = ['S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S20', 'S21', 'S22',
          'S23', 'S24', 'S25', 'S26', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33', 'S34', 'S35',
          'S36', 'S37', 'S38', 'S39', 'S40', 'S41', 'S42', 'S43', 'S44', 'S45', 'S46', 'S47', 'S48',
          'S49', 'S50', 'S51', 'S52', 'S53', 'S54', 'S55', 'S56', 'S57']

# hypotheses each named operation needs (my reading of the operation, step ids of the ledger)
REQUIRED = {
    'S9': {'S2': 'balance', 'S3': 'e in terms of T', 'S6': 'T_t continuous up to the ends'},
    'S10': {'S4': 'j in terms of T', 'S6': 'T_xx continuous up to the ends'},
    'S11': {'S2': 'balance', 'S9': 'Leibniz', 'S10': 'FTC'},
    'S12': {'S11': 'integral form', 'S6': 'continuity of e_t, j_x', 'S3': 'e_t', 'S4': 'j_x',
            'S5': 'q continuous (constant) for the localization lemma'},
    'S13': {'S12': 'local balance', 'S3': 'caloric law', 'S1': 'rho, c_p (and T_ref) constant'},
    'S14': {'S4': 'Fourier', 'S1': 'k constant'},
    'S15': {'S13': '', 'S14': ''},
    'S16': {'S1': 'rho c_p > 0'},
    'S17': {'S15': '', 'S16': ''},
    'S18': {'S8': 'no grains', 'S15': ''},
    'S20': {'S1': 'L, k > 0', 'S5': 'q constant', 'S16': 'kappa'},
    'S21': {'S20': 'change of variables', 'S6': 'T differentiable (chain rule needs it)'},
    'S22': {'S20': 'change of variables', 'S6': 'T twice differentiable in x (chain rule)'},
    'S23': {'S17': '', 'S21': '', 'S22': ''},
    'S24': {'S23': '', 'S16': 'kappa > 0', 'S20': 'DeltaT_q != 0'},
    'S25': {'S24': '', 'S16': '', 'S20': ''},
    'S26': {'S7': 'data', 'S20': 'map'},
    'S27': {'S7': 'T_0 > 0', 'S16': ''},
    'S28': {'S17': '', 'S27': '', 'S6': 'chain rule'},
    'S29': {'S28': '', 'S16': '', 'S27': '', 'S7': ''},
    'S30': {'S20': '', 'S27': ''},
    'S31': {'S19': 'reference values', 'S16': '', 'S20': '', 'S27': ''},
    'S33': {'S17': 'linear equation', 'S7': 'shared data'},
    'S34': {'S33': '', 'S6': 'w_t continuous up to the ends'},
    'S35': {'S34': '', 'S33': 'w = 0 at the ends', 'S6': 'w_x, w_xx up to the ends', 'S16': 'kappa > 0'},
    'S36': {'S35': '', 'S34': '', 'S33': 'w(x,0) = 0', 'S6': 'continuity at t = 0'},
    'S37': {'S36': '', 'S17': '', 'S7': ''},
    'S48': {'S6': 'continuity on the closed set', 'S20': 'theta defined (q != 0)'},
    'S49': {'S48': '', 'S26': 'theta = 0 on three sides', 'S6': 'C^2 in xi, T_t at interior/top points'},
    'S50': {'S49': '', 'S25': 'source > 0', 'S48': ''},
    'S38': {'S36': 'uniqueness', 'S20': '', 'S1': 'k > 0', 'S50': 'theta >= 0', 'S25': '', 'S26': ''},
    'S39': {'S4': '', 'S7': 'Dirichlet ends', 'S6': 'T_x at the ends'},
    'S40': {'S8': ''},
    'S41': {'S25': '', 'S26': ''},
    'S42': {'S41': ''},
    'S43': {'S42': '', 'S41': ''},
    'S44': {'S43': '', 'S20': ''},
    'S45': {'S44': '', 'S39': ''},
    'S51': {'S25': '', 'S26': '', 'S41': '', 'S43': ''},
    'S52': {},
    'S53': {'S43': '', 'S52': ''},
    'S54': {'S51': '', 'S53': '', 'S38': 'uniqueness of theta', 'S6': 'class'},
    'S55': {'S54': ''},
    'S56': {'S55': ''},
    'S57': {'S54': '', 'S55': ''},
    'S32': {'S25': '', 'S26': '', 'S6': '', 'S56': '', 'S57': '', 'S29': '', 'S31': ''},
    'S46': {'S8': ''},
    'S47': {'S8': '', 'S15': ''},
}


def cited_in_operation(op, sid):
    op_wo = re.sub(r'as in S\d+(?:\s*-\s*S\d+)?', '', op)   # "as in Sx" names a pattern, not a premise
    ids = set()
    for a, b in re.findall(r'\bS(\d+)(?:\s*-\s*S(\d+))?', op_wo):
        if b:
            ids |= {f'S{i}' for i in range(int(a), int(b) + 1)}
        else:
            ids.add(f'S{a}')
    return ids - {sid}


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_assumptions.py', L)
    x, t = C.S('x'), C.S('t')

    def from_lists():
        bad = {}
        for sid in L.derived_ids():
            st = L.by_id[sid]
            missing = sorted(cited_in_operation(st.get('operation', ''), sid) - set(st.get('from', [])),
                             key=lambda v: int(v[1:]))
            if missing:
                bad[sid] = {'operation cites': missing, 'from': st.get('from')}
        as_in = {sid: re.findall(r'as in S\d+(?:\s*-\s*S\d+)?', L.op(sid)) for sid in L.derived_ids()
                 if re.search(r'as in S\d+', L.op(sid))}
        return not bad, {'missing from from-lists': bad,
                         'informative: "as in Sx" pattern citations (not required in from)': as_in}
    rep.check('from_lists', "each step named as a premise by an operation is listed in that step's from", from_lists)

    def required():
        missing = {}
        for sid in L.derived_ids():
            need = REQUIRED.get(sid)
            if need is None:
                missing[sid] = 'no entry in my REQUIRED table (new step?)'
                continue
            have = set(L.by_id[sid].get('from', [])) | L.closure(sid)
            lack = {p: why for p, why in need.items() if p not in have}
            if lack:
                missing[sid] = lack
        return not missing, {'hypotheses not in the step\'s premises or their closure': missing}
    rep.check('required_premises', 'every hypothesis an operation needs is among its premises (directly or through '
              'the steps it cites)', required)

    def tref_constant():
        T = C.Fn('T')(x, t)
        rho, cp = C.S('rho'), C.S('c_p')
        Tr = sp.Function('Tref')(t)
        extra = sp.simplify(sp.diff(rho * cp * (T - Tr), t) - rho * cp * sp.diff(T, t))
        # Only an assumption STATEMENT counts (the symbols table is not an assumption), and the
        # assumption that states it must be among S13's premises (directly or by closure).
        have = set(L.by_id['S13'].get('from', [])) | L.closure('S13')
        assumptions = [s['id'] for s in L.steps if s.get('type') == 'assumption']
        stating = [sid for sid in assumptions
                   if re.search(r'T_\{?ref\}?[^;]*constant[^;]*\bt\b|constant[^;]*T_\{?ref\}?', L.stmt(sid))]
        declared = any(sid in have for sid in stating)
        tref_sym = [s['meaning'] for s in L.data['symbols'] if s['symbol'] in ('T_{ref}', 'T_ref')]
        op_claims = 'T_{ref}' in L.op('S13') and '(S1)' in L.op('S13')
        return declared, {'extra term in S13 if T_ref depends on t': str(extra),
                          'assumptions stating T_ref constant in t': stating,
                          'stating assumption among S13 premises': declared,
                          'S13 operation attributes T_ref constancy to S1': op_claims,
                          'T_ref in symbols table (informative only)': tref_sym,
                          'finding': None if declared else
                          'no assumption among S13 premises states that T_ref is constant in t'}
    rep.check('tref_constant_for_S13', 'S13 needs T_ref constant in time, and an assumption states it', tref_constant)

    def regularity_at_ends():
        s6 = L.stmt('S6').replace(' ', '')
        ok = 'C^{2,1}\\big([0,L]\\times(0,\\infty)\\big)' in s6
        return ok, {'S6 has C^{2,1} on [0,L] x (0,inf)': ok,
                    'needed by': ['S9-S11 with a = 0 or b = L', 'S34', 'S35', 'S39', 'S45', 'S32 (theta_tau at xi = 0)']}
    rep.check('regularity_at_ends', 'S6 makes T, T_x, T_xx, T_t continuous up to x = 0 and x = L for t > 0',
              regularity_at_ends)

    def continuity_t0():
        s6 = L.stmt('S6').replace(' ', '')
        ok = 'C\\big([0,L]\\times[0,\\infty)\\big)' in s6
        return ok, {'S6 has C on [0,L] x [0,inf)': ok, 'needed by': ['S36 (E -> 0 as t -> 0+)', 'S48 (compact set)']}
    rep.check('continuity_at_t0', 'S6 makes T continuous on [0,L] x [0,inf), including t = 0', continuity_t0)

    def theta_sign():
        p50 = [p for p in C.pieces(L.stmt('S50')) if p.ok]
        stated = bool(p50) and p50[0].exprs == [C.field('theta'), 0] and p50[0].rels == ['>=']
        cited = 'S50' in L.by_id['S38'].get('from', [])
        return stated and cited, {'S50 states theta >= 0': stated, 'S38 cites S50': cited}
    rep.check('theta_nonnegative', "the ledger supplies theta >= 0 and S38's sign clause uses it", theta_sign)

    def k_constant():
        kx = sp.Function('kx')(x)
        T = C.Fn('T')(x, t)
        extra = sp.simplify(sp.diff(-kx * sp.diff(T, x), x) - (-kx * sp.diff(T, x, 2)))
        listed = 'S1' in L.by_id['S14']['from'] and 'constant' in L.stmt('S1') and 'k' in L.stmt('S1')
        return extra != 0 and listed, {'extra term if k = k(x)': str(extra), 'S1 (k constant) in S14 from': listed}
    rep.check('k_constant_for_S14', 'S14 needs k constant, and cites S1 for it', k_constant)

    def q_nonzero():
        users = [sid for sid in L.derived_ids() if '\\theta' in L.stmt(sid) and "\\theta'" not in L.stmt(sid)]
        out = {}
        for sid in users:
            out[sid] = sid == 'S20' or 'S20' in (L.closure(sid) | set(L.by_id[sid].get('from', [])))
        s20_says = '\\ne 0' in L.stmt('S20')
        return all(out.values()) and s20_says, {'uses theta and rests on S20 (q != 0)': out,
                                                'S20 restricts to q != 0': s20_says}
    rep.check('q_nonzero', 'every step that uses theta rests on S20, which is restricted to q != 0', q_nonzero)

    def source_positive_for_S50():
        r25 = C.piece(L.stmt('S25'), 0).residual()
        th = C.field('theta')
        c = -r25.subs({sp.Derivative(th, C.S('tau')): 0, sp.Derivative(th, (C.S('xi'), 2)): 0})
        return bool(sp.sympify(c).is_positive), {
            'S25 source': str(c),
            'why': 'S50 contradicts theta_tau - theta_xixi <= 0 at a negative interior minimum only if the source > 0'}
    rep.check('source_positive_for_S50', 'the minimum principle S50 has the strictly positive source it relies on',
              source_positive_for_S50)

    rep.finish()


if __name__ == '__main__':
    main()
