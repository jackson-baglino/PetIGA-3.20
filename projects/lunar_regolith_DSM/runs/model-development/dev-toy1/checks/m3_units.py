"""M3 check: units of every term (Pint).

Units of each symbol come from the ledger's own symbols table (--ledger). For every
relation parsed from every step (assumptions and derived), each additive term and
each member of an equality or inequality chain must have the same dimension.
Derivatives divide by the units of the variable, integrals multiply by it, and the
arguments of sin/exp/log must be dimensionless. A literal 0 matches any unit.
S19 and S31 carry explicit units (\\mathrm{...}); those are compared with the table.
"""
import os
import re
import sys

import pint
import sympy as sp
from sympy.core.function import AppliedUndef

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

ureg = pint.UnitRegistry()
SKIP_GENERIC = {'S19', 'S31'}   # numbers with explicit units: dedicated check below


class Mismatch(Exception):
    pass


def _split_top(s):
    parts, depth, cur = [], 0, ''
    for ch in s:
        if ch in '({[':
            depth += 1
        elif ch in ')}]':
            depth -= 1
        if ch == ',' and depth == 0:
            parts.append(cur)
            cur = ''
        else:
            cur += ch
    parts.append(cur)
    return [p.strip() for p in parts if p.strip()]


def _name_of(raw):
    if re.fullmatch(r'[A-Za-z]{2,}(_[A-Za-z0-9]+)?', raw):
        return raw
    toks = C.tokenize(raw)
    if not toks or toks[0][0] != 'name':
        return None
    name = toks[0][1]
    if name == 'theta' and len(toks) > 1 and toks[1] == ('op', "'"):
        name = 'thetap'
    return name


def _unit(u):
    u = u.strip()
    if u in ('1', ''):
        return ureg.dimensionless
    return ureg('*'.join(u.replace('^', '**').split())).units


def unit_table(L):
    tab = {}
    for s in L.data['symbols']:
        parts = _split_top(s['symbol'])
        us = [u.strip() for u in s['units'].split(',')]
        if len(us) == 1:
            us = us * len(parts)
        for p, u in zip(parts, us):
            n = _name_of(p)
            if n:
                tab[n] = _unit(u)
    return tab


def units_of(e, tab):
    """Return pint units of expression e, or None for a literal zero (matches anything)."""
    if e == 0:
        return None
    if e.is_number:
        return ureg.dimensionless
    if isinstance(e, sp.Symbol):
        if e.name not in tab:
            raise KeyError(f'symbol {e.name} has no units in the ledger symbols table')
        return tab[e.name]
    if isinstance(e, AppliedUndef):
        n = e.func.__name__
        if n.startswith('Max_') or n.startswith('Min_'):
            return units_of(e.args[0], tab)
        if n not in tab:
            raise KeyError(f'field {n} has no units in the ledger symbols table')
        return tab[n]
    if isinstance(e, sp.Derivative):
        u = units_of(e.expr, tab)
        for v, n in e.variable_count:
            u = u / units_of(v, tab) ** int(n)
        return u
    if isinstance(e, sp.Integral):
        u = units_of(e.function, tab)
        for lim in e.limits:
            u = u * units_of(lim[0], tab)
        return u
    if isinstance(e, sp.Subs):
        return units_of(e.expr, tab)
    if isinstance(e, sp.Sum):
        return units_of(e.function, tab)
    if isinstance(e, sp.Abs):
        return units_of(e.args[0], tab)
    if isinstance(e, sp.Add):
        us = [units_of(a, tab) for a in e.args]
        us = [u for u in us if u is not None]
        for u in us[1:]:
            if (1 * u).dimensionality != (1 * us[0]).dimensionality:
                raise Mismatch(f'terms of {e} have units {[str(x_) for x_ in us]}')
        return us[0] if us else None
    if isinstance(e, sp.Mul):
        u = ureg.dimensionless
        for a in e.args:
            ua = units_of(a, tab)
            if ua is None:
                return None
            u = u * ua
        return u
    if isinstance(e, sp.Pow):
        ub = units_of(e.base, tab)
        if (1 * ub).dimensionless:
            return ureg.dimensionless
        if not e.exp.is_number:
            raise Mismatch(f'non-numeric exponent on a dimensional base: {e}')
        return ub ** float(e.exp)
    if isinstance(e, (sp.sin, sp.cos, sp.tan, sp.exp, sp.log)):
        ua = units_of(e.args[0], tab)
        if ua is not None and not (1 * ua).dimensionless:
            raise Mismatch(f'argument of {e.func} is not dimensionless: {ua}')
        return ureg.dimensionless
    raise Mismatch(f'cannot assign units to {e} ({type(e).__name__})')


def check_piece(p, tab):
    us = [units_of(e, tab) for e in p.exprs]
    known = [u for u in us if u is not None]
    for u in known[1:]:
        if (1 * u).dimensionality != (1 * known[0]).dimensionality:
            raise Mismatch(f'members have units {[str(x_) for x_ in us]}')
    return str((1 * known[0]).to_base_units().units) if known else 'any (all zero)'


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    tab = unit_table(L)
    rep = C.Reporter('m3_units.py', L, {'unit_table': {kk: str(v) for kk, v in sorted(tab.items())}})
    for sid in L.order:
        if sid in SKIP_GENERIC:
            continue
        ps = [p for p in C.pieces(L.stmt(sid)) if p.ok]
        if not ps:
            continue

        def one(ps=ps):
            rows = []
            for p in ps:
                rows.append((p.text, check_piece(p, tab)))
            return True, rows
        rep.check(f'units_{sid}', f'every term of {sid} has the same units', one)

    def explicit_units():
        res = {}
        for sid in sorted(SKIP_GENERIC):
            st = L.stmt(sid)
            for m in re.finditer(r'(\\Delta\s*T_q|\\[A-Za-z]+|[A-Za-z](?:_\{?[A-Za-z0-9]+\}?)?)\s*=\s*[0-9.]+\s*'
                                 r'(\\,\\mathrm\{)?', st):
                name = _name_of(m.group(1).strip())
                if name is None:
                    continue
                if m.group(2):
                    content, _ = C._balanced(st, m.end() - 1)
                    u = content.replace('\\,', ' ').replace('{', '').replace('}', '')
                    u = re.sub(r'\^(-?\d+)', r'^\1', u)
                    u = re.sub(r'([A-Za-z])(?=[A-Z])', r'\1 ', u)
                    stated = _unit(u)
                else:
                    stated = ureg.dimensionless
                if name not in tab:
                    res[f'{sid}:{name}'] = 'not in table'
                    continue
                match = (1 * stated).dimensionality == (1 * tab[name]).dimensionality
                res[f'{sid}:{name}'] = {'stated': str(stated), 'table': str(tab[name]), 'match': match}
        ok = bool(res) and all(isinstance(v, dict) and v['match'] for v in res.values())
        return ok, res
    rep.check('units_S19_S31', 'the units stated in S19 and S31 match the symbols table', explicit_units)
    rep.finish()


COVERS = ['S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S20', 'S21', 'S22', 'S23',
          'S24', 'S25', 'S26', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33', 'S34', 'S35', 'S36', 'S37',
          'S38', 'S39', 'S41', 'S42', 'S43', 'S44', 'S45', 'S48', 'S49', 'S50', 'S51', 'S52', 'S53', 'S54',
          'S55', 'S56', 'S57']

if __name__ == '__main__':
    main()
