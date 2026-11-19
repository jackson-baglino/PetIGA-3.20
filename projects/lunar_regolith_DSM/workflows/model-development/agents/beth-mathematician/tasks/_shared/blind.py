"""Shared machinery for the blind-derivation graders (T04 onward). Never shown to the agent.

In a blind task the agent is told the problem, not the equation, and names its own symbols, each with units
and a typical SI value. These helpers turn its answer into something checks can run on:
- parse:     only the work order's names, the agent's declared names, Derivative, numbers and arithmetic;
- dims:      the units of every term (pint), failing on a sum that mixes units;
- concrete:  each declared parameter becomes its typical value times a random factor in [1/2, 2], each
             declared function its typical value times a smooth positive variation of its arguments, and each
             field a given test field; derivatives are then evaluated exactly;
- helpers for the checks: solve for a rate, average over a period.
"""
import json
import random
import re

import pint
import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.parsing.sympy_parser import parse_expr

UREG = pint.UnitRegistry()
MAX_LEN = 6000
RESERVED = {"Derivative", "Integral", "Subs", "Function", "Symbol", "Lambda", "lambda", "import"}

x, y, z, t = sp.symbols("x y z t", real=True)
COORDS = (x, y, z)
T, phi, rho_v = sp.Function("T"), sp.Function("phi"), sp.Function("rho_v")
KNOWN = {"x": x, "y": y, "z": z, "t": t, "T": T, "phi": phi, "rho_v": rho_v}
# Typical size of each field, used to scale a declared function's variation with an argument that holds it.
T_SCALE = 200                      # K, the temperature the typical values refer to
RHOV_SCALE = sp.Rational(1, 10**6)  # kg m^-3, saturated vapor near 200 K is about 1.6e-6
SCALES = ((T, T_SCALE), (rho_v, RHOV_SCALE))


class Reject(Exception):
    pass


def quantity(text):
    try:
        q = UREG(str(text).replace("^", "**"))
    except Exception as e:
        raise Reject(f"units {text!r} do not parse: {e}") from None
    return q if isinstance(q, UREG.Quantity) else UREG.Quantity(q, "dimensionless")


def _task(task_dir):
    return json.loads((task_dir / "task.json").read_text())


def fixed_symbols(task_dir):
    """The names the work order defines, as SymPy objects."""
    out = {}
    for s in _task(task_dir)["symbols"]:
        name = s["name"].split("(")[0].strip()
        out[name] = KNOWN.get(name) or (sp.Function(name) if "(" in s["name"] else sp.Symbol(name, real=True))
    return out


def fixed_units(task_dir):
    return {s["name"].split("(")[0].strip(): quantity(s["units"]) for s in _task(task_dir)["symbols"]}


def declared(symbols, fixed):
    """{name: (sympy object, units, typical value)} for the agent's own symbols."""
    out = {}
    for s in symbols or []:
        m = re.fullmatch(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*(\(.*\))?\s*", str(s.get("name", "")))
        if not m:
            raise Reject(f"symbol name {s.get('name')!r} is not a name or a name with arguments")
        name, is_func = m.group(1), bool(m.group(2))
        if name in RESERVED or name.startswith("_"):
            raise Reject(f"symbol name {name!r} is reserved")
        if name in fixed:
            continue                          # the work order defines it
        try:
            typical = sp.Rational(str(s.get("typical")))
        except (TypeError, ValueError):
            raise Reject(f"symbol {name!r} has no numeric typical value") from None
        # A typical value of 0 is refused only for a symbol the result uses (concrete); the agent may declare
        # intermediate symbols with no natural size (a residual R, an enthalpy with an arbitrary zero).
        # Found in runs 20261118-161936 (h_v) and 20261118-162313 (R).
        out[name] = (sp.Function(name) if is_func else sp.Symbol(name, real=True), quantity(s.get("units")), typical)
    return out


def parse(text, decl, fixed):
    if not isinstance(text, str) or not text.strip():
        raise Reject("empty result_expr")
    if len(text) > MAX_LEN:
        raise Reject(f"result_expr longer than {MAX_LEN} characters")
    if not re.fullmatch(r"[A-Za-z0-9_+\-*/(),. \t\n]*", text):
        raise Reject("characters outside names, numbers, + - * / ** ( ) , and spaces")
    if re.search(r"(?<![0-9])\.|\.(?![0-9])", text):
        raise Reject("'.' is allowed only inside a decimal number")
    names = set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", text))
    unknown = names - set(fixed) - set(decl) - {"Derivative"}
    if unknown:
        raise Reject(f"names neither in the work order nor declared in symbols: {sorted(unknown)}")
    local = dict(fixed, Derivative=sp.Derivative, **{n: d[0] for n, d in decl.items()})
    try:
        return parse_expr(text, local_dict=local, evaluate=True)
    except Exception as e:      # a parse error is the agent's fault, not the grader's
        raise Reject(f"does not parse: {type(e).__name__}: {e}") from None


def dims(e, units):
    """The units of e as a pint quantity; raises Reject when a sum mixes units."""
    if e.is_Number:
        return UREG.Quantity(1, "dimensionless")
    if isinstance(e, sp.Symbol):
        return units[e.name]
    if isinstance(e, AppliedUndef):
        return units[e.func.__name__]
    if isinstance(e, sp.Derivative):
        q = dims(e.expr, units)
        for v, n in e.variable_count:
            q = q / dims(v, units) ** int(n)
        return q
    if isinstance(e, sp.Add):
        qs = [(a, dims(a, units)) for a in e.args]
        first = qs[0][1].dimensionality
        bad = [(str(a), str(q.dimensionality)) for a, q in qs if q.dimensionality != first]
        if bad:
            raise Reject(f"terms with different units: {str(qs[0][0])[:80]} is {first}; {bad[:3]}")
        return qs[0][1]
    if isinstance(e, sp.Mul):
        q = UREG.Quantity(1, "dimensionless")
        for a in e.args:
            q = q * dims(a, units)
        return q
    if isinstance(e, sp.Pow) and e.exp.is_Number:
        return dims(e.base, units) ** float(e.exp)
    raise Reject(f"cannot take the units of {type(e).__name__}: {str(e)[:80]}")


def _scale(arg):
    for f, s in SCALES:
        if arg.has(f):
            return s
    return 1


def concrete(expr, decl, fields, seed):
    """expr with each declared parameter replaced by its typical value times a random factor in [1/2, 2], each
    declared function by its typical value times a random smooth positive variation of its arguments (each
    argument scaled by the typical size of the field it holds), and each field f by fields[f] (a Lambda of
    x, y, z, t); derivatives evaluated. Typical values keep the sizes physical, so a sign that depends on
    sizes (rho_ice - rho_v) is judged at physical sizes."""
    rng = random.Random(seed)
    e = expr
    used = expr.free_symbols | {a.func for a in expr.atoms(AppliedUndef)}
    for name, (obj, _, typical) in decl.items():
        factor = sp.Rational(rng.randint(2, 8), 4)
        if obj not in used:
            continue
        if typical == 0:
            raise Reject(f"symbol {name!r} is used in the result but has a typical value of 0")
        if isinstance(obj, sp.Symbol):
            e = e.subs(obj, typical * factor)
        else:
            coef = [sp.Rational(rng.randint(1, 10), 50) for _ in range(8)]
            def poly(*args, c=coef, base=typical * factor):
                return base * (1 + sum(c[i] * (a / _scale(a))**2 for i, a in enumerate(args[:8])))
            e = e.replace(lambda a, f=obj: isinstance(a, AppliedUndef) and a.func == f,
                          lambda a, p=poly: p(*a.args))
    for f, body in fields.items():
        e = e.replace(f, body)
    return e.doit()


def lam(body):
    return sp.Lambda((x, y, z, t), body)


def solve_rate(e, a):
    sol = [s for s in sp.solve(sp.Eq(e, 0), a) if s.is_real]
    if len(sol) != 1:
        raise Reject(f"expected one real rate from R = 0, found {sol}")
    return sol[0]


def periodic_mean(e, n=48):
    """(mean over x in [0, 1), max |value|) of a concrete expression, at y = z = t = 0. The trapezoid rule is
    exact for trigonometric polynomials of low degree, so this is the period average to round-off."""
    f = sp.lambdify(x, e.subs({y: 0, z: 0, t: 0}), "math")
    vals = [float(f(j / n)) for j in range(n)]
    return sum(vals) / n, max(abs(v) for v in vals)


def at_origin(e):
    return sp.simplify(e.subs({x: 0, y: 0, z: 0, t: 0}))
