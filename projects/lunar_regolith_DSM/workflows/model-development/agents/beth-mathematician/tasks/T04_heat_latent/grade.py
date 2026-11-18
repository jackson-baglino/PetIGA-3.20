#!/usr/bin/env python3
"""Grader for T04: the energy equation for T with the latent heat of sublimation, derived blind.

Usage: grade.py <result.json>

The agent never sees this file. It was told the problem, not the equation, and it names its own symbols, so
there is no single answer key. The grade has two parts.

1. Checks that need no key (all must pass), on the agent's own declared symbols and units.
   Each declared symbol carries a typical SI value; the draws stay within a factor of 2 of it, so a sign
   that depends on sizes (rho_ice - rho_v) is judged at physical sizes, near T = 200 K.
   - names:        R uses only T, phi, x, y, z, t, Derivative and the symbols the agent declared;
   - units:        every term of R has the same units (pint, from the declared units);
   - equilibrium:  uniform T and a static phi give R = 0 (no source out of nothing);
   - latent sign:  uniform T(t) and uniform phi(t) with phi_t > 0 (deposition) give T_t > 0 (warming);
   - conduction:   a hot spot (T peaked at a point, phi uniform and static) cools: T_t < 0 at the peak;
   - conservation: with phi and T static and periodic in x, R integrates to zero over a period, so
                   conduction is the divergence of a flux and moves energy without making it.
2. A comparison with the current model (src/assembly.c:173-175, with the time scaling xi_T = 1 and the
   phase-change switch on): R_ref = rho(phi) cp(phi) T_t - div(k(phi) grad T) - rho_i L phi_t. The agent's
   terms are grouped by role (storage, conduction, phase change, other) beside the reference's, with the
   effects the agent left out, for the user to judge. Terms in "other" are listed, not failed: a term the
   current model lacks may be right.

The checks themselves are proven by key_check, with no LLM in it: R_ref passes all of them, and each of five
corrupted references fails at least one.

Prints one JSON line {"check", "pass", "details"}; exits 0 on pass, 1 on fail, 2 on an unreadable file.
"""
import hashlib
import json
import math
import random
import re
import sys
from pathlib import Path

import pint
import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.parsing.sympy_parser import parse_expr

CHECK = "T04_heat_latent"
HERE = Path(__file__).resolve().parent
UREG = pint.UnitRegistry()
MAX_LEN = 4000
RESERVED = {"Derivative", "Integral", "Subs", "Function", "Symbol", "Lambda", "lambda", "import"}

x, y, z, t = sp.symbols("x y z t", real=True)
COORDS = (x, y, z)
T, phi = sp.Function("T"), sp.Function("phi")
FIXED = {"x": x, "y": y, "z": z, "t": t, "T": T, "phi": phi}


class Reject(Exception):
    pass


def quantity(text):
    try:
        q = UREG(str(text).replace("^", "**"))
    except Exception as e:
        raise Reject(f"units {text!r} do not parse: {e}") from None
    return q if isinstance(q, UREG.Quantity) else UREG.Quantity(q, "dimensionless")


def fixed_units():
    task = json.loads((HERE / "task.json").read_text())
    return {s["name"].split("(")[0].strip(): quantity(s["units"]) for s in task["symbols"]}


def declared(symbols):
    """{name: (sympy object, units, typical value)} for the agent's own symbols."""
    out = {}
    for s in symbols or []:
        m = re.fullmatch(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*(\(.*\))?\s*", str(s.get("name", "")))
        if not m:
            raise Reject(f"symbol name {s.get('name')!r} is not a name or a name with arguments")
        name, is_func = m.group(1), bool(m.group(2))
        if name in RESERVED or name.startswith("_"):
            raise Reject(f"symbol name {name!r} is reserved")
        if name in FIXED:
            continue                          # the work order defines it
        try:
            typical = sp.Rational(str(s.get("typical")))
        except (TypeError, ValueError):
            raise Reject(f"symbol {name!r} has no numeric typical value") from None
        if typical == 0:
            raise Reject(f"symbol {name!r} has a typical value of 0")
        out[name] = (sp.Function(name) if is_func else sp.Symbol(name, real=True), quantity(s.get("units")), typical)
    return out


def parse(text, decl):
    if not isinstance(text, str) or not text.strip():
        raise Reject("empty result_expr")
    if len(text) > MAX_LEN:
        raise Reject(f"result_expr longer than {MAX_LEN} characters")
    if not re.fullmatch(r"[A-Za-z0-9_+\-*/(),. \t\n]*", text):
        raise Reject("characters outside names, numbers, + - * / ** ( ) , and spaces")
    if re.search(r"(?<![0-9])\.|\.(?![0-9])", text):
        raise Reject("'.' is allowed only inside a decimal number")
    names = set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", text))
    unknown = names - set(FIXED) - set(decl) - {"Derivative"}
    if unknown:
        raise Reject(f"names neither in the work order nor declared in symbols: {sorted(unknown)}")
    local = dict(FIXED, Derivative=sp.Derivative, **{n: d[0] for n, d in decl.items()})
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


T_SCALE = 200        # K: the fields below are near 200 K, the temperature the typical values refer to


def concrete(expr, decl, Tlam, philam, seed):
    """expr with each declared parameter replaced by its typical value times a random factor in [1/2, 2], each
    declared function by its typical value times a random smooth positive variation of its arguments
    (scaled by 200 K for an argument holding T, by 1 otherwise), and T, phi by the given fields; derivatives
    evaluated. Typical values keep the sizes physical, so a sign that depends on sizes (rho_ice - rho_v)
    is judged at physical sizes."""
    rng = random.Random(seed)
    e = expr
    for obj, _, typical in decl.values():
        factor = sp.Rational(rng.randint(2, 8), 4)
        if isinstance(obj, sp.Symbol):
            e = e.subs(obj, typical * factor)
        else:
            coef = [sp.Rational(rng.randint(1, 10), 50) for _ in range(8)]
            def poly(*args, c=coef, base=typical * factor):
                return base * (1 + sum(c[i] * (a / (T_SCALE if a.has(T) else 1))**2
                                       for i, a in enumerate(args[:8])))
            e = e.replace(lambda a, f=obj: isinstance(a, AppliedUndef) and a.func == f,
                          lambda a, p=poly: p(*a.args))
    return e.replace(T, Tlam).replace(phi, philam).doit()


def lam(body):
    return sp.Lambda((x, y, z, t), body)


def check_equilibrium(expr, decl):
    T0 = sp.Integer(T_SCALE)
    field = lam(sp.Rational(1, 2) + (x - y + 2 * z) / 10)
    for seed in range(3):
        e = concrete(expr, decl, lam(T0 + 0 * x), field, seed)
        for pt in ((0, 0, 0), (sp.Rational(1, 3), sp.Rational(-1, 5), sp.Rational(1, 7))):
            v = sp.simplify(e.subs(dict(zip(COORDS, pt))).subs(t, 0))
            if v != 0:
                return False, f"uniform T and a static phi give R = {v}, not 0: a source out of nothing"
    return True, "R = 0 for uniform T and a static phi (3 draws, 2 points)"


def solve_rate(e, a):
    sol = [s for s in sp.solve(sp.Eq(e, 0), a) if s.is_real]
    if len(sol) != 1:
        raise Reject(f"expected one real T_t from R = 0, found {sol}")
    return sol[0]


def check_latent_sign(expr, decl):
    a = sp.Symbol("a")
    for seed in range(3):
        e = concrete(expr, decl, lam(T_SCALE + a * t), lam(sp.Rational(1, 2) + t), seed)
        rate = solve_rate(sp.simplify(e.subs({x: 0, y: 0, z: 0, t: 0})), a)
        if not rate > 0:
            return False, f"deposition (phi_t > 0) at uniform T gives T_t = {rate}: it should warm (draw {seed})"
    return True, "deposition warms and, by linearity, sublimation cools (3 draws)"


def check_conduction(expr, decl):
    a = sp.Symbol("a")
    for seed in range(3):
        e = concrete(expr, decl, lam(T_SCALE - (x**2 + y**2 + z**2) + a * t),
                     lam(sp.Rational(1, 2) + 0 * x), seed)
        rate = solve_rate(sp.simplify(e.subs({x: 0, y: 0, z: 0, t: 0})), a)
        if not rate < 0:
            return False, f"a hot spot gives T_t = {rate} at its peak: heat should flow out (draw {seed})"
    return True, "a hot spot cools at its peak (3 draws)"


def check_conservation(expr, decl):
    n = 48
    for seed in range(3):
        e = concrete(expr, decl, lam(T_SCALE + 10 * sp.sin(2 * sp.pi * x) + 5 * sp.cos(4 * sp.pi * x)),
                     lam(sp.Rational(1, 2) + sp.sin(2 * sp.pi * x) / 4), seed)
        f = sp.lambdify(x, e.subs({y: 0, z: 0, t: 0}), "math")
        vals = [float(f(j / n)) for j in range(n)]
        scale = max(abs(v) for v in vals)
        if scale == 0:
            return False, "R vanishes for a non-uniform T: no conduction term"
        mean = sum(vals) / n
        if abs(mean) > 1e-9 * scale or not math.isfinite(mean):
            return False, (f"with static, periodic phi and T, R averages {mean:.3g} over a period (max |R| "
                           f"{scale:.3g}): conduction is not the divergence of a flux, so it makes energy")
    return True, "conduction integrates to zero over a period (3 draws): it is the divergence of a flux"


def run_checks(expr, decl, units):
    results = []
    try:
        q = dims(expr, units)
        results.append(("units", True, f"every term has units {q.units:~}"))
    except Reject as e:
        results.append(("units", False, str(e)))
    for name, fn in (("equilibrium", check_equilibrium), ("latent sign", check_latent_sign),
                     ("conduction", check_conduction), ("conservation", check_conservation)):
        try:
            results.append((name, *fn(expr, decl)))
        except Reject as e:
            results.append((name, False, str(e)))
        except Exception as e:      # the agent's expression broke a check: report it, do not crash
            results.append((name, False, f"could not evaluate: {type(e).__name__}: {str(e)[:160]}"))
    return results


def roles(expr):
    """The additive terms of R, grouped by what they do."""
    groups = {"storage": [], "conduction": [], "phase change": [], "other": []}
    for term in sp.Add.make_args(sp.expand(expr)):
        kinds = set()
        for d in term.atoms(sp.Derivative):
            wrt = {v for v, _ in d.variable_count}
            if d.expr.has(T) and wrt & set(COORDS):
                kinds.add("conduction")
            elif d.expr.has(T) and t in wrt:
                kinds.add("storage")
            elif d.expr.has(phi) and t in wrt:
                kinds.add("phase change")
            else:
                kinds.add("other")
        groups[kinds.pop() if len(kinds) == 1 else "other"].append(term)
    return {k: [str(v) for v in vs] for k, vs in groups.items()}


# ---- the reference: the current model's heat row --------------------------------------------------------
REF_SYMBOLS = [
    {"name": "rho(phi)", "meaning": "mixture density", "units": "kg m^-3", "typical": 500},
    {"name": "cp(phi)", "meaning": "mixture specific heat", "units": "J kg^-1 K^-1", "typical": 1500},
    {"name": "k(phi)", "meaning": "mixture thermal conductivity", "units": "W m^-1 K^-1", "typical": 1},
    {"name": "rho_i", "meaning": "density of ice", "units": "kg m^-3", "typical": 917},
    {"name": "L", "meaning": "latent heat of sublimation", "units": "J kg^-1", "typical": 2.8e6},
]
_P, _T = "phi(x, y, z, t)", "T(x, y, z, t)"
REF_STORAGE = f"rho({_P})*cp({_P})*Derivative({_T}, t)"
REF_CONDUCTION = "(" + " + ".join(f"Derivative(k({_P})*Derivative({_T}, {c}), {c})" for c in "xyz") + ")"
REF_LATENT = f"rho_i*L*Derivative({_P}, t)"
REFERENCE = f"{REF_STORAGE} - {REF_CONDUCTION} - {REF_LATENT}"
CORRUPTED = {
    "latent sign flipped": f"{REF_STORAGE} - {REF_CONDUCTION} + {REF_LATENT}",
    "conduction sign flipped": f"{REF_STORAGE} + {REF_CONDUCTION} - {REF_LATENT}",
    "k outside the divergence": f"{REF_STORAGE} - k({_P})*(" + " + ".join(
        f"Derivative({_T}, ({c}, 2))" for c in "xyz") + f") - {REF_LATENT}",
    "rho_i missing (units)": f"{REF_STORAGE} - {REF_CONDUCTION} - L*Derivative({_P}, t)",
    "a constant source": f"{REF_STORAGE} - {REF_CONDUCTION} - {REF_LATENT} - q0",
}
EXTRA_SYMBOLS = [{"name": "q0", "meaning": "a heat source with no physical cause (planted)", "units": "W m^-3",
                  "typical": 1e3}]


def evaluate(text, symbols):
    """(all checks passed, check results, roles) for one answer."""
    units = dict(fixed_units())
    decl = declared(symbols)
    units.update({n: d[1] for n, d in decl.items()})
    expr = parse(text, decl)
    if not (expr.has(T) and expr.has(phi)):
        raise Reject("R must contain both T and phi")
    results = run_checks(expr, decl, units)
    return all(ok for _, ok, _ in results), results, roles(expr)


def key_check():
    """The key's pedigree, with no LLM in it: the current model's heat row (src/assembly.c:173-175, read by
    hand; you can check its three terms there) passes every check, and each of five corrupted copies fails
    at least one. So the checks tell a right energy equation from these wrong ones."""
    ok, results, _ = evaluate(REFERENCE, REF_SYMBOLS)
    if not ok:
        return False, f"the reference fails its own checks: {[r for r in results if not r[1]]}"
    for name, text in CORRUPTED.items():
        try:
            bad_ok, _, _ = evaluate(text, REF_SYMBOLS + EXTRA_SYMBOLS)
        except Reject:
            bad_ok = False
        if bad_ok:
            return False, f"the checks do not catch the corrupted reference '{name}'"
    return True, f"the reference passes all checks; {len(CORRUPTED)} corrupted references each fail one"


def grade(answer):
    """Return (pass, details). Raises nothing."""
    try:
        key_ok, why = key_check()
        if not key_ok:
            return False, f"grader fault: {why}"
        ok, results, groups = evaluate(answer.get("result_expr"), answer.get("symbols"))
        details = {
            "checks": [{"check": n, "pass": p, "details": d} for n, p, d in results],
            "comparison (for the user)": {
                "agent's terms by role": groups,
                "current model (src/assembly.c:173-175)": {
                    "storage": [REF_STORAGE], "conduction": [f"-{REF_CONDUCTION}"],
                    "phase change": [f"-{REF_LATENT}"], "other": []},
                "agent's symbols": answer.get("symbols"),
                "agent left out": answer.get("neglected"),
            },
        }
        if groups["other"]:
            details["for the user"] = (f"{len(groups['other'])} term(s) the current model does not have; "
                                       "judge whether each is right")
        return ok, details
    except Reject as e:
        return False, str(e)


def main(argv):
    if len(argv) != 2:
        print("usage: grade.py <result.json>", file=sys.stderr)
        return 2
    path = argv[1]
    try:
        raw = open(path, "rb").read()
        result = json.loads(raw)
    except (OSError, ValueError) as e:
        print(json.dumps({"check": CHECK, "pass": False, "details": f"cannot read {path}: {e}"}))
        return 2
    print(f"grading {path} (sha256 {hashlib.sha256(raw).hexdigest()[:12]}), sympy {sp.__version__}, "
          f"pint {pint.__version__}", file=sys.stderr)
    if result.get("status") != "done":
        ok, details = False, f"status is {result.get('status')!r}, not 'done'"
    else:
        ok, details = grade(result)
    print(json.dumps({"check": CHECK, "pass": ok, "details": details}, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
