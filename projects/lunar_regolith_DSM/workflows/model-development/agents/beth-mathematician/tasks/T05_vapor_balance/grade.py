#!/usr/bin/env python3
"""Grader for T05: the water-mass balance for the vapor density rho_v, derived blind.

Usage: grade.py <result.json>

The agent never sees this file. It was told the problem, not the equation, and names its own symbols, each
with units and a typical SI value (shared machinery: ../_shared/blind.py). Two parts:

1. Checks that need no key (all must pass), near T = 200 K and rho_v = 1e-6 kg m^-3:
   - names, units:  as in T04;
   - equilibrium:   uniform rho_v and T with a static phi give R = 0 (no vapor out of nothing);
   - exchange:      uniform phi(t) = phi0 + t and rho_v(t) = rho0 + a t. Water is conserved only if the vapor
                    lost, d/dt[(1 - phi) rho_v] = (1 - phi0) a - rho0, equals the ice gained, rho* times
                    phi_t, for one ice density rho*. So rho* = rho0 - (1 - phi0) a must be the same (relative
                    1e-6) for three choices of (phi0, rho0), and ice-like (100 to 5000 kg m^-3). This checks
                    the storage weight (1 - phi) and the exchange term together;
   - diffusion:     a vapor peak (phi uniform and static) decays: rho_v,t < 0 at the peak;
   - conservation:  with phi, rho_v and T static and periodic in x, R averages to zero over a period, so the
                    flux is a divergence and moves vapor without making it.
2. A comparison with the current model (src/assembly.c:188-191, with xi_v = 1, the phase-change switch on,
   and the numerical floor air_lim left out): R_ref = (1 - phi) rho_v,t - div(D(T) (1 - phi) grad rho_v)
   + (rho_i - rho_v) phi_t, grouped by role (storage, diffusion, phase change, other) for the user.

key_check, with no LLM in it: R_ref passes every check, and each of six corrupted copies fails one.

Prints one JSON line {"check", "pass", "details"}; exits 0 on pass, 1 on fail, 2 on an unreadable file.
"""
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_shared"))
import blind  # noqa: E402
from blind import (COORDS, RHOV_SCALE, T_SCALE, Reject, T, at_origin, lam, periodic_mean,  # noqa: E402
                   phi, rho_v, solve_rate, t, x, y, z)

CHECK = "T05_vapor_balance"
FIXED = blind.fixed_symbols(HERE)
R0 = RHOV_SCALE


def fields(T_body, phi_body, rv_body):
    return {T: lam(T_body), phi: lam(phi_body), rho_v: lam(rv_body)}


def check_equilibrium(expr, decl):
    for seed in range(3):
        e = blind.concrete(expr, decl, fields(T_SCALE + 0 * x, sp.Rational(1, 2) + (x - y + 2 * z) / 10, R0 + 0 * x),
                           seed)
        for pt in ((0, 0, 0), (sp.Rational(1, 3), sp.Rational(-1, 5), sp.Rational(1, 7))):
            v = sp.simplify(e.subs(dict(zip(COORDS, pt))).subs(t, 0))
            if v != 0:
                return False, f"uniform rho_v and T with a static phi give R = {v}, not 0: vapor out of nothing"
    return True, "R = 0 for uniform rho_v and T and a static phi (3 draws, 2 points)"


def check_exchange(expr, decl):
    a = sp.Symbol("a")
    cases = ((sp.Rational(3, 10), R0), (sp.Rational(6, 10), 2 * R0), (sp.Rational(1, 2), R0 / 2))
    for seed in range(3):
        implied = []
        for phi0, rho0 in cases:
            e = blind.concrete(expr, decl, fields(T_SCALE + 0 * x, phi0 + t, rho0 + a * t), seed)
            rate = solve_rate(at_origin(e), a)
            implied.append(rho0 - (1 - phi0) * rate)
        lo, hi = min(implied), max(implied)
        if not (100 <= lo and hi <= 5000):
            return False, (f"the ice density implied by the exchange is {[float(v) for v in implied]} kg m^-3 "
                           f"(draw {seed}): water is not conserved, or the exchange has the wrong sign or size")
        if (hi - lo) > sp.Rational(1, 10**6) * hi:
            return False, (f"the implied ice density changes with phi and rho_v: {[float(v) for v in implied]} "
                           f"(draw {seed}): the vapor storage or the exchange term does not conserve water")
    return True, (f"deposition consumes vapor at one ice-like density (draw 2: {float(implied[0]):.6g} kg m^-3) "
                  "for every phi and rho_v: water is conserved")


def check_diffusion(expr, decl):
    a = sp.Symbol("a")
    for seed in range(3):
        e = blind.concrete(expr, decl, fields(T_SCALE + 0 * x, sp.Rational(1, 2) + 0 * x,
                                              R0 - R0 / 10 * (x**2 + y**2 + z**2) + a * t), seed)
        rate = solve_rate(at_origin(e), a)
        if not rate < 0:
            return False, f"a vapor peak gives rho_v,t = {rate} at its peak: vapor should spread out (draw {seed})"
    return True, "a vapor peak decays (3 draws)"


def check_conservation(expr, decl):
    for seed in range(3):
        e = blind.concrete(expr, decl, fields(T_SCALE + 5 * sp.sin(2 * sp.pi * x),
                                              sp.Rational(1, 2) + sp.sin(2 * sp.pi * x) / 4,
                                              # shares a mode with phi, so a misplaced (1 - phi) shows up
                                              R0 * (1 + sp.sin(2 * sp.pi * x) / 2 + sp.cos(4 * sp.pi * x) / 5)),
                           seed)
        mean, scale = periodic_mean(e)
        if scale == 0:
            return False, "R vanishes for a non-uniform rho_v: no transport term"
        if abs(mean) > 1e-9 * scale:
            return False, (f"with static, periodic fields R averages {mean:.3g} over a period (max |R| "
                           f"{scale:.3g}): the vapor flux is not a divergence, so it makes or destroys vapor")
    return True, "transport integrates to zero over a period (3 draws): the flux is a divergence"


def run_checks(expr, decl, units):
    results = []
    try:
        q = blind.dims(expr, units)
        results.append(("units", True, f"every term has units {q.units:~}"))
    except Reject as e:
        results.append(("units", False, str(e)))
    for name, fn in (("equilibrium", check_equilibrium), ("exchange", check_exchange),
                     ("diffusion", check_diffusion), ("conservation", check_conservation)):
        try:
            results.append((name, *fn(expr, decl)))
        except Reject as e:
            results.append((name, False, str(e)))
        except Exception as e:      # the agent's expression broke a check: report it, do not crash
            results.append((name, False, f"could not evaluate: {type(e).__name__}: {str(e)[:160]}"))
    return results


def roles(expr):
    groups = {"storage": [], "diffusion": [], "phase change": [], "other": []}
    for term in sp.Add.make_args(sp.expand(expr)):
        kinds = set()
        for d in term.atoms(sp.Derivative):
            wrt = {v for v, _ in d.variable_count}
            if d.expr.has(rho_v) and wrt & set(COORDS):
                kinds.add("diffusion")
            elif d.expr.has(rho_v) and t in wrt:
                kinds.add("storage")
            elif d.expr.has(phi) and t in wrt:
                kinds.add("phase change")
            elif d.expr.has(phi) and wrt & set(COORDS) and term.has(rho_v):
                kinds.add("diffusion")          # (D grad phi) . grad rho_v from expanding the divergence
            else:
                kinds.add("other")
        groups[kinds.pop() if len(kinds) == 1 else "other"].append(term)
    return {k: [str(v) for v in vs] for k, vs in groups.items()}


# ---- the reference: the current model's vapor row --------------------------------------------------------
REF_SYMBOLS = [
    {"name": "D(T)", "meaning": "vapor diffusivity in air", "units": "m^2 s^-1", "typical": 1e-5},
    {"name": "rho_i", "meaning": "density of ice", "units": "kg m^-3", "typical": 917},
]
_P, _R, _T = "phi(x, y, z, t)", "rho_v(x, y, z, t)", "T(x, y, z, t)"
REF_STORAGE = f"(1 - {_P})*Derivative({_R}, t)"
REF_DIFFUSION = "(" + " + ".join(f"Derivative(D({_T})*(1 - {_P})*Derivative({_R}, {c}), {c})" for c in "xyz") + ")"
REF_EXCHANGE = f"(rho_i - {_R})*Derivative({_P}, t)"
REFERENCE = f"{REF_STORAGE} - {REF_DIFFUSION} + {REF_EXCHANGE}"
CORRUPTED = {
    "exchange sign flipped": f"{REF_STORAGE} - {REF_DIFFUSION} - {REF_EXCHANGE}",
    "diffusion sign flipped": f"{REF_STORAGE} + {REF_DIFFUSION} + {REF_EXCHANGE}",
    "gas fraction outside the divergence": f"{REF_STORAGE} - (1 - {_P})*(" + " + ".join(
        f"Derivative(D({_T})*Derivative({_R}, {c}), {c})" for c in "xyz") + f") + {REF_EXCHANGE}",
    "storage not weighted by the gas fraction": f"Derivative({_R}, t) - {REF_DIFFUSION} + {REF_EXCHANGE}",
    "exchange with the wrong units": f"{REF_STORAGE} - {REF_DIFFUSION} + D({_T})*Derivative({_P}, t)",
    "a constant source": f"{REF_STORAGE} - {REF_DIFFUSION} + {REF_EXCHANGE} - s0",
}
EXTRA_SYMBOLS = [{"name": "s0", "meaning": "a vapor source with no physical cause (planted)",
                  "units": "kg m^-3 s^-1", "typical": 1e-9}]


def evaluate(text, symbols):
    """(all checks passed, check results, roles) for one answer."""
    units = blind.fixed_units(HERE)
    decl = blind.declared(symbols, FIXED)
    units.update({n: d[1] for n, d in decl.items()})
    expr = blind.parse(text, decl, FIXED)
    if not expr.has(rho_v):
        raise Reject("R must contain rho_v")
    results = run_checks(expr, decl, units)
    return all(ok for _, ok, _ in results), results, roles(expr)


def key_check():
    """The key's pedigree, with no LLM in it: the current model's vapor row (src/assembly.c:188-191, read by
    hand; you can check its three terms there) passes every check, and each of six corrupted copies fails at
    least one. So the checks tell a right water balance from these wrong ones."""
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
                "current model (src/assembly.c:188-191)": {
                    "storage": [REF_STORAGE], "diffusion": [f"-{REF_DIFFUSION}"],
                    "phase change": [REF_EXCHANGE], "other": []},
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
          f"pint {blind.pint.__version__}", file=sys.stderr)
    if result.get("status") != "done":
        ok, details = False, f"status is {result.get('status')!r}, not 'done'"
    else:
        ok, details = grade(result)
    print(json.dumps({"check": CHECK, "pass": ok, "details": details}, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
