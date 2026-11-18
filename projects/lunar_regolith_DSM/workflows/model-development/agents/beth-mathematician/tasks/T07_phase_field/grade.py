#!/usr/bin/env python3
"""Grader for T07: the evolution equation of the phase field, derived blind from the T06 free energy.

Usage: grade.py <result.json>

The agent never sees this file. It was given F (from T06) and told the physics of the driving force (ice
deposits where rho_v > rho_vs(T), sublimates where rho_v < rho_vs(T), no nucleation), not the equation. Two
parts:

1. Checks that need no key (all must pass). gamma, eps and the agent's parameters get values near their
   typical sizes; rho_vs(T) is a positive increasing test function near 1e-6 kg m^-3 at 200 K. phi_t is
   found by solving R = 0 for it.
   - names, units:  as in T04;
   - bulk:          uniform phi = 0 and phi = 1 stay put for vapor above and below saturation (no
                    nucleation in gas, no sublimation from inside ice; phi stays in [0, 1]);
   - growth:        at phi = 1/2 (where dF/dphi = 0 for a flat field), supersaturation gives phi_t > 0 and
                    undersaturation phi_t < 0;
   - relaxation:    at saturation (rho_v = rho_vs(T)), phi_t = -M mu with mu = dF/dphi = h'(phi) - kappa
                    lap(phi) and M > 0: for two fields with the same phi but different gradients and
                    curvature, phi_t / mu is the same negative number. So F can only decrease, and the
                    curvature term and the well enter in the ratio F fixes.
2. A comparison with the current model's phase-field row (src/assembly.c:177-181, phase-change switch on,
   clamps left out): R_ref = phi_t - 3 M eps lap(phi) + (3 M / eps) phi (1 - phi)(1 - 2 phi)
   - (alpha / rho_i) phi^2 (1 - phi)^2 (rho_v - rho_vs(T)), with M = mob_sub constant
   (src/lunar_main.c:1103). Its relaxation part is (M / (2 gamma)) mu for the F of T06.

key_check, with no LLM in it: R_ref passes every check, and each of six corrupted copies fails one.

Prints one JSON line {"check", "pass", "details"}; exits 0 on pass, 1 on fail, 2 on an unreadable file.
"""
import hashlib
import json
import random
import sys
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_shared"))
import blind  # noqa: E402
from blind import RHOV_SCALE, T_SCALE, Reject, T, at_origin, lam, phi, rho_v, solve_rate, t, x, y, z  # noqa: E402

CHECK = "T07_phase_field"
FIXED = blind.fixed_symbols(HERE)
GAMMA, EPS, RHO_VS = FIXED["gamma"], FIXED["eps"], FIXED["rho_vs"]
_T = sp.Symbol("_T", real=True)
SAT = sp.Lambda(_T, RHOV_SCALE * (1 + (_T - T_SCALE) / 20))     # positive and increasing near 200 K
p, L = sp.symbols("p lap", real=True)
A = sp.Symbol("a", real=True)


def mu(vals):
    """dF/dphi of the given F at phi = p with Laplacian lap, for the given gamma and eps."""
    h = 3 * GAMMA / EPS * p**2 * (1 - p)**2
    kappa = 6 * GAMMA * EPS
    return (sp.diff(h, p) - kappa * L).subs(vals)


def scales(seed):
    rng = random.Random(seed)
    return {GAMMA: sp.Rational(rng.randint(5, 20), 100), EPS: sp.Rational(rng.randint(1, 100), 10**8)}


def phi_rate(expr, decl, seed, phi_body, rv_factor, T_body=T_SCALE + 0 * x):
    """phi_t at the origin for phi = phi_body + a t, T = T_body and rho_v = rv_factor * rho_vs(T)."""
    e = expr.subs(scales(seed))
    rv = rv_factor * SAT(T_body)
    e = blind.concrete(e, decl, {phi: lam(phi_body + A * t), T: lam(T_body), rho_v: lam(rv), RHO_VS: SAT}, seed)
    return solve_rate(at_origin(e), A)


def check_bulk(expr, decl):
    for seed in range(3):
        for p0 in (0, 1):
            for f in (sp.Rational(3, 2), sp.Rational(1, 2)):
                rate = phi_rate(expr, decl, seed, sp.Integer(p0) + 0 * x, f)
                if rate != 0:
                    return False, (f"uniform phi = {p0} with rho_v = {f} rho_vs gives phi_t = {rate}: "
                                   f"{'ice appears in gas' if p0 == 0 else 'ice changes inside bulk ice'}, or "
                                   "phi leaves [0, 1]")
    return True, "pure gas and pure ice stay put above and below saturation (3 draws)"


def check_growth(expr, decl):
    for seed in range(3):
        up = phi_rate(expr, decl, seed, sp.Rational(1, 2) + 0 * x, sp.Rational(6, 5))
        down = phi_rate(expr, decl, seed, sp.Rational(1, 2) + 0 * x, sp.Rational(4, 5))
        if not (up > 0 and down < 0):
            return False, (f"at phi = 1/2, rho_v = 1.2 rho_vs gives phi_t = {up} and 0.8 rho_vs gives {down}: "
                           "supersaturation must deposit ice, undersaturation sublimate it")
    return True, "supersaturation deposits ice and undersaturation sublimates it (3 draws)"


def check_relaxation(expr, decl):
    shapes = ((1, 3, -1, 2), (-2, -4, 5, 1))          # (grad x, curvature x, y, z) of two test fields
    for seed in range(3):
        vals = scales(seed)
        e0 = vals[EPS]
        for p0 in (sp.Rational(3, 10), sp.Rational(7, 10)):
            ratios = []
            for c1, c2, c3, c4 in shapes:
                body = p0 + c1 * x / e0 / 10 + (c2 * x**2 + c3 * y**2 + c4 * z**2) / e0**2 / 100
                rate = phi_rate(expr, decl, seed, body, 1)
                m = mu(vals).subs({p: p0, L: 2 * sp.Rational(c2 + c3 + c4) / e0**2 / 100})
                if m == 0:
                    raise Reject("grader fault: a relaxation test field has mu = 0")
                ratios.append(sp.nsimplify(rate / m))
            if not ratios[0] < 0:
                return False, (f"at saturation phi_t / (dF/dphi) = {float(ratios[0]):.4g} at phi = {p0}: "
                               "the interface must relax downhill in F (ratio < 0)")
            if ratios[0] != ratios[1]:
                return False, (f"at saturation phi_t / (dF/dphi) is {[float(r) for r in ratios]} at phi = {p0} "
                               "for two fields with the same phi: phi_t is not -M dF/dphi (curvature and well "
                               "are not in the ratio F gives)")
    return True, "at saturation phi_t = -M dF/dphi with M > 0, so F decreases (3 draws, phi = 0.3 and 0.7)"


def run_checks(expr, decl, units):
    results = []
    try:
        q = blind.dims(expr, units)
        results.append(("units", True, f"every term has units {q.units:~}"))
    except Reject as e:
        results.append(("units", False, str(e)))
    for name, fn in (("bulk", check_bulk), ("growth", check_growth), ("relaxation", check_relaxation)):
        try:
            results.append((name, *fn(expr, decl)))
        except Reject as e:
            results.append((name, False, str(e)))
        except Exception as e:      # the agent's expression broke a check: report it, do not crash
            results.append((name, False, f"could not evaluate: {type(e).__name__}: {str(e)[:160]}"))
    return results


def roles(expr):
    groups = {"rate": [], "relaxation": [], "driving force": [], "other": []}
    for term in sp.Add.make_args(expr):          # as written, not expanded, so each term stays readable
        if term.has(rho_v) or term.has(RHO_VS):
            groups["driving force"].append(term)
        elif any(d.expr.has(phi) and t in {v for v, _ in d.variable_count} for d in term.atoms(sp.Derivative)):
            groups["rate"].append(term)
        elif term.has(phi):
            groups["relaxation"].append(term)
        else:
            groups["other"].append(term)
    return {k: [str(v) for v in vs] for k, vs in groups.items()}


REF_SYMBOLS = [
    {"name": "M", "meaning": "interface mobility (mob_sub = eps / (3 tau))", "units": "m s^-1", "typical": 1e-4},
    {"name": "alpha", "meaning": "phase-change rate (alph_sub)", "units": "s^-1", "typical": 1e3},
    {"name": "rho_i", "meaning": "density of ice", "units": "kg m^-3", "typical": 917},
]
_P, _T4 = "phi(x, y, z, t)", "T(x, y, z, t)"
_LAP = "(" + " + ".join(f"Derivative({_P}, ({c}, 2))" for c in "xyz") + ")"
REF_RATE = f"Derivative({_P}, t)"
REF_RELAX = f"- 3*M*eps*{_LAP} + 3*M/eps*{_P}*(1 - {_P})*(1 - 2*{_P})"
REF_DRIVE = f"alpha/rho_i*{_P}**2*(1 - {_P})**2*(rho_v(x, y, z, t) - rho_vs({_T4}))"
REFERENCE = f"{REF_RATE} {REF_RELAX} - {REF_DRIVE}"
CORRUPTED = {
    "driving force with the wrong sign": f"{REF_RATE} {REF_RELAX} + {REF_DRIVE}",
    "driving force everywhere (nucleation)": f"{REF_RATE} {REF_RELAX} - alpha/rho_i*(rho_v(x, y, z, t) - rho_vs({_T4}))",
    "curvature term doubled": f"{REF_RATE} - 6*M*eps*{_LAP} + 3*M/eps*{_P}*(1 - {_P})*(1 - 2*{_P}) - {REF_DRIVE}",
    "relaxation uphill": f"{REF_RATE} + 3*M*eps*{_LAP} - 3*M/eps*{_P}*(1 - {_P})*(1 - 2*{_P}) - {REF_DRIVE}",
    "no relaxation": f"{REF_RATE} - {REF_DRIVE}",
    "wrong units": f"{REF_RATE} {REF_RELAX} - alpha*{_P}**2*(1 - {_P})**2*(rho_v(x, y, z, t) - rho_vs({_T4}))",
}


def evaluate(text, symbols):
    units = blind.fixed_units(HERE)
    decl = blind.declared(symbols, FIXED)
    units.update({n: d[1] for n, d in decl.items()})
    expr = blind.parse(text, decl, FIXED)
    if not expr.has(phi):
        raise Reject("R must contain phi")
    results = run_checks(expr, decl, units)
    return all(ok for _, ok, _ in results), results, roles(expr)


def key_check():
    """The key's pedigree, with no LLM in it: the current model's phase-field row (src/assembly.c:177-181, read
    by hand; you can check its four terms there) passes every check, and each of six corrupted copies fails
    at least one."""
    ok, results, _ = evaluate(REFERENCE, REF_SYMBOLS)
    if not ok:
        return False, f"the reference fails its own checks: {[r for r in results if not r[1]]}"
    for name, text in CORRUPTED.items():
        try:
            bad_ok, _, _ = evaluate(text, REF_SYMBOLS)
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
            "checks": [{"check": n, "pass": q, "details": d} for n, q, d in results],
            "comparison (for the user)": {
                "agent's terms by role": groups,
                "current model (src/assembly.c:177-181)": {
                    "rate": [REF_RATE], "relaxation": [REF_RELAX], "driving force": [f"-{REF_DRIVE}"],
                    "other": []},
                "agent's symbols": answer.get("symbols"),
                "agent left out": answer.get("neglected"),
            },
        }
        if groups["other"]:
            details["for the user"] = f"{len(groups['other'])} term(s) the current model does not have"
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
    print(f"grading {path} (sha256 {hashlib.sha256(raw).hexdigest()[:12]}), sympy {sp.__version__}",
          file=sys.stderr)
    if result.get("status") != "done":
        ok, details = False, f"status is {result.get('status')!r}, not 'done'"
    else:
        ok, details = grade(result)
    print(json.dumps({"check": CHECK, "pass": ok, "details": details}, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
