#!/usr/bin/env python3
"""Grader for T06: the free-energy density f of the ice-gas phase field, derived blind.

Usage: grade.py <result.json>

The agent never sees this file. It was told what the free energy must do (stable bulk ice and gas, an
interface of finite width with energy gamma per area, width proportional to eps), not its form. Two parts:

1. Checks that need no key (all must pass). h(p) is f with phi = p and no gradients; the gradient part is
   f - h. gamma and eps get random positive values per draw (0.05 to 0.2 J m^-2, 1e-8 to 1e-6 m).
   - names, units:   f is an energy per volume (pint); first derivatives of phi only;
   - wells:          h'(0) = h'(1) = 0, h''(0) > 0, h''(1) > 0, h(0) = h(1) (coexistence), and h > h(0) on
                     (0, 1): bulk gas and bulk ice are the two stable states, separated by a barrier;
   - gradient:       f - h = kappa(p) |grad phi|^2 / 2 with kappa > 0: the same for every direction of the
                     gradient, and four times larger for a gradient twice as steep (isotropic, quadratic,
                     penalizing);
   - interface:      the energy per area of the flat equilibrium interface of this f,
                     sigma = integral_0^1 sqrt(2 kappa(p) (h(p) - h(0))) dp, equals gamma to 1e-8.
   Reported, not failed: the 10-90 % width of that interface in units of eps, and the barrier height.
2. A comparison with the free energy that the current model's phase-field row implies (src/assembly.c:177-181
   with the well comment at :7-12): f_ref = 3 gamma eps |grad phi|^2 + 3 gamma phi^2 (1 - phi)^2 / eps. This
   is an inference, not a line of code: the code has only the evolution equation, with -3 M eps lap(phi) +
   (3 M / eps) phi (1 - phi)(1 - 2 phi); f_ref is the energy whose variational derivative has that shape and
   whose flat interface has energy gamma. Its 10-90 % width is 2 ln 9 eps = 4.39 eps.

key_check, with no LLM in it: f_ref passes every check, and each of six corrupted copies fails one.

Prints one JSON line {"check", "pass", "details"}; exits 0 on pass, 1 on fail, 2 on an unreadable file.
"""
import hashlib
import json
import random
import sys
from pathlib import Path

import mpmath
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_shared"))
import blind  # noqa: E402
from blind import COORDS, Reject, phi, t  # noqa: E402

CHECK = "T06_free_energy"
FIXED = blind.fixed_symbols(HERE)
GAMMA, EPS = FIXED["gamma"], FIXED["eps"]
p = sp.Symbol("p", real=True)
G = sp.symbols("gx gy gz", real=True)
PROBES = [sp.Rational(k, 10) for k in range(1, 10)]


def split(expr):
    """(h(p), gradient part as a function of p and gx, gy, gz) from f."""
    for d in expr.atoms(sp.Derivative):
        if d.expr != phi(*COORDS, t) or sum(n for _, n in d.variable_count) != 1 or d.variables[0] not in COORDS:
            raise Reject(f"f may hold phi and its first space derivatives only; found {d}")
    with_g = expr.subs({sp.Derivative(phi(*COORDS, t), c): g for c, g in zip(COORDS, G)}).subs(phi(*COORDS, t), p)
    if with_g.has(*COORDS, t):
        raise Reject("f depends on position or time other than through phi")
    h = with_g.subs({g: 0 for g in G})
    return sp.simplify(h), sp.simplify(with_g - h)


def draws(decl):
    """Three sets of values: gamma, eps and the agent's own parameters (typical value times 1/2 to 2)."""
    out = []
    for seed in range(3):
        rng = random.Random(seed)
        vals = {GAMMA: sp.Rational(rng.randint(5, 20), 100), EPS: sp.Rational(rng.randint(1, 100), 10**8)}
        for obj, _, typical in decl.values():
            if not isinstance(obj, sp.Symbol):
                raise Reject(f"f may not hold a declared function ({obj}): at coexistence its coefficients are "
                             "constants")
            vals[obj] = typical * sp.Rational(rng.randint(2, 8), 4)
        out.append(vals)
    return out


def check_wells(h, vals_list):
    for vals in vals_list:
        hv = h.subs(vals)
        d1, d2 = sp.diff(hv, p), sp.diff(hv, p, 2)
        problems = []
        for end in (0, 1):
            if sp.simplify(d1.subs(p, end)) != 0:
                problems.append(f"h'({end}) = {sp.simplify(d1.subs(p, end))}, not 0")
            if not sp.simplify(d2.subs(p, end)) > 0:
                problems.append(f"h''({end}) = {sp.simplify(d2.subs(p, end))}, not > 0")
        if sp.simplify(hv.subs(p, 1) - hv.subs(p, 0)) != 0:
            problems.append(f"h(1) - h(0) = {sp.simplify(hv.subs(p, 1) - hv.subs(p, 0))}: not at coexistence")
        low = [q for q in PROBES if not hv.subs(p, q) > hv.subs(p, 0)]
        if low:
            problems.append(f"h is not above h(0) at p = {[str(q) for q in low]}: no barrier between the phases")
        if problems:
            return False, "; ".join(problems)
    return True, "minima at 0 and 1 with equal depths and a barrier between (3 draws)"


def kappa_of(grad, vals):
    """kappa(p) if grad = kappa(p) |g|^2 / 2 for every direction and size of g, else Reject."""
    gv = grad.subs(vals)
    unit = sp.simplify(2 * gv.subs({G[0]: 1, G[1]: 0, G[2]: 0}))
    for vec in ((0, 1, 0), (0, 0, 1), (1, 2, -2), (3, 0, 4)):
        n2 = sum(c * c for c in vec)
        if sp.simplify(2 * gv.subs(dict(zip(G, vec))) - unit * n2) != 0:
            raise Reject(f"the gradient part is not kappa |grad phi|^2 / 2: direction {vec} differs")
    return unit


def check_gradient(grad, vals_list):
    for vals in vals_list:
        kappa = kappa_of(grad, vals)
        bad = [q for q in PROBES + [0, 1] if not kappa.subs(p, q) > 0]
        if bad:
            return False, f"kappa(p) is not positive at p = {[str(q) for q in bad]}: gradients would lower F"
    return True, "isotropic, quadratic and penalizing: kappa(p) |grad phi|^2 / 2 with kappa > 0 (3 draws)"


def interface(h, grad, vals):
    """(energy per area, 10-90 % width) of the flat equilibrium interface, numerically."""
    kappa = kappa_of(grad, vals)
    hv = h.subs(vals)
    dh = sp.lambdify(p, hv - hv.subs(p, 0), "mpmath")
    kf = sp.lambdify(p, kappa, "mpmath")
    mpmath.mp.dps = 30
    sigma = mpmath.quad(lambda q: mpmath.sqrt(2 * kf(q) * dh(q)), [0, 0.5, 1])
    width = mpmath.quad(lambda q: mpmath.sqrt(kf(q) / (2 * dh(q))), [0.1, 0.5, 0.9])
    return sigma, width


def check_interface(h, grad, vals_list):
    widths = []
    for vals in vals_list:
        sigma, width = interface(h, grad, vals)
        gamma = mpmath.mpf(str(float(vals[GAMMA])))
        if abs(sigma / gamma - 1) > 1e-8:
            return False, f"the flat interface has energy {float(sigma / gamma):.6g} gamma, not gamma"
        widths.append(float(width / mpmath.mpf(str(float(vals[EPS])))))
    return True, f"the flat interface has energy gamma (3 draws); its 10-90 % width is {widths[0]:.4g} eps"


def run_checks(expr, decl, units):
    results = []
    try:
        q = blind.dims(expr, units)
        results.append(("units", q.dimensionality == blind.UREG("J/m**3").dimensionality,
                        f"f has units {q.units:~} (needs an energy per volume)"))
    except Reject as e:
        results.append(("units", False, str(e)))
    try:
        h, grad = split(expr)
        vals_list = draws(decl)
    except Reject as e:
        return results + [("form", False, str(e))]
    for name, fn in (("wells", lambda: check_wells(h, vals_list)),
                     ("gradient", lambda: check_gradient(grad, vals_list)),
                     ("interface", lambda: check_interface(h, grad, vals_list))):
        try:
            results.append((name, *fn()))
        except Reject as e:
            results.append((name, False, str(e)))
        except Exception as e:      # the agent's expression broke a check: report it, do not crash
            results.append((name, False, f"could not evaluate: {type(e).__name__}: {str(e)[:160]}"))
    return results


_P = "phi(x, y, z, t)"
_G2 = " + ".join(f"Derivative({_P}, {c})**2" for c in "xyz")
REFERENCE = f"3*gamma*eps*({_G2}) + 3*gamma*{_P}**2*(1 - {_P})**2/eps"
CORRUPTED = {
    "wells of unequal depth": f"{REFERENCE} + gamma/eps*{_P}/10",
    "single well": f"3*gamma*eps*({_G2}) + 3*gamma*{_P}**2/eps",
    "gradient term with the wrong sign": f"-3*gamma*eps*({_G2}) + 3*gamma*{_P}**2*(1 - {_P})**2/eps",
    "anisotropic gradient": f"3*gamma*eps*Derivative({_P}, x)**2 + 3*gamma*{_P}**2*(1 - {_P})**2/eps",
    "interface energy 2 gamma": f"6*gamma*eps*({_G2}) + 6*gamma*{_P}**2*(1 - {_P})**2/eps",
    "wrong units": f"3*gamma*({_G2}) + 3*gamma*{_P}**2*(1 - {_P})**2/eps",
}


def evaluate(text, symbols):
    units = blind.fixed_units(HERE)
    decl = blind.declared(symbols, FIXED)
    units.update({n: d[1] for n, d in decl.items()})
    expr = blind.parse(text, decl, FIXED)
    if not expr.has(phi):
        raise Reject("f must contain phi")
    results = run_checks(expr, decl, units)
    return all(ok for _, ok, _ in results), results


def key_check():
    """The key's pedigree, with no LLM in it: the free energy the current model's phase-field row implies
    (src/assembly.c:177-181; see the module docstring for why it is an inference) passes every check, its
    10-90 % width is 2 ln 9 eps as the closed form says, and each of six corrupted copies fails a check."""
    ok, results = evaluate(REFERENCE, [])
    if not ok:
        return False, f"the reference fails its own checks: {[r for r in results if not r[1]]}"
    h, grad = split(blind.parse(REFERENCE, {}, FIXED))
    _, width = interface(h, grad, draws({})[0])
    if abs(float(width) / float(draws({})[0][EPS]) - 2 * float(sp.log(9))) > 1e-8:
        return False, "the reference's 10-90 % width is not 2 ln 9 eps"
    for name, text in CORRUPTED.items():
        try:
            bad_ok, _ = evaluate(text, [])
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
        ok, results = evaluate(answer.get("result_expr"), answer.get("symbols"))
        return ok, {
            "checks": [{"check": n, "pass": q, "details": d} for n, q, d in results],
            "comparison (for the user)": {
                "agent's f": answer.get("result_expr"),
                "implied by the current model (src/assembly.c:177-181)": REFERENCE,
                "agent's symbols": answer.get("symbols"),
                "agent left out": answer.get("neglected"),
            },
        }
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
