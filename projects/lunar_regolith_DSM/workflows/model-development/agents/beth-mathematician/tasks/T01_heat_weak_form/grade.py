#!/usr/bin/env python3
"""Grader for T01: is the agent's R a weak form of the 1D transient heat equation?

Usage: grade.py <result.json>

The agent never sees this file. Its answer is untrusted text, so only the allowed names, numbers and
arithmetic are accepted before SymPy parses it (SymPy's parser evaluates what it reads).

The key itself is checked first, with no LLM in the loop (key_check): it must equal the integral of w times
the given strong form, in exact arithmetic. If it does not, the grader refuses to grade.

Two checks of the answer:
1. Structure: no second x-derivative of T remains, and no derivative is applied to an expression that
   already holds a derivative of T (the conduction term really was integrated by parts).
2. Value: with random polynomial T, w, k, Q and random rho, c, L, t, the answer equals a constant,
   nonzero multiple of the reference weak form. A multiple is accepted because R = 0 and -R = 0 are the
   same weak form; a wrong sign on one term, a missing factor or a missing boundary term are not.

Prints one JSON line {"check", "pass", "details"}; exits 0 on pass, 1 on fail, 2 on an unreadable file.
"""
import hashlib
import json
import random
import re
import sys

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

CHECK = "T01_heat_weak_form"
ALLOWED = {"x", "t", "L", "T", "w", "k", "Q", "rho", "c", "Integral", "Derivative", "Subs"}
MAX_LEN = 4000
N_SAMPLES = 6

x, t, L, rho, c = sp.symbols("x t L rho c")
T, w, k, Q = (sp.Function(n) for n in ("T", "w", "k", "Q"))
LOCAL = {"x": x, "t": t, "L": L, "rho": rho, "c": c, "T": T, "w": w, "k": k, "Q": Q,
         "Integral": sp.Integral, "Derivative": sp.Derivative, "Subs": sp.Subs}


class Reject(Exception):
    pass


def reference():
    """R = ∫ ρc w T_t dx + ∫ w_x k T_x dx − ∫ w Q dx − [w k T_x] from 0 to L."""
    tx = sp.Derivative(T(x, t), x)
    flux = w(x) * k(x) * tx
    return (sp.Integral(rho * c * w(x) * sp.Derivative(T(x, t), t), (x, 0, L))
            + sp.Integral(sp.Derivative(w(x), x) * k(x) * tx, (x, 0, L))
            - sp.Integral(w(x) * Q(x), (x, 0, L))
            - (sp.Subs(flux, x, L) - sp.Subs(flux, x, 0)))


def key_check():
    """The key's pedigree, with no LLM in it: for concrete polynomial T, w, k, Q and numbers, the reference
    equals the integral of w times the given strong-form residual, in exact rational arithmetic. So the
    integration by parts and the boundary term in the key are right. Returns (ok, details)."""
    strong = (rho * c * sp.Derivative(T(x, t), t) - sp.Derivative(k(x) * sp.Derivative(T(x, t), x), x)
              - Q(x))
    weighted = sp.Integral(w(x) * strong, (x, 0, L))
    for seed in range(100, 103):
        funcs, nums = sample(seed)
        if evaluate(reference(), funcs, nums) != evaluate(weighted, funcs, nums):
            return False, f"the key differs from the integral of w times the strong form (sample {seed})"
    return True, "key equals the integral of w times the strong-form residual on 3 exact samples"


def parse(text):
    if not isinstance(text, str) or not text.strip():
        raise Reject("empty result_expr")
    if len(text) > MAX_LEN:
        raise Reject(f"result_expr longer than {MAX_LEN} characters")
    if not re.fullmatch(r"[A-Za-z0-9_+\-*/(),. \t\n]*", text):
        raise Reject("characters outside names, numbers, + - * / ** ( ) , and spaces")
    if re.search(r"(?<![0-9])\.|\.(?![0-9])", text):
        raise Reject("'.' is allowed only inside a decimal number")
    names = set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", text))
    if names - ALLOWED:
        raise Reject(f"names not allowed: {sorted(names - ALLOWED)}")
    try:
        return parse_expr(text, local_dict=dict(LOCAL), evaluate=True)
    except Exception as e:      # a parse error is the agent's fault, not the grader's
        raise Reject(f"does not parse: {type(e).__name__}: {e}") from None


def structure(expr):
    if not expr.has(w):
        raise Reject("the test function w does not appear")
    for d in expr.atoms(sp.Derivative):
        if not d.expr.has(T):
            continue
        if d.expr.atoms(sp.Derivative):
            raise Reject(f"a derivative is applied to a derivative of T: {d}; integrate by parts")
        if sum(n for v, n in d.variable_count if v == x) > 1:
            raise Reject(f"second x-derivative of T remains: {d}")


def sample(seed):
    rng = random.Random(seed)
    def r():
        return sp.Rational(rng.randint(-9, 9), rng.randint(1, 5))
    funcs = {
        T: sp.Lambda((x, t), sum(r() * x**i * t**j for i in range(4) for j in range(2))),
        w: sp.Lambda(x, 1 + sum(r() * x**i for i in range(3))),
        k: sp.Lambda(x, 2 + rng.randint(0, 3) * x + rng.randint(0, 2) * x**2),
        Q: sp.Lambda(x, sum(r() * x**i for i in range(3))),
    }
    nums = {L: sp.Rational(rng.randint(1, 4), rng.randint(1, 3)), rho: rng.randint(1, 5),
            c: rng.randint(1, 5), t: sp.Rational(rng.randint(0, 5), 2)}
    return funcs, nums


def evaluate(expr, funcs, nums):
    e = expr
    for f, lam in funcs.items():
        e = e.replace(f, lam)
    e = sp.simplify(e.doit().subs(nums))
    if e.free_symbols or not e.is_number:
        raise Reject(f"does not reduce to a number for concrete T, w, k, Q (left: {e})")
    return e


def grade(text):
    """Return (pass, details). Raises nothing."""
    try:
        key_ok, why = key_check()
        if not key_ok:
            return False, f"grader fault: {why}"
        expr = parse(text)
        structure(expr)
        ref = reference()
        ratios = []
        for seed in range(N_SAMPLES):
            funcs, nums = sample(seed)
            v_ref = evaluate(ref, funcs, nums)
            if v_ref == 0:
                continue
            ratios.append(sp.nsimplify(evaluate(expr, funcs, nums) / v_ref))
        if len(ratios) < 4:
            return False, "grader fault: too few usable samples"
        if ratios[0] == 0:
            return False, "the answer is zero where the weak form is not"
        if any(q != ratios[0] for q in ratios[1:]):
            return False, f"not a constant multiple of the weak form (ratios {[str(q) for q in ratios]})"
        return True, f"equals {ratios[0]} times the weak form on {len(ratios)} random samples"
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
        ok, details = grade(result.get("result_expr"))
    print(json.dumps({"check": CHECK, "pass": ok, "details": details}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
