#!/usr/bin/env python3
"""Grader for T02: is the agent's R the nondimensional 1D heat equation with a source?

Usage: grade.py <result.json>

The agent never sees this file. Its answer is untrusted text, so only the allowed names, numbers and
arithmetic are accepted before SymPy parses it (SymPy's parser evaluates what it reads).

Reference: with x = L xi, t = (rho c L^2 / k0) s and T = T0 + DT theta,
    R = theta_s - theta_xixi - Q0 L^2 / (k0 DT).
Check: with a random polynomial theta, random positive parameters and a random point (xi, s), the answer
equals a constant multiple of R, the same pure number on every sample. Because the parameters change from
sample to sample, an answer that still carries a dimensional factor (for example the equation multiplied by
k0 DT / L^2) gives a different ratio each time and fails, and so does a wrong group (L/k0 in place of L^2/k0)
or a missing or flipped term.

Prints one JSON line {"check", "pass", "details"}; exits 0 on pass, 1 on fail, 2 on an unreadable file.
"""
import hashlib
import json
import random
import re
import sys

import sympy as sp
from sympy.parsing.sympy_parser import parse_expr

CHECK = "T02_heat_nondimensionalize"
ALLOWED = {"xi", "s", "theta", "L", "k0", "rho", "c", "Q0", "T0", "DT", "Derivative"}
MAX_LEN = 2000
N_SAMPLES = 6

xi, s, L, k0, rho, c, Q0, T0, DT = sp.symbols("xi s L k0 rho c Q0 T0 DT")
theta = sp.Function("theta")
LOCAL = {"xi": xi, "s": s, "L": L, "k0": k0, "rho": rho, "c": c, "Q0": Q0, "T0": T0, "DT": DT,
         "theta": theta, "Derivative": sp.Derivative}


class Reject(Exception):
    pass


def reference():
    th = theta(xi, s)
    return sp.Derivative(th, s) - sp.Derivative(th, (xi, 2)) - Q0 * L**2 / (k0 * DT)


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


def sample(seed):
    rng = random.Random(seed)
    def r():
        return sp.Rational(rng.randint(-9, 9), rng.randint(1, 5))
    lam = sp.Lambda((xi, s), sum(r() * xi**i * s**j for i in range(4) for j in range(3)))
    nums = {L: sp.Rational(rng.randint(1, 9), rng.randint(1, 4)), k0: sp.Rational(rng.randint(1, 9), 2),
            rho: rng.randint(1, 9), c: rng.randint(1, 9), Q0: sp.Rational(rng.randint(1, 9), rng.randint(1, 3)),
            T0: rng.randint(100, 300), DT: sp.Rational(rng.randint(1, 9), 2)}
    point = {xi: sp.Rational(rng.randint(0, 8), 8), s: sp.Rational(rng.randint(0, 6), 3)}
    return lam, nums, point


def evaluate(expr, lam, nums, point):
    e = expr.replace(theta, lam).doit()
    e = sp.simplify(e.subs(nums).subs(point))
    if e.free_symbols or not e.is_number:
        raise Reject(f"does not reduce to a number for a concrete theta and parameters (left: {e})")
    return e


def grade(text):
    """Return (pass, details). Raises nothing."""
    try:
        expr = parse(text)
        if not expr.has(theta):
            raise Reject("theta does not appear")
        ref = reference()
        ratios = []
        for seed in range(N_SAMPLES):
            lam, nums, point = sample(seed)
            v_ref = evaluate(ref, lam, nums, point)
            if v_ref == 0:
                continue
            ratios.append(sp.nsimplify(evaluate(expr, lam, nums, point) / v_ref))
        if len(ratios) < 4:
            return False, "grader fault: too few usable samples"
        if ratios[0] == 0:
            return False, "the answer is zero where the reference is not"
        if any(q != ratios[0] for q in ratios[1:]):
            return False, f"not a constant multiple of the reference (ratios {[str(q) for q in ratios]})"
        return True, f"equals {ratios[0]} times the reference on {len(ratios)} random samples"
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
