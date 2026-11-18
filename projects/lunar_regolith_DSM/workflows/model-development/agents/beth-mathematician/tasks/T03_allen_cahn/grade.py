#!/usr/bin/env python3
"""Grader for T03: is the agent's R the 2D Allen-Cahn equation from the given Ginzburg-Landau energy?

Usage: grade.py <result.json>

The agent never sees this file. Its answer is untrusted text, so only the allowed names, numbers and
arithmetic are accepted before SymPy parses it (SymPy's parser evaluates what it reads).

Reference: F = int kappa/2 |grad phi|^2 + W phi^2 (1 - phi)^2, so
    mu = -kappa (phi_xx + phi_yy) + 2 W phi (1 - phi)(1 - 2 phi),
    R  = phi_t + M mu = phi_t - M kappa (phi_xx + phi_yy) + 2 M W phi (1 - phi)(1 - 2 phi).
The reference mu is also checked here against SymPy's own Euler-Lagrange operator, so the key does not
rest on hand algebra alone.
Check: with a random polynomial phi, random positive kappa, W, M and a random point, the answer equals a
constant multiple of R, the same pure number on every sample. A flipped sign, a lost factor 2 or 1/2, a
missing y-derivative or a wrong double-well derivative gives a ratio that changes from sample to sample.

Prints one JSON line {"check", "pass", "details"}; exits 0 on pass, 1 on fail, 2 on an unreadable file.
"""
import hashlib
import json
import random
import re
import sys

import sympy as sp
from sympy.calculus.euler import euler_equations
from sympy.parsing.sympy_parser import parse_expr

CHECK = "T03_allen_cahn"
ALLOWED = {"x", "y", "t", "phi", "kappa", "W", "M", "Derivative"}
MAX_LEN = 2000
N_SAMPLES = 6

x, y, t, kappa, W, M = sp.symbols("x y t kappa W M")
phi = sp.Function("phi")
LOCAL = {"x": x, "y": y, "t": t, "kappa": kappa, "W": W, "M": M, "phi": phi, "Derivative": sp.Derivative}


class Reject(Exception):
    pass


def mu_reference():
    p = phi(x, y, t)
    return (-kappa * (sp.Derivative(p, (x, 2)) + sp.Derivative(p, (y, 2)))
            + 2 * W * p * (1 - p) * (1 - 2 * p))


def mu_from_euler_lagrange():
    """SymPy's Euler-Lagrange equation for the energy density, in x and y only (t is a parameter)."""
    u = sp.Function("u")(x, y)
    density = kappa / 2 * (u.diff(x)**2 + u.diff(y)**2) + W * u**2 * (1 - u)**2
    (eq,) = euler_equations(density, u, [x, y])
    return eq.lhs.subs(u, phi(x, y, t)).doit()


def reference():
    return sp.Derivative(phi(x, y, t), t) + M * mu_reference()


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
    lam = sp.Lambda((x, y, t), sum(r() * x**i * y**j * t**m
                                   for i in range(3) for j in range(3) for m in range(2)))
    nums = {kappa: sp.Rational(rng.randint(1, 9), 2), W: sp.Rational(rng.randint(1, 9), rng.randint(1, 3)),
            M: sp.Rational(rng.randint(1, 9), rng.randint(1, 4))}
    point = {x: sp.Rational(rng.randint(0, 8), 8), y: sp.Rational(rng.randint(0, 8), 8),
             t: sp.Rational(rng.randint(0, 6), 3)}
    return lam, nums, point


def evaluate(expr, lam, nums, point):
    e = expr.replace(phi, lam).doit()
    e = sp.simplify(e.subs(nums).subs(point))
    if e.free_symbols or not e.is_number:
        raise Reject(f"does not reduce to a number for a concrete phi and parameters (left: {e})")
    return e


def grade(text):
    """Return (pass, details). Raises nothing."""
    try:
        if sp.simplify(mu_reference() - mu_from_euler_lagrange()) != 0:
            return False, "grader fault: the reference mu disagrees with SymPy's Euler-Lagrange operator"
        expr = parse(text)
        if not expr.has(phi):
            raise Reject("phi does not appear")
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
