"""Shared helpers for the M3 derivation checks of dev-toy1 (not a check by itself).

Provides
  * ledger loading (``--ledger`` argument, sha256 of what was read),
  * a small LaTeX -> SymPy parser for the notation used in the ledger statements,
  * the JSON-line reporter (one line per claim, exit 0 only if all pass),
  * an independent finite-difference solver for d_tau theta = d_xi^2 theta + c.

Nothing in this module copies a ledger statement: every check parses the statement
text it tests from the ledger file given on the command line, so a planted error in
the ledger changes what the check sees.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import sympy as sp

DEFAULT_LEDGER = (
    "/Users/jacksonbaglino/PetIGA-3.20/projects/lunar_regolith_DSM/runs/"
    "model-development/dev-toy1/derivation/ledger.json"
)
DEFAULT_SRC = "/Users/jacksonbaglino/PetIGA-3.20/projects/lunar_regolith_DSM/src"


# --------------------------------------------------------------------------- CLI
def cli(description: str, with_src: bool = False):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--ledger", default=DEFAULT_LEDGER, help="path to ledger.json")
    if with_src:
        ap.add_argument("--src", default=DEFAULT_SRC, help="src/ folder of main")
    return ap.parse_args()


class Ledger:
    def __init__(self, path):
        self.path = Path(path)
        raw = self.path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        self.data = json.loads(raw)
        self.steps = self.data["steps"]
        self.by_id = {s["id"]: s for s in self.steps}
        self.order = [s["id"] for s in self.steps]

    def stmt(self, sid):
        return self.by_id[sid]["statement"]

    def op(self, sid):
        return self.by_id[sid].get("operation", "")

    def derived_ids(self):
        return [s["id"] for s in self.steps if s.get("type") == "derived"]

    def closure(self, sid):
        seen, stack = set(), list(self.by_id[sid].get("from", []))
        while stack:
            p = stack.pop()
            if p in seen or p not in self.by_id:
                continue
            seen.add(p)
            stack.extend(self.by_id[p].get("from", []))
        return seen

    def symbol_meaning(self, latex_symbol):
        for s in self.data.get("symbols", []):
            if s["symbol"] == latex_symbol:
                return s["meaning"]
        return None


class Reporter:
    """Prints one JSON line per claim; exit code 0 only if every claim passes."""

    def __init__(self, script, ledger: Ledger, extra_inputs=None):
        self.script = script
        self.results = []
        used = {"script": script, "ledger": str(ledger.path), "ledger_sha256": ledger.sha256}
        if extra_inputs:
            used.update(extra_inputs)
        print("inputs used: " + json.dumps(used), file=sys.stderr)

    def check(self, check, claim, fn):
        try:
            ok, details = fn()
        except Exception as exc:  # noqa: BLE001 - a crash is a failed check
            ok, details = False, f"check raised {type(exc).__name__}: {exc}"
        rec = {"check": check, "claim": claim, "pass": bool(ok), "details": details}
        self.results.append(rec)
        print(json.dumps(rec, default=str))
        sys.stdout.flush()
        return ok

    def finish(self):
        n_fail = sum(1 for r in self.results if not r["pass"])
        print(f"{self.script}: {len(self.results) - n_fail} pass, {n_fail} fail", file=sys.stderr)
        sys.exit(0 if n_fail == 0 and self.results else 1)


# ---------------------------------------------------------------- symbols/fields
POSITIVE = {"L", "k", "rho", "c_p", "T_0", "kappa", "x_c", "t_c"}
FIELD_VARS = {
    "T": ("x", "t"), "e": ("x", "t"), "j": ("x", "t"), "w": ("x", "t"),
    "T_1": ("x", "t"), "T_2": ("x", "t"),
    "theta": ("xi", "tau"), "thetap": ("xi", "tau"), "v": ("xi", "tau"),
    "theta_s": ("xi",), "T_s": ("x",), "E": ("t",), "g": ("tau",),
}
_SYMS: dict = {}


def S(name):
    if name == "pi":
        return sp.pi
    if name == "oo":
        return sp.oo
    if name not in _SYMS:
        if name in POSITIVE:
            _SYMS[name] = sp.Symbol(name, positive=True)
        elif name == "r":
            _SYMS[name] = sp.Symbol(name, integer=True, nonnegative=True)
        elif name == "ell":
            _SYMS[name] = sp.Symbol(name, integer=True, positive=True)
        else:
            _SYMS[name] = sp.Symbol(name, real=True)
    return _SYMS[name]


def Fn(name):
    return sp.Function(name)


def field(name):
    return Fn(name)(*[S(v) for v in FIELD_VARS[name]])


def plug(expr, name, lam):
    """Replace field ``name`` by the Lambda ``lam`` and evaluate derivatives/integrals."""
    return expr.replace(Fn(name), lam).doit()


# ------------------------------------------------------------------- tokenizer
class ParseError(Exception):
    pass


FUNCS = {"sin": sp.sin, "cos": sp.cos, "exp": sp.exp, "ln": sp.log, "log": sp.log, "tan": sp.tan}
RELCMD = {"le": "<=", "leq": "<=", "ge": ">=", "geq": ">=", "ne": "!=", "neq": "!=", "equiv": "="}
SEPCMD = {"Rightarrow", "Leftrightarrow", "implies", "SEP", "in", "to", "mapsto", "forall"}
LIMCMD = {"int", "sum", "sumodd", "max", "min", "lim"}
RELOPS = {"=", "<=", ">=", "<", ">", "!="}


def _balanced(s, i):
    assert s[i] == "{"
    depth, j = 0, i
    while j < len(s):
        if s[j] == "{" and (j == 0 or s[j - 1] != "\\"):
            depth += 1
        elif s[j] == "}" and (j == 0 or s[j - 1] != "\\"):
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    raise ParseError("unbalanced braces")


def strip_group_cmd(s, cmd, repl):
    out, i, pat = [], 0, "\\" + cmd
    while True:
        j = s.find(pat, i)
        if j < 0:
            out.append(s[i:])
            break
        k = j + len(pat)
        if k < len(s) and s[k].isalpha():
            out.append(s[i:k])
            i = k
            continue
        out.append(s[i:j])
        while k < len(s) and s[k] == " ":
            k += 1
        if k < len(s) and s[k] == "{":
            content, k = _balanced(s, k)
        else:
            content = ""
        out.append(repl(content) if callable(repl) else repl)
        i = k
    return "".join(out)


def preprocess(s):
    sep = " \\SEP "
    s = re.sub(r"=\s*\\begin\{cases\}.*?\\end\{cases\}", lambda m: " ", s, flags=re.S)
    s = re.sub(r"\\begin\{cases\}.*?\\end\{cases\}", lambda m: sep, s, flags=re.S)
    s = re.sub(r"\\sum_\{\\ell\s*\\text\{\s*odd\s*\}\}", lambda m: "\\sumodd ", s)
    s = strip_group_cmd(s, "text", sep)
    s = strip_group_cmd(s, "mathrm", " ")
    s = strip_group_cmd(s, "mathbb", " ")
    s = re.sub(r"\\(left|right)\.", lambda m: "", s)
    s = re.sub(r"\\(left|right|Bigg|bigg|Big|big|displaystyle)(?![A-Za-z])", lambda m: "", s)
    s = re.sub(r"\\(qquad|quad)(?![A-Za-z])", lambda m: sep, s)
    s = re.sub(r"\\[,;:! ]", lambda m: " ", s)
    s = re.sub(r"\\frac\s*(\d)\s*(\d)", lambda m: "\\frac{%s}{%s}" % (m.group(1), m.group(2)), s)
    return s


def _sanitize(raw):
    raw = raw.replace("\\", "").replace(" ", "").replace("{", "").replace("}", "")
    return raw.replace("*", "star")


def _read_script(s, i):
    """s[i] is '_' or '^'; return (raw, next index)."""
    i += 1
    while i < len(s) and s[i] == " ":
        i += 1
    if i >= len(s):
        raise ParseError("dangling script")
    if s[i] == "{":
        return _balanced(s, i)
    if s[i] == "\\":
        m = re.match(r"\\[A-Za-z]+", s[i:])
        return m.group(0), i + len(m.group(0))
    return s[i], i + 1


def tokenize(s):
    s = preprocess(s)
    toks, i, n = [], 0, len(s)
    while i < n:
        c = s[i]
        if c.isspace():
            i += 1
            continue
        if c == "\\":
            m = re.match(r"\\([A-Za-z]+)", s[i:])
            if not m:
                raise ParseError(f"bad escape at {s[i:i+5]!r}")
            name = m.group(1)
            i += len(m.group(0))
            if name in RELCMD:
                toks.append(("op", RELCMD[name]))
            elif name in ("cdot", "times"):
                toks.append(("op", "*"))
            elif name in SEPCMD:
                toks.append(("sep",))
            elif name == "infty":
                toks.append(("name", "oo"))
            elif name == "frac":
                toks.append(("frac",))
            elif name in FUNCS:
                toks.append(("func", name))
            elif name == "partial":
                sub = None
                if i < n and s[i] == "_":
                    raw, i = _read_script(s, i)
                    sub = _sanitize(raw)
                toks.append(("partial", sub))
            elif name in LIMCMD:
                sub = sup = None
                for _ in range(2):
                    if i < n and s[i] == "_":
                        sub, i = _read_script(s, i)
                    elif i < n and s[i] == "^":
                        sup, i = _read_script(s, i)
                toks.append(("lim", name, sub, sup))
            elif name == "Delta":
                j = i
                while j < n and s[j] == " ":
                    j += 1
                if j < n and s[j].isalpha():
                    base, i = s[j], j + 1
                    if i < n and s[i] == "_":
                        raw, i = _read_script(s, i)
                        base += "_" + _sanitize(raw)
                    toks.append(("name", "Delta" + base))
                else:
                    toks.append(("name", "Delta"))
            elif name.isalpha() and name not in ("begin", "end", "lim"):
                base = name
                if i < n and s[i] == "_":
                    raw, i = _read_script(s, i)
                    base += "_" + _sanitize(raw)
                toks.append(("name", base))
            else:
                toks.append(("bad", name))
            continue
        if c.isdigit() or (c == "." and i + 1 < n and s[i + 1].isdigit()):
            m = re.match(r"\d*\.?\d+", s[i:])
            toks.append(("num", m.group(0)))
            i += len(m.group(0))
            continue
        if c.isalpha():
            base = c
            i += 1
            if i < n and s[i] == "_":
                raw, i = _read_script(s, i)
                base += "_" + _sanitize(raw)
            toks.append(("name", base))
            continue
        if c in "{}()[]^+-=<>,:;|*/'!&":
            toks.append(("op", c))
            i += 1
            continue
        raise ParseError(f"unexpected character {c!r}")
    return toks


def _is_op(tok, *vals):
    return tok is not None and tok[0] == "op" and tok[1] in vals


def tok_text(toks):
    out = []
    for t in toks:
        if t[0] in ("op", "name", "num", "func"):
            out.append(t[1])
        elif t[0] == "partial":
            out.append("d_" + (t[1] or ""))
        elif t[0] == "lim":
            out.append(t[1])
        else:
            out.append(t[0])
    return " ".join(out)


# ---------------------------------------------------------------------- parser
class _P:
    def __init__(self, toks):
        self.t, self.i = toks, 0

    def peek(self, k=0):
        j = self.i + k
        return self.t[j] if j < len(self.t) else None

    def next(self):
        tok = self.peek()
        self.i += 1
        return tok

    def expect(self, val):
        tok = self.next()
        if not _is_op(tok, val):
            raise ParseError(f"expected {val!r}, got {tok}")

    def at_end(self):
        return self.i >= len(self.t)

    def full_expr(self):
        v = self.expr()
        if not self.at_end():
            raise ParseError(f"trailing tokens: {tok_text(self.t[self.i:])}")
        return v

    def full_chain(self):
        exprs, rels = [self.expr()], []
        while _is_op(self.peek(), *RELOPS):
            rels.append(self.next()[1])
            exprs.append(self.expr())
        if not self.at_end():
            raise ParseError(f"trailing tokens: {tok_text(self.t[self.i:])}")
        return exprs, rels

    def expr(self):
        sign = 1
        if _is_op(self.peek(), "+"):
            self.next()
        elif _is_op(self.peek(), "-"):
            self.next()
            sign = -1
        val = sign * self.term()
        while True:
            if _is_op(self.peek(), "+"):
                self.next()
                val = val + self.term()
            elif _is_op(self.peek(), "-"):
                self.next()
                val = val - self.term()
            else:
                return val

    @staticmethod
    def starts_factor(tok):
        if tok is None:
            return False
        if tok[0] in ("name", "num", "func", "frac", "partial", "lim"):
            return True
        return _is_op(tok, "(", "{")

    def term(self):
        val = self.factor()
        while True:
            tok = self.peek()
            if _is_op(tok, "*"):
                self.next()
                val = val * self.factor()
            elif _is_op(tok, "/"):
                self.next()
                val = val / self.factor()
            elif self.starts_factor(tok):
                val = val * self.factor()
            else:
                return val

    def factor(self):
        base, vars_ = self.primary()
        base, vars_ = self.postfix(base, vars_)
        if _is_op(self.peek(), "^"):
            self.next()
            if _is_op(self.peek(), "-"):
                self.next()
                expo = -self.primary()[0]
            else:
                expo = self.primary()[0]
            base = base ** expo
        return base

    def postfix(self, base, vars_):
        if vars_ is not None and _is_op(self.peek(), "("):
            self.next()
            args = [self.expr()]
            while _is_op(self.peek(), ","):
                self.next()
                args.append(self.expr())
            self.expect(")")
            if len(args) != len(vars_):
                raise ParseError("wrong number of arguments for field")
            base = base.subs({S(v): a for v, a in zip(vars_, args)}, simultaneous=True)
            vars_ = None
        return base, vars_

    def group_tokens(self):
        tok = self.next()
        if _is_op(tok, "{"):
            depth, start = 1, self.i
            while depth:
                t = self.next()
                if t is None:
                    raise ParseError("unbalanced group")
                if _is_op(t, "{"):
                    depth += 1
                elif _is_op(t, "}"):
                    depth -= 1
            return self.t[start:self.i - 1]
        if tok is None:
            raise ParseError("missing group")
        return [tok]

    def frac(self):
        num, den = self.group_tokens(), self.group_tokens()
        dnum = bool(num) and (num[0][0] == "partial" or num[0] == ("name", "d"))
        dden = bool(den) and (den[0][0] == "partial" or den[0] == ("name", "d"))
        if dnum and dden:
            j, order = 1, 1
            if len(num) > 2 and _is_op(num[1], "^") and num[2][0] == "num":
                order, j = int(num[2][1]), 3
            if len(den) < 2 or den[1][0] != "name":
                raise ParseError("bad derivative denominator")
            var = S(den[1][1])
            dorder = 1
            if len(den) == 4 and _is_op(den[2], "^") and den[3][0] == "num":
                dorder = int(den[3][1])
            elif len(den) != 2:
                raise ParseError("bad derivative denominator")
            if dorder != order:
                raise ParseError("derivative order mismatch between numerator and denominator")
            operand_toks = num[j:]
            if operand_toks:
                sub = _P(operand_toks)
                op, vars_ = sub.primary()
                if not sub.at_end():
                    raise ParseError("derivative operand is not a single primary")
            else:
                op, vars_ = self.factor(), None
            return sp.Derivative(op, (var, order)), vars_
        n = _P(num).full_expr()
        d = _P(den).full_expr()
        return n / d, None

    def primary(self):
        tok = self.next()
        if tok is None:
            raise ParseError("unexpected end of expression")
        kind = tok[0]
        if kind == "num":
            v = tok[1]
            return (sp.Integer(v) if v.isdigit() else sp.Rational(v)), None
        if kind == "name":
            n = tok[1]
            if n == "d":
                raise ParseError("stray differential d")
            primes = 0
            while _is_op(self.peek(), "'"):
                self.next()
                primes += 1
            if n == "theta" and primes >= 1:
                n, primes = "thetap", primes - 1
            if n in FIELD_VARS:
                f = field(n)
                vars_ = FIELD_VARS[n]
                if primes:
                    if len(vars_) != 1:
                        raise ParseError("prime on a field of two variables")
                    f = sp.Derivative(f, (S(vars_[0]), primes))
                return f, vars_
            if primes:
                raise ParseError(f"prime on non-field {n}")
            return S(n), None
        if kind == "op" and tok[1] in "({[":
            close = {"(": ")", "{": "}", "[": "]"}[tok[1]]
            v = self.expr()
            self.expect(close)
            return v, None
        if kind == "op" and tok[1] == "|":
            v = self.expr()
            self.expect("|")
            return sp.Abs(v), None
        if kind == "func":
            arg = self.factor()
            return FUNCS[tok[1]](arg), None
        if kind == "frac":
            return self.frac()
        if kind == "partial":
            if tok[1] is None:
                raise ParseError("bare \\partial outside a fraction")
            order = 1
            if _is_op(self.peek(), "^"):
                self.next()
                order = int(self.next()[1])
            op, vars_ = self.primary()
            return sp.Derivative(op, (S(tok[1]), order)), vars_
        if kind == "lim":
            _, which, sub, sup = tok
            if which == "int":
                lo, hi = parse_expr_str(sub), parse_expr_str(sup)
                depth, j = 0, self.i
                while j < len(self.t):
                    tt = self.t[j]
                    if _is_op(tt, "(", "{", "["):
                        depth += 1
                    elif _is_op(tt, ")", "}", "]"):
                        if depth == 0:
                            raise ParseError("integral without differential")
                        depth -= 1
                    elif depth == 0 and tt == ("name", "d") and j + 1 < len(self.t) \
                            and self.t[j + 1][0] == "name":
                        break
                    j += 1
                else:
                    raise ParseError("integral without differential")
                body = _P(self.t[self.i:j]).full_expr()
                var = S(self.t[j + 1][1])
                self.i = j + 2
                return sp.Integral(body, (var, lo, hi)), None
            if which == "sumodd":
                body = self.term()
                return sp.Sum(body.subs(S("ell"), 2 * S("r") + 1), (S("r"), 0, sp.oo)), None
            if which in ("max", "min"):
                op = self.factor()
                tag = ("Max_" if which == "max" else "Min_") + _sanitize(sub or "")
                return Fn(tag)(op), None
            raise ParseError(f"unsupported \\{which}")
        if kind == "bad":
            raise ParseError(f"unsupported command \\{tok[1]}")
        raise ParseError(f"unexpected token {tok}")


def parse_expr_str(s):
    if s is None:
        raise ParseError("missing limit")
    return _P(tokenize(s)).full_expr()


def _match_paren(toks, i):
    depth = 0
    for j in range(i, len(toks)):
        if _is_op(toks[j], "(", "{", "["):
            depth += 1
        elif _is_op(toks[j], ")", "}", "]"):
            depth -= 1
            if depth == 0:
                return j
    raise ParseError("unbalanced parentheses")


def _has_top_rel(toks):
    depth = 0
    for t in toks:
        if _is_op(t, "(", "{", "["):
            depth += 1
        elif _is_op(t, ")", "}", "]"):
            depth -= 1
        elif depth == 0 and _is_op(t, *RELOPS):
            return True
    return False


def split_pieces(toks):
    pieces, cur, depth, i = [], [], 0, 0
    while i < len(toks):
        tk = toks[i]
        if _is_op(tk, "(") and depth == 0:
            j = _match_paren(toks, i)
            if _has_top_rel(toks[i + 1:j]):  # "(t>0)" style condition: a separator
                pieces.append(cur)
                cur = []
                i = j + 1
                continue
        if _is_op(tk, "(", "{", "["):
            depth += 1
        elif _is_op(tk, ")", "}", "]"):
            depth -= 1
        if depth == 0 and (tk[0] == "sep" or _is_op(tk, ",", ":", ";")):
            pieces.append(cur)
            cur = []
            i += 1
            continue
        cur.append(tk)
        i += 1
    pieces.append(cur)
    return [p for p in pieces if p]


class Piece:
    def __init__(self, toks):
        self.text = tok_text(toks)
        self.has_rel = _has_top_rel(toks)
        self.ok, self.exprs, self.rels, self.error = False, [], [], None
        try:
            self.exprs, self.rels = _P(toks).full_chain()
            self.ok = True
        except (ParseError, Exception) as exc:  # noqa: BLE001
            self.error = f"{type(exc).__name__}: {exc}"

    @property
    def lhs(self):
        return self.exprs[0]

    @property
    def rhs(self):
        return self.exprs[-1]

    def residual(self):
        if not (self.ok and len(self.exprs) == 2 and self.rels == ["="]):
            raise ParseError(f"not a single equation: {self.text} ({self.error})")
        return sp.expand(self.exprs[0] - self.exprs[1])

    def pairs(self):
        """Equalities a = b from a chain a = b = c (each against the last)."""
        if not self.ok or any(r != "=" for r in self.rels):
            raise ParseError(f"not an equality chain: {self.text} ({self.error})")
        return [(e, self.exprs[-1]) for e in self.exprs[:-1]]


class _FailedPiece:
    def __init__(self, stmt, err):
        self.text, self.has_rel, self.ok = stmt[:80], True, False
        self.exprs, self.rels, self.error = [], [], f"tokenize: {err}"


def pieces(stmt, only_relations=True):
    try:
        toks = tokenize(stmt)
    except (ParseError, Exception) as exc:  # noqa: BLE001
        return [_FailedPiece(stmt, exc)]
    out = [Piece(p) for p in split_pieces(toks)]
    return [p for p in out if p.has_rel] if only_relations else out


def piece(stmt, idx):
    ps = pieces(stmt)
    if idx >= len(ps):
        raise ParseError(f"statement has {len(ps)} relation pieces, wanted index {idx}")
    p = ps[idx]
    if not p.ok:
        raise ParseError(f"piece {idx} did not parse: {p.text} ({p.error})")
    return p


def find_piece(stmt, pred, what):
    for p in pieces(stmt):
        if p.ok and pred(p):
            return p
    bad = [f"{p.text} ({p.error})" for p in pieces(stmt) if not p.ok]
    raise ParseError(f"no piece with {what}; unparsed pieces: {bad}")


def lhs_is(name):
    """Predicate: piece is an equation whose left side is symbol/field ``name``."""
    def pred(p):
        if len(p.exprs) < 2 or p.rels[0] != "=":
            return False
        lhs = p.exprs[0]
        if name in FIELD_VARS:
            return isinstance(lhs, sp.core.function.AppliedUndef) and lhs.func.__name__ == name
        return lhs == S(name)
    return pred


def same_up_to_sign(a, b):
    d1 = sp.simplify(sp.expand(a - b))
    if d1 == 0:
        return True, "equal"
    d2 = sp.simplify(sp.expand(a + b))
    if d2 == 0:
        return True, "equal up to overall sign"
    return False, f"difference {sp.sstr(d1)}"


def short(x, n=300):
    s = sp.sstr(x) if isinstance(x, sp.Basic) else str(x)
    return s if len(s) <= n else s[:n] + "..."


# ------------------------------------------------- independent numerical solver
def solve_theta_fd(source, taus, n_int=400):
    """Method of lines for d_tau th = d_xi^2 th + source on (0,1), th=0 at the ends
    and at tau=0. Second-order centered differences, stiff BDF integrator.
    Returns (xi grid incl. ends, array [len(taus), n_int+2])."""
    import numpy as np
    from scipy.integrate import solve_ivp
    from scipy.sparse import diags

    h = 1.0 / (n_int + 1)
    A = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n_int, n_int)) / h**2
    A = A.tocsr()

    def rhs(_t, y):
        return A @ y + source

    sol = solve_ivp(rhs, (0.0, max(taus)), np.zeros(n_int), method="BDF", t_eval=sorted(taus),
                    jac=A, rtol=1e-10, atol=1e-13)
    xi = np.linspace(0.0, 1.0, n_int + 2)
    out = np.zeros((len(sol.t), n_int + 2))
    out[:, 1:-1] = sol.y.T
    return xi, sol.t, out


# ------------------------------------------------------- structure and coverage
def _covers_of_scripts(folder):
    out = {}
    for p in sorted(Path(folder).glob("m3_*.py")):
        m = re.search(r"^COVERS\s*=\s*\[(.*?)\]", p.read_text(), re.S | re.M)
        if m:
            out[p.name] = re.findall(r"'(S\d+)'|\"(S\d+)\"", m.group(1))
            out[p.name] = [a or b for a, b in out[p.name]]
    return out


def main():
    ap = argparse.ArgumentParser(description="M3 structure, parse and coverage checks")
    ap.add_argument("--ledger", default=DEFAULT_LEDGER)
    ap.add_argument("--dump", default=None, help="write every parsed piece to this file")
    args = ap.parse_args()
    L = Ledger(args.ledger)
    rep = Reporter("m3_common.py", L, {"checks_dir": str(Path(__file__).resolve().parent)})

    def structure():
        bad = {}
        pos = {sid: i for i, sid in enumerate(L.order)}
        for sid in L.derived_ids():
            st = L.by_id[sid]
            probs = []
            if not st.get("operation"):
                probs.append("no operation")
            for p in st.get("from", []):
                if p not in L.by_id:
                    probs.append(f"premise {p} does not exist")
                elif pos[p] >= pos[sid]:
                    probs.append(f"premise {p} comes after {sid} in derivation order")
            if not st.get("from"):
                probs.append("empty from list")
            if probs:
                bad[sid] = probs
        dup = len(L.order) != len(set(L.order))
        return not bad and not dup, {"problems": bad, "duplicate ids": dup,
                                     "n_derived": len(L.derived_ids())}
    rep.check("structure", "every derived step names an operation and cites existing earlier steps",
              structure)

    def parse_all():
        failed, n_ok, lines = {}, 0, []
        for sid in L.order:
            for p in pieces(L.stmt(sid)):
                lines.append(f"{sid} | {p.text} | {[str(e) for e in p.exprs]} {p.rels} {p.error}")
                if p.ok:
                    n_ok += 1
                else:
                    failed.setdefault(sid, []).append(f"{p.text} ({p.error})")
        if args.dump:
            Path(args.dump).write_text("\n".join(lines) + "\n")
        # pieces that are text or set conditions are allowed to fail; report them
        return True, {"relation pieces parsed": n_ok, "not parsed (informative)": failed}
    rep.check("parse_report", "statements parse into SymPy relations (informative)", parse_all)

    def coverage():
        cov = _covers_of_scripts(Path(__file__).resolve().parent)
        where = {sid: [s for s, ids in cov.items() if sid in ids] for sid in L.derived_ids()}
        missing = [sid for sid, w in where.items() if not w]
        return not missing, {"missing": missing, "n_derived": len(where), "covered_by": where}
    rep.check("coverage", "every derived step is listed in COVERS of at least one m3 check", coverage)
    rep.finish()


if __name__ == "__main__":
    main()
