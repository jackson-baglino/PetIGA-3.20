"""Shared helpers for the M4 discretization checks of dev-toy1 (not a check by itself).

Every check reads the text it tests from the files given on the command line:
  --spec    derivation/spec.md       (display equations M4-1 ... M4-9, tables, text formulas)
  --ledger  derivation/ledger.json   (strong form S29 and the scales)
so a planted error in either file changes what the checks see.

The LaTeX of the spec is turned into SymPy with the M3 parser (m3_common), after a
spec-specific preprocessing step (U^{n+1} -> Unew, \\dot U -> Udot, \\Delta\\tau -> dtau, ...).
Plain-text (Unicode) formulas in tables are turned into SymPy by ``uni``.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

RUN = "/Users/jacksonbaglino/PetIGA-3.20/projects/lunar_regolith_DSM/runs/model-development/dev-toy1"
DEFAULT_SPEC = RUN + "/derivation/spec.md"
DEFAULT_LEDGER = RUN + "/derivation/ledger.json"
DEFAULT_MMS = RUN + "/derivation/mms/mms_sources.py"
DEFAULT_M1 = RUN + "/checks/m1_inputs.json"

# fields used by the spec equations
C.FIELD_VARS["f"] = ("xi", "tau")
C.FIELD_VARS["w"] = ("xi",)
C.FIELD_VARS["thetam"] = ("xi", "tau")

Reporter = C.Reporter
Ledger = C.Ledger


def cli(description, extra=()):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--spec", default=DEFAULT_SPEC, help="path to spec.md")
    ap.add_argument("--ledger", default=DEFAULT_LEDGER, help="path to ledger.json")
    for name, default in extra:
        ap.add_argument(name, default=default)
    return ap.parse_args()


class Spec:
    def __init__(self, path):
        self.path = Path(path)
        raw = self.path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        self.text = raw.decode("utf-8")
        self.blocks = {}
        prev_untagged = []
        for m in re.finditer(r"\$\$(.*?)\$\$", self.text, re.S):
            body = m.group(1)
            t = re.search(r"\\tag\{([^}]*)\}", body)
            if t:
                tag = t.group(1)
                self.blocks[tag] = body
                for i, b in enumerate(prev_untagged):
                    self.blocks[f"{tag}_pre{i if i else ''}"] = b
                prev_untagged = []
            else:
                prev_untagged.append(body)

    def block(self, tag):
        if tag not in self.blocks:
            raise KeyError(f"display equation {tag} not found in spec")
        return self.blocks[tag]

    def section(self, title_regex):
        m = re.search(r"^(#+)\s*" + title_regex + r".*?$", self.text, re.M)
        if not m:
            raise KeyError(f"section {title_regex} not found")
        level = len(m.group(1))
        rest = self.text[m.end():]
        n = re.search(r"^#{1,%d}\s" % level, rest, re.M)
        return rest[: n.start()] if n else rest

    def table_row(self, first_cell_regex, within=None):
        txt = within if within is not None else self.text
        for line in txt.splitlines():
            if line.startswith("|"):
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                if cells and re.fullmatch(first_cell_regex, cells[0]):
                    return cells
        raise KeyError(f"table row {first_cell_regex!r} not found")


# ----------------------------------------------------------- LaTeX preprocessing
def _balanced(s, i):
    return C._balanced(s, i)


def strip_underbrace(s):
    out, i = [], 0
    while True:
        j = s.find("\\underbrace", i)
        if j < 0:
            out.append(s[i:])
            return "".join(out)
        out.append(s[i:j])
        k = j + len("\\underbrace")
        while s[k] == " ":
            k += 1
        content, k = _balanced(s, k)
        while k < len(s) and s[k] == " ":
            k += 1
        if k < len(s) and s[k] == "_":
            k += 1
            while s[k] == " ":
                k += 1
            if s[k] == "{":
                _, k = _balanced(s, k)
        out.append(" +(" + content + ") ")
        i = k


IDX = r"(\{[^{}]*\}|[A-Za-z0-9])"


def latex_pre(s):
    s = re.sub(r"\\tag\{[^}]*\}", " ", s)
    s = strip_underbrace(s)
    s = s.replace("\\theta'_m", "\\thetam ")
    s = re.sub(r"H\^1_0\(0,\s*1\)", " ", s)
    s = re.sub(r"\\dot\s*U_", r"\\Udot_", s)
    s = re.sub(r"\\dot\s*U(?![A-Za-z_])", r"\\Udot ", s)
    s = re.sub(r"(?<![A-Za-z\\])U_" + IDX + r"\^\{n\+1\}", r"\\Unew_\1", s)
    s = re.sub(r"(?<![A-Za-z\\])U_" + IDX + r"\^n", r"\\Uold_\1", s)
    s = s.replace("U^{n+1}", "\\Unew ").replace("U^n", "\\Uold ")
    s = s.replace("\\tau^{n+1}", "\\taunew ")
    s = re.sub(r"\\Delta\s*\\tau", r"\\dtau", s)
    s = re.sub(r"(?<![A-Za-z\\])e\^\{", r"\\exp{", s)
    s = re.sub(r"\\sin\^2\\?!?\s*\\frac\{([^{}]*)\}\{([^{}]*)\}", r"(\\sin(\\frac{\1}{\2}))^{2}", s)
    s = re.sub(r"\\varepsilon_\{\\theta'\}", r"\\epsth ", s)
    s = re.sub(r"(?<![A-Za-z\\])dt(?![A-Za-z])", r"\\dtstep ", s)
    s = s.replace("\\lceil", "(").replace("\\rceil", ")")
    s = s.replace("\\dots", " ")
    s = re.sub(r"\\(big|Big|bigg|Bigg)\[", "(", s)
    s = re.sub(r"\\(big|Big|bigg|Bigg)\]", ")", s)
    return s


def parse_latex_expr(s):
    return C._P(C.tokenize(latex_pre(s))).full_expr()


def latex_pieces(s, only_relations=True):
    return C.pieces(latex_pre(s), only_relations=only_relations)


def piece_eq(s, pred, what):
    for p in latex_pieces(s):
        if p.ok and pred(p):
            return p
    bad = [f"{p.text} ({p.error})" for p in latex_pieces(s) if not p.ok]
    raise C.ParseError(f"no piece with {what}; unparsed: {bad}")


# ------------------------------------------------------------ index normalising
_IDXRE = re.compile(r"^(U|Udot|Unew|Uold|xi)_([A-Za-z])([+-]1)?$")


def normalise_row(expr, idx=None):
    """Map U_{j-1} -> U_m, U_j -> U_c, U_{j+1} -> U_p (same for Udot/Unew/Uold/xi).
    Returns (expr, index letter used)."""
    rep, letters = {}, set()
    for s in expr.free_symbols:
        m = _IDXRE.match(s.name)
        if m:
            base, letter, off = m.groups()
            letters.add(letter)
            tag = {None: "c", "-1": "m", "+1": "p"}[off]
            rep[s] = sp.Symbol(f"{base}_{tag}", real=True)
    if idx is not None and letters - {idx}:
        raise ValueError(f"unexpected row index letters {letters}")
    return expr.xreplace(rep), letters


def sym(name):
    return sp.Symbol(name, real=True)


# -------------------------------------------------------------- unicode formulas
_UNI = [
    ("−", "-"), ("–", "-"), ("·", "*"), ("×", "*"), ("²", "**2"), ("³", "**3"), ("⁻¹", "**(-1)"),
    ("Δτ", " dtau "), ("τ", " tau "), ("ξ", " xi "), ("κ", " kappa "), ("ρ", " rho "),
    ("π", " pi "), ("λ", " lam "), ("θ′", " thp "), ("′", "p"), ("ε", " eps "),
]
_NAMES = ["shift", "dtau", "kappa", "rho", "cp", "c_p", "T0", "T_0", "pi", "lam", "eps", "tau",
          "xi", "dt", "g_L", "g_R", "U_0", "U_1", "U_N", "L", "k", "q", "h", "N", "Q", "e"]


def uni(s, extra_names=()):
    """Unicode text formula -> SymPy, with implicit multiplication between known names."""
    for a, b in _UNI:
        s = s.replace(a, b)
    s = re.sub(r"\s+_", "_", s)
    s = s.replace("^", "**")
    names = sorted(set(_NAMES) | set(extra_names), key=len, reverse=True)
    toks, i = [], 0
    while i < len(s):
        c = s[i]
        if c.isspace():
            i += 1
            continue
        if c.isdigit() or c == ".":
            m = re.match(r"\d*\.?\d+(e-?\d+)?", s[i:])
            toks.append(("n", m.group(0)))
            i += len(m.group(0))
            continue
        if c.isalpha() or c == "_":
            for nm in names:
                if s.startswith(nm, i) and not (i + len(nm) < len(s) and (s[i + len(nm)] == "_")):
                    toks.append(("v", nm))
                    i += len(nm)
                    break
            else:
                raise ValueError(f"unknown name at {s[i:i+12]!r} in {s!r}")
            continue
        if s.startswith("**", i):
            toks.append(("o", "**"))
            i += 2
            continue
        toks.append(("o", c))
        i += 1
    out = []
    for k, (kind, val) in enumerate(toks):
        if out:
            pk, pv = toks[k - 1]
            left = pk in ("n", "v") or pv == ")"
            right = kind in ("n", "v") or val == "("
            if left and right:
                out.append("*")
        out.append(val)
    expr = "".join(out)
    loc = {nm: sp.Symbol(nm.replace("c_p", "cp").replace("T_0", "T0"), real=True) for nm in names}
    loc["pi"] = sp.pi
    loc["e"] = sp.E
    return sp.sympify(expr, locals=loc)


def U(name):
    return sp.Symbol(name, real=True)


def num_in(text, pattern, what):
    m = re.search(pattern, text, re.S)
    if not m:
        raise KeyError(f"{what}: pattern {pattern!r} not found in spec")
    return float(m.group(1).replace("−", "-").replace("×10⁻", "e-"))


def rel_ok(a, b, tol):
    return abs(a - b) <= tol * max(abs(b), 1e-300)
