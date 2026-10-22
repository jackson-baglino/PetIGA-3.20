"""Parse units, printed values and symbols so that S3 compares meaning, not strings.

- parse_unit("J/(m3 K)") -> Unit(scale=1.0, dims=(-1, 1, -2, -1, 0)) over (m, kg, s, K, mol).
  "J m^-3 K^-1", "J m⁻³ K⁻¹", r"\\mathrm{J\\,m^{-3}\\,K^{-1}}" and "J/(m3 K)" all give the same Unit.
  °C is read as a temperature difference, so W/(m °C) == W m^-1 K^-1 (WORKFLOW.md S3).
  After "/", every juxtaposed factor up to the next "/" or ")" is in the denominator: W/m K == W/(m K).
- parse_numbers("1.8 × 10⁶, 1.4×10^{3}") -> [1.8e6, 1400.0]
- flatten(latex) turns Marker's LaTeX into plain text (Greek letters as Unicode, ^ and _ kept).
- symbol_regex(printed) finds a printed symbol (e.g. "ρi", "vn,j") in flattened text.
"""
import math
import re
from typing import NamedTuple

DIM_NAMES = ("m", "kg", "s", "K", "mol")


def _d(m=0, kg=0, s=0, K=0, mol=0):
    return (m, kg, s, K, mol)


# Base and derived units. "h", "d", "L", "C" and "min" are left out on purpose: in these papers they are
# symbols (grid spacing h, pore size L, heat capacity C), not units.
UNITS = {
    "m": (1.0, _d(m=1)), "g": (1e-3, _d(kg=1)), "s": (1.0, _d(s=1)), "K": (1.0, _d(K=1)),
    "degC": (1.0, _d(K=1)),  # temperature difference; an absolute °C value differs by an offset, not compared
    "mol": (1.0, _d(mol=1)), "Hz": (1.0, _d(s=-1)),
    "N": (1.0, _d(m=1, kg=1, s=-2)), "Pa": (1.0, _d(m=-1, kg=1, s=-2)),
    "J": (1.0, _d(m=2, kg=1, s=-2)), "W": (1.0, _d(m=2, kg=1, s=-3)),
}
ALIASES = {"Kg": "kg", "KG": "kg", "sec": "s"}
# Micro is accepted only as μ or µ: an ASCII "u" prefix would make the word "us" a unit.
PREFIXES = {"G": 1e9, "M": 1e6, "k": 1e3, "c": 1e-2, "m": 1e-3, "μ": 1e-6, "n": 1e-9, "p": 1e-12}
PROSE_ONLY_IN_PARSE = {"g"}  # bare "g" is a unit to parse_unit, but in running text it is usually "e.g."
DIMENSIONLESS = {"", "1", "-", "--", "dimensionless", "none", "unitless", "(-)", "(1)"}


class UnitError(ValueError):
    pass


class Unit(NamedTuple):
    scale: float
    dims: tuple

    def __mul__(self, o):
        return Unit(self.scale * o.scale, tuple(a + b for a, b in zip(self.dims, o.dims)))

    def __pow__(self, n):
        return Unit(self.scale ** n, tuple(a * n for a in self.dims))

    def inverse(self):
        return self ** -1

    def canonical(self):
        """Readable SI form, e.g. 'kg m^-1 s^-2', with a leading scale if not 1."""
        parts = [n if e == 1 else f"{n}^{e}" for n, e in zip(DIM_NAMES, self.dims) if e]
        s = " ".join(sorted(parts, key=lambda p: -int(p.split("^")[1]) if "^" in p else -1)) or "1"
        return s if math.isclose(self.scale, 1.0, rel_tol=1e-12) else f"{self.scale:g} {s}"


ONE = Unit(1.0, _d())


def same_unit(a, b, rel=1e-9):
    return a.dims == b.dims and math.isclose(a.scale, b.scale, rel_tol=rel)


# ---------------------------------------------------------------- LaTeX and Unicode normalisation

GREEK = {
    "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "varepsilon": "ε", "zeta": "ζ",
    "eta": "η", "theta": "θ", "vartheta": "θ", "iota": "ι", "kappa": "κ", "lambda": "λ", "mu": "μ", "nu": "ν",
    "xi": "ξ", "pi": "π", "rho": "ρ", "varrho": "ρ", "sigma": "σ", "tau": "τ", "upsilon": "υ", "phi": "φ",
    "varphi": "φ", "chi": "χ", "psi": "ψ", "omega": "ω", "Gamma": "Γ", "Delta": "Δ", "Theta": "Θ",
    "Lambda": "Λ", "Xi": "Ξ", "Pi": "Π", "Sigma": "Σ", "Phi": "Φ", "Psi": "Ψ", "Omega": "Ω",
}
SYMBOLS = {"times": "×", "cdot": "·", "sim": "∼", "approx": "≈", "simeq": "≃", "leq": "≤", "le": "≤",
           "geq": "≥", "ge": "≥", "lt": "<", "gt": ">", "pm": "±", "propto": "∝", "infty": "∞",
           "partial": "∂", "nabla": "∇", "circ": "°", "degree": "°", "%": "%", "ll": "≪", "gg": "≫",
           "to": "→", "rightarrow": "→", "ldots": "...", "dots": "...", "cdots": "..."}
WRAPPERS = ("mathrm", "text", "textrm", "mathit", "mathbf", "boldsymbol", "operatorname", "textit", "mathsf",
            "bm", "mathnormal", "textup", "mbox", "hbox", "textnormal", "mathbb", "mathcal", "pmb")
SUPERSCRIPT = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
MINUSES = str.maketrans({"−": "-", "–": "-", "‒": "-", "—": "-", "﹣": "-", "－": "-"})


def flatten(s):
    """Plain-text form of a LaTeX or Unicode string: Greek as Unicode, ^X and _X kept, braces dropped."""
    s = s.replace("$", " ")
    s = re.sub(r"<sup>(.*?)</sup>", r"^{\1}", s)
    s = re.sub(r"<sub>(.*?)</sub>", r"_{\1}", s)
    s = re.sub(r"<br\s*/?>", " ", s)
    s = re.sub(r"</?(?:i|b|em|strong|span)[^>]*>", "", s)
    s = re.sub(r"\\tag\*?\{[^{}]*\}", " ", s)
    s = re.sub(r"\\(?:rm|it|bf|sf|left|right|big|Big|bigg|Bigg|displaystyle|textstyle)(?![A-Za-z])\s*", "", s)
    for _ in range(6):  # innermost first; nested constructs need several passes
        s0 = s
        s = re.sub(r"\\(%s)(?![A-Za-z])" % "|".join(list(GREEK) + [k for k in SYMBOLS if k.isalpha()]),
                   lambda m: GREEK.get(m.group(1)) or SYMBOLS[m.group(1)], s)
        s = re.sub(r"\\frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}", r"(\1)/(\2)", s)
        s = re.sub(r"\\(?:%s)\s*\{([^{}]*)\}" % "|".join(WRAPPERS), r"\1", s)
        s = re.sub(r"\{\}\s*\^", "^", s)
        s = re.sub(r"([_^])\s*\{\s*([^{}]*?)\s*\}", lambda m: m.group(1) + m.group(2).replace(" ", ""), s)
        s = re.sub(r"(?<![A-Za-z])\{([^{}]*)\}", r"\1", s)  # a brace after a letter belongs to a command
        if s == s0:
            break
    s = re.sub(r"\\(?:,|;|:|>|\s|quad|qquad)", " ", s)
    s = s.replace(r"\!", "").replace("~", " ")
    s = re.sub(r"\\([A-Za-z]+)", lambda m: GREEK.get(m.group(1), SYMBOLS.get(m.group(1), m.group(1))), s)
    s = s.replace("\\%", "%").replace("\\", "").replace("{", "").replace("}", "")
    s = s.replace("◦", "°").replace("º", "°").replace("˚", "°").replace("^°", "°").replace("℃", "°C")
    s = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+", lambda m: "^" + m.group(0).translate(SUPERSCRIPT), s)
    s = s.translate(MINUSES)
    return re.sub(r"[ \t]+", " ", s).strip()


# ---------------------------------------------------------------- units

_TOKEN = re.compile(r"\s*(?:(?P<deg>°\s*C)(?:\s*\^\s*(?P<dexp>[-+]?\d+)|(?P<dtexp>-?\d+))?|(?P<atom>[A-Za-zµμ]+)(?:\s*\^\s*(?P<exp>[-+]?\d+)|(?P<texp>-?\d+))?"
                    r"|(?P<one>1)(?![\d.])|(?P<op>[/()]))")


def resolve_atom(name):
    """Unit for one atom such as 'kg', 'μm', 'degC'; UnitError if it is not a unit."""
    name = ALIASES.get(name, name).replace("µ", "μ")
    if name in UNITS:
        return Unit(*UNITS[name])
    if len(name) > 1 and name[0] in PREFIXES and name[1:] in UNITS and name[1:] != "degC":
        sc, dims = UNITS[name[1:]]
        return Unit(PREFIXES[name[0]] * sc, dims)
    raise UnitError(f"not a unit: {name!r}")


def _tokens(s):
    pos, out = 0, []
    s = s.rstrip()
    while pos < len(s):
        m = _TOKEN.match(s, pos)
        if not m or m.end() == pos:
            raise UnitError(f"cannot read unit at {s[pos:pos + 12]!r} in {s!r}")
        pos = m.end()
        if m.group("deg"):
            out.append(("atom", "degC", int(m.group("dexp") or m.group("dtexp") or "1")))
        elif m.group("atom"):
            e = m.group("exp") or m.group("texp") or "1"
            out.append(("atom", m.group("atom"), int(e)))
        elif m.group("one"):
            out.append(("one", "1", 1))
        else:
            out.append(("op", m.group("op"), 1))
    return out


def _normalise_unit_text(s):
    s = flatten(s)
    s = re.sub(r"\(\s*dimensionless\s*\)", "", s).strip()
    s = s.replace("·", " ").replace("⋅", " ").replace("*", " ").replace("×", " ")
    s = re.sub(r"°\s*C", "°C", s)
    return s.strip()


def parse_unit(text):
    """Unit for a printed unit string. UnitError if any part is not understood (never guesses)."""
    if text is None:
        raise UnitError("no unit")
    s = _normalise_unit_text(text)
    if s.lower() in DIMENSIONLESS:
        return ONE
    toks = _tokens(s)
    unit, i = _parse_expr(toks, 0)
    if i != len(toks):
        raise UnitError(f"unbalanced or trailing tokens in {text!r}")
    return unit


def _parse_expr(toks, i):
    """expr := product ('/' product)*  -- a product binds every juxtaposed factor (W/m K = W/(m K))."""
    u, i = _parse_product(toks, i)
    while i < len(toks) and toks[i] == ("op", "/", 1):
        d, i = _parse_product(toks, i + 1)
        u = u * d.inverse()
    return u, i


def _parse_product(toks, i):
    u, n = ONE, 0
    while i < len(toks):
        kind, val, e = toks[i]
        if kind == "atom":
            u, i = u * resolve_atom(val) ** e, i + 1
        elif kind == "one":
            i += 1
        elif (kind, val) == ("op", "("):
            inner, i = _parse_expr(toks, i + 1)
            if i >= len(toks) or toks[i][1] != ")":
                raise UnitError("missing ')'")
            u, i = u * inner, i + 1
        else:
            break
        n += 1
    if n == 0:
        raise UnitError("empty factor")
    return u, i


# ---------------------------------------------------------------- printed values

_NUM = re.compile(r"(?<![\w.^_+-])(?:(?P<pow>10\s*\^\s*(?P<pexp>[-+]?\d+))"
                  r"|(?P<man>[-+]?\d+(?:\.\d+)?)(?:\s*[×x*·]\s*10\s*\^?\s*(?P<exp>[-+]?\d+))?)(?![\d.])")


def parse_numbers(text):
    """Every number printed in the text, with powers of ten applied, in order of appearance."""
    if not text:
        return []
    s = flatten(text)
    out = []
    for m in _NUM.finditer(s):
        if m.group("pow"):
            out.append(10.0 ** int(m.group("pexp")))
        else:
            v = float(m.group("man"))
            out.append(v * 10.0 ** int(m.group("exp")) if m.group("exp") else v)
    return out


def numbers_contained(read, marker, rel=1e-9):
    """True if every number in `read` is also among `marker` (both from parse_numbers)."""
    return all(any(math.isclose(a, b, rel_tol=rel, abs_tol=0.0) for b in marker) for a in read)


# ---------------------------------------------------------------- symbols and unit runs in flattened text

_LETTER = "A-Za-z0-9α-ωΑ-Ωϵϑϕ"


def symbol_regex(printed):
    """Regex for a printed symbol in flattened text: 'ρi' matches 'ρ_i', 'ρi', 'ρ_ i'; 'vn,j' matches 'v_n,j'.

    The first glyph is the base; the rest are subscript glyphs, which may be separated by '_' or spaces.
    A symbol with no subscript may not be followed by '_' (β must not match β_j).
    """
    glyphs = [c for c in printed if c not in "_ {}"]
    if not glyphs:
        raise ValueError("empty symbol")
    parts = [re.escape(glyphs[0])]
    for c in glyphs[1:]:
        parts.append(r"\s*,?\s*" if c == "," else r"[\s_]*" + re.escape(c))
    tail = rf"(?![{_LETTER}])" + ("" if len(glyphs) > 1 else r"(?!\s*_)")
    return re.compile(rf"(?<![{_LETTER}_\\]){''.join(parts)}{tail}")


def needs_math(printed):
    """Short ASCII symbols ('a', 'k', 'Ci') are ordinary words or units in prose; only count them in math."""
    glyphs = [c for c in printed if c not in "_ {},"]
    return all(c.isascii() for c in glyphs) and len(glyphs) <= 2


class Flat(NamedTuple):
    text: str   # one line per Markdown line; a table row stays on one line
    math: list  # [(start, end)] spans that were $...$ math, or the first cell of a table row (the symbol column)


def _flatten_line(line, offset, spans):
    out, pos = [], 0

    def emit(piece, is_math):
        start = offset + sum(len(p) for p in out) + len(out)
        out.append(piece)
        if is_math:
            spans.append((start, start + len(piece)))

    for m in re.finditer(r"\$\$(.+?)\$\$|\$(.+?)\$", line):
        emit(flatten(line[pos:m.start()]), False)
        emit(flatten(m.group(1) or m.group(2)), True)
        pos = m.end()
    emit(flatten(line[pos:]), False)
    return " ".join(out)


def flatten_markdown(md):
    """Flatten a Marker page line by line, remembering which parts were math or a table's symbol column."""
    lines, spans, offset = [], [], 0
    for line in md.splitlines():
        cells = line.strip().strip("|").split("|") if line.lstrip().startswith("|") else None
        if cells and len(cells) > 1 and not re.fullmatch(r"[\s:|-]*", line):
            first = _flatten_line(cells[0], offset, spans)
            spans.append((offset, offset + len(first)))
            rest = _flatten_line("|".join(cells[1:]), offset + len(first) + 3, spans)
            flat = first + " | " + rest
        else:
            flat = _flatten_line(line, offset, spans)
        lines.append(flat)
        offset += len(flat) + 1
    return Flat("\n".join(lines), spans)


def find_symbol(flat, printed):
    """[(start, end)] occurrences of the printed symbol in a flattened page."""
    rx = symbol_regex(printed)
    hits = [(m.start(), m.end()) for m in rx.finditer(flat.text)]
    if needs_math(printed):
        hits = [h for h in hits if any(a <= h[0] and h[1] <= b for a, b in flat.math)]
    return hits


_RUN_TOKEN = re.compile(r"°\s*C(?:\s*\^\s*[-+]?\d+|-?\d+)?|[A-Za-zµμ]+(?:\s*\^\s*[-+]?\d+|-?\d+)?|[/()]|1(?![\d.])|\S")


def unit_runs(text):
    """[(start, end, raw, Unit or None)] maximal stretches of the text made only of unit atoms, '/', '(' and ')'.

    An atom glued to a letter or written as a subscript ('x_s', 'ρa') is not a unit.
    """
    runs, cur = [], []

    def close():
        while cur and cur[-1][2] in ("/", "("):
            cur.pop()
        while cur and cur[0][2] in (")", "/"):
            cur.pop(0)
        while cur and cur[-1][2] == ")" and sum(t[2] == ")" for t in cur) > sum(t[2] == "(" for t in cur):
            cur.pop()
        while cur and cur[0][2] == "(" and sum(t[2] == "(" for t in cur) > sum(t[2] == ")" for t in cur):
            cur.pop(0)
        if cur and any(t[3] for t in cur):
            raw = text[cur[0][0]:cur[-1][1]]
            try:
                u = parse_unit(raw)
            except UnitError:
                u = None
            runs.append((cur[0][0], cur[-1][1], raw, u))
        cur.clear()

    for m in _RUN_TOKEN.finditer(text):
        tok = m.group(0)
        before = text[m.start() - 1] if m.start() > 0 else " "
        if tok in PROSE_ONLY_IN_PARSE:
            close()
            continue
        is_atom = False
        if tok.startswith("°"):
            is_atom = True
        elif re.match(r"[A-Za-zµμ]", tok) and before not in "_'’" and not re.match(rf"[{_LETTER}]", before):
            name = re.match(r"[A-Za-zµμ]+", tok).group(0)
            try:
                resolve_atom(name)
                is_atom = True
            except UnitError:
                pass
        if is_atom or tok in "/()" or (tok == "1" and all(t[2] == "(" for t in cur) and text[m.end():m.end() + 1] == "/"):
            if tok == "1" or tok in "/()":
                cur.append((m.start(), m.end(), tok, False))
            else:
                cur.append((m.start(), m.end(), tok, True))
        else:
            close()
    close()
    return runs
