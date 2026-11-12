"""Turn plain-text symbols and equations into $LaTeX$ segments.

mathify(text) returns the text with every symbol, variable, relation and power
of ten wrapped as $...$ (LaTeX). Existing $...$ segments are left alone.
Rules run in order: literal expressions, then relations (VAR = values unit),
then powers of ten, then single tokens.
"""
import re

SUP = {"⁻": "-", "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7",
       "⁸": "8", "⁹": "9"}
SUPC = "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"


def sup(s):
    return "".join(SUP[c] for c in s)


def rm(s):
    return r"\mathrm{" + s + "}"


# ---------------------------------------------------------------- literal expressions
LIT = [
    ("Π = β0·v·ε/d0", r"\Pi = \beta_0 v \varepsilon / d_0"),
    ("β0 = 7.94×10⁵ s/m", r"\beta_0 = 7.94 \times 10^{5}\ \mathrm{s/m}"),
    ("v = |σ|/(β0 + K·l), K = ρ_ice/(D_v·ρvs)",
     r"v = |\sigma|/(\beta_0 + K l),\ K = \rho_{\mathrm{ice}}/(D_v \rho_{vs})"),
    ("v ≈ D_v·dρvs/dT·G/ρ_ice", r"v \approx D_v\,(d\rho_{vs}/dT)\,G/\rho_{\mathrm{ice}}"),
    ("D_v(T) = 2.178 × 10⁻⁵ (T/273.15)^1.81", r"D_v(T) = 2.178 \times 10^{-5}\,(T/273.15)^{1.81}"),
    ("φ²(1 − φ)²", r"\phi^{2}(1 - \phi)^{2}"),
    ("∫ −D_v(1 − φ)∇ρv·n", r"\int -D_v (1 - \phi)\, \nabla \rho_v \cdot \mathbf{n}"),
    ("∫φ dA", r"\int \phi\, dA"),
    ("−(cos θ ± sin ω)/(r·sin ω)", r"-(\cos\theta \pm \sin\omega)/(r \sin\omega)"),
    ("90° − ω < θ < 90° + ω", r"90^{\circ} - \omega < \theta < 90^{\circ} + \omega"),
    ("2·d0·cos θ/H", r"2 d_0 \cos\theta / H"),
    ("ln(1/h_open)/(s·G)", r"\ln(1/h_{\mathrm{open}})/(s G)"),
    ("σ = d0χ + β0·v_n", r"\sigma = d_0\chi + \beta_0 v_n"),
    ("ρv = ρvs(T0)(1 + d0χ)", r"\rho_v = \rho_{vs}(T_0)(1 + d_0\chi)"),
    ("ρvs(T0)(1 + d0χ)", r"\rho_{vs}(T_0)(1 + d_0\chi)"),
    ("ρvs(1 + d0χ)", r"\rho_{vs}(1 + d_0\chi)"),
    ("× (1 + d0χ)", r"\times (1 + d_0\chi)"),
    ("ρv = h_open·ρvs(T0)", r"\rho_v = h_{\mathrm{open}}\,\rho_{vs}(T_0)"),
    ("ρ_open = h_open·ρvs(T0)", r"\rho_{\mathrm{open}} = h_{\mathrm{open}}\,\rho_{vs}(T_0)"),
    ("h_lo·ρvs(T0), h_hi·ρvs(T0)", r"h_{\mathrm{lo}}\,\rho_{vs}(T_0),\ h_{\mathrm{hi}}\,\rho_{vs}(T_0)"),
    ("h·ρvs(T0)", r"h\,\rho_{vs}(T_0)"),
    ("d0(1/20 − 1/50 µm⁻¹) = 3.0 × 10⁻⁵", r"d_0(1/20 - 1/50\ \mu\mathrm{m}^{-1}) = 3.0 \times 10^{-5}"),
    ("h = 1 + d0χ_avg = 1 + 3.55 × 10⁻⁵", r"h = 1 + d_0\chi_{\mathrm{avg}} = 1 + 3.55 \times 10^{-5}"),
    ("h = 1 + d0/R0", r"h = 1 + d_0/R_0"),
    ("h = 1 ± Δh/2", r"h = 1 \pm \Delta h/2"),
    ("T0 ∓ G·L/2", r"T_0 \mp G L/2"),
    ("0.5× β0/σ", r"0.5 \times \beta_0/\sigma"),
    ("dtmax = 2τ_sub", r"\mathrm{dtmax} = 2\tau_{\mathrm{sub}}"),
    ("dtmax/τ_sub", r"\mathrm{dtmax}/\tau_{\mathrm{sub}}"),
    ("h = ρv/ρvs", r"h = \rho_v/\rho_{vs}"),
    ("σ = h − 1", r"\sigma = h - 1"),
    ("ε/R = 1/25", r"\varepsilon/R = 1/25"),
    ("R_min/ε = 29", r"R_{\mathrm{min}}/\varepsilon = 29"),
    ("R_min/25", r"R_{\mathrm{min}}/25"),
    ("ρv/ρvs(T)", r"\rho_v/\rho_{vs}(T)"),
    ("G* ≈ 1.6 K/m", r"G^{*} \approx 1.6\ \mathrm{K/m}"),
    ("20 × G*", r"20 \times G^{*}"),
    ("d0χ ≈ 0.8–6 × 10⁻⁵", r"d_0\chi \approx 0.8\text{–}6 \times 10^{-5}"),
    ("ρvs(T_end)", r"\rho_{vs}(T_{\mathrm{end}})"),
    ("(8.2 × 10⁻¹⁰ m/s)", r"(8.2 \times 10^{-10}\ \mathrm{m/s})"),
    ("1/H", r"1/H"),
    ("v_n > 0", r"v_n > 0"),
    ("D_v(T)", r"D_v(T)"),
    ("1/radius", r"1/\mathrm{radius}"),
]

# ---------------------------------------------------------------- tokens
TOK = [
    (r"ρvs\(T0\)", r"\rho_{vs}(T_0)"), (r"ρvs\(T\)", r"\rho_{vs}(T)"), (r"ρvs", r"\rho_{vs}"),
    (r"ρv", r"\rho_v"), (r"ρ_ice", r"\rho_{\mathrm{ice}}"), (r"ρ_open", r"\rho_{\mathrm{open}}"),
    (r"β0", r"\beta_0"), (r"d0χ", r"d_0\chi"), (r"\bd0\b", r"d_0"), (r"τ_sub", r"\tau_{\mathrm{sub}}"),
    (r"α_c", r"\alpha_c"), (r"ξ_v", r"\xi_v"), (r"ξ_T", r"\xi_T"), (r"ξ", r"\xi"),
    (r"θ_ref", r"\theta_{\mathrm{ref}}"), (r"θ\(x\)", r"\theta(x)"), (r"θ", r"\theta"),
    (r"γ_ia", r"\gamma_{ia}"), (r"γ_as", r"\gamma_{as}"), (r"γ_is", r"\gamma_{is}"),
    (r"χ_avg", r"\chi_{\mathrm{avg}}"), (r"χ", r"\chi"), (r"ΔT", r"\Delta T"), (r"Δh", r"\Delta h"),
    (r"(?<![\w.])(\d+(?:\.\d+)?)ε", lambda m: m.group(1) + r"\varepsilon"), (r"ε", r"\varepsilon"),
    (r"Π\(t\)", r"\Pi(t)"), (r"Π", r"\Pi"), (r"φ", r"\phi"), (r"χ", r"\chi"), (r"σ", r"\sigma"),
    (r"ω", r"\omega"),
    (r"\bh_open\b", r"h_{\mathrm{open}}"), (r"\bh_res\b", r"h_{\mathrm{res}}"), (r"\bT0\b", r"T_0"),
    (r"\bR0\b", r"R_0"), (r"\bR_min\b", r"R_{\mathrm{min}}"), (r"\bD_v\b", r"D_v"), (r"\bw_t\b", r"w_t"),
    (r"\bw_c\b", r"w_c"), (r"\bv_n\b", r"v_n"), (r"\bl_T\b", r"l_T"), (r"\bl1\b", r"l_1"),
    (r"(?<=gradient )G\b", "G"), (r"(?<=Gradient )G\b", "G"), (r"(?<=vs )G\b", "G"),
    (r"(?<=temperature )T\b", "T"), (r"(?<=width )H\b", "H"), (r"(?<=speeds )v\b", "v"),
    (r"(?<=with )r(?= the distance)", "r"), (r"(?<=length )l\b", "l"),
    # multipliers and products of integers
    (r"(?<![\w.])(\d+) × (\d+)(?= elements)", lambda m: m.group(1) + r" \times " + m.group(2)),
    (r"(?<![\w.])(\d+(?:\.\d+)?)×", lambda m: m.group(1) + r"\times"),
    (r"≈ ([−-]?\d+(?:,\d{3}(?!\d))*(?:\.\d+)?)", lambda m: r"\approx " + m.group(1).replace("−", "-").replace(",", "{,}")),
    (r"≥ (\d+(?:,\d{3}(?!\d))*)", lambda m: r"\geq " + m.group(1).replace(",", "{,}")),
    (r"≤ (\d[\d.]*)", lambda m: r"\leq " + m.group(1)),
]

VARS = {
    "h_open": r"h_{\mathrm{open}}", "h_res": r"h_{\mathrm{res}}", "Δh": r"\Delta h", "α_c": r"\alpha_c",
    "ξ_v": r"\xi_v", "ξ_T": r"\xi_T", "β0": r"\beta_0", "R0": "R_0", "w_t": "w_t", "w_c": "w_c",
    "T0": "T_0", "ε": r"\varepsilon", "Π": r"\Pi", "θ": r"\theta", "ω": r"\omega", "σ": r"\sigma",
    "G": "G", "H": "H", "h": "h", "R": "R", "dz": "dz",
}
UNITS = {"µm": r"\mu\mathrm{m}", "mm": r"\mathrm{mm}", "K/m": r"\mathrm{K/m}", "K": r"\mathrm{K}",
         "s/m": r"\mathrm{s/m}", "s": r"\mathrm{s}", "°": None}
INT = r"\d+(?:,\d{3}(?!\d))*"
NUM = r"[−-]?10[" + SUPC + r"]+|[−-]?" + INT + r"(?:\.\d+)?(?:\s*×\s*10[" + SUPC + r"]+)?"
REL = r"=|≈|≤|≥|<|>"
VAR_RE = "|".join(sorted((re.escape(v) for v in VARS), key=len, reverse=True))
REL_RE = re.compile(
    rf"(?<![\w/])({VAR_RE})\s*({REL})\s*((?:{NUM})°?(?:(?:, | and | or |–)(?:{NUM})°?)*)"
    rf"(?:\s(µm|mm|K/m|s/m)\b|\s(K)\b(?!/))?")


def num_tex(s):
    s = s.replace("−", "-")
    m = re.fullmatch(r"(-?[\d,]*\.?\d*)\s*×\s*10([" + SUPC + r"]+)", s)
    if m:
        return m.group(1).replace(",", "{,}") + r" \times 10^{" + sup(m.group(2)) + "}"
    m = re.fullmatch(r"(-?)10([" + SUPC + r"]+)", s)
    if m:
        return m.group(1) + "10^{" + sup(m.group(2)) + "}"
    return s.replace(",", "{,}")


def values_tex(v):
    parts = re.split(r"(, | and | or |–)", v)
    out = []
    for p in parts:
        if p == ", ":
            out.append(r",\ ")
        elif p in (" and ", " or "):
            # Close the equation around the word: LibreOffice reads "and"/"or"
            # inside an equation as its logic keywords and draws a wedge/vee.
            out.append("$" + p + "$")
        elif p == "–":
            out.append(r"\text{–}")
        else:
            deg = p.endswith("°")
            out.append(num_tex(p.rstrip("°")) + (r"^{\circ}" if deg else ""))
    return "".join(out)


RELTEX = {"=": "=", "≈": r"\approx", "≤": r"\leq", "≥": r"\geq", "<": "<", ">": ">"}


def mathify(text):
    """Return text with symbols and equations wrapped in $...$."""
    store = []

    def keep(latex):
        store.append(latex)
        return f"\x00{len(store) - 1}\x00"

    # protect existing math
    text = re.sub(r"\$([^$]+)\$", lambda m: keep(m.group(1)), text)
    for plain, latex in LIT:
        if plain in text:
            text = text.replace(plain, keep(latex))

    def rel(m):
        var, op, vals, u1, u2 = m.groups()
        unit = u1 or u2
        tex = f"{VARS[var]} {RELTEX[op]} {values_tex(vals)}"
        if unit:
            tex += r"\ " + UNITS[unit]
        return keep(tex)

    text = REL_RE.sub(rel, text)
    # powers of ten, alone or with a mantissa
    text = re.sub(rf"([−-]?\d[\d.]*)\s*×\s*10([{SUPC}]+)",
                  lambda m: keep(num_tex(m.group(0))), text)
    text = re.sub(rf"(?<![\w.])([−-]?)10([{SUPC}]+)", lambda m: keep(num_tex(m.group(0))), text)
    for pat, rep in TOK:
        text = re.sub(pat, (lambda m, r=rep: keep(r(m) if callable(r) else r)), text)
    # a relation sign directly before an existing math piece joins it
    rels = {"≈": r"\approx ", "≥": r"\geq ", "≤": r"\leq "}
    text = re.sub(r"([≈≥≤]) \x00(\d+)\x00", lambda m: keep(rels[m.group(1)] + store[int(m.group(2))]), text)
    # merge adjacent placeholders separated only by nothing
    out = re.sub(r"\x00(\d+)\x00", lambda m: "$" + store[int(m.group(1))] + "$", text)
    return out.replace("$$", "")
