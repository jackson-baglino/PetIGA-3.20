#!/usr/bin/env python3
"""Write derivation.md from a structured answer, in one fixed layout, and lint it.

    render_derivation.py <run_dir>            re-render <run_dir>/workspace/derivation.md from answer.json
    render_derivation.py <run_dir> --stdout   print it instead

The launcher calls render() after a call whose agent's guard.json has "render": "derivation". The agent never
writes derivation.md itself: it fills the answer's structured fields (result.schema.json), and this script turns
them into the document the user reads, so every derivation has the same layout (user, 2026-11-18):

    title, run line, work order (verbatim), symbols (LaTeX, meaning, units, typical value, source),
    given equations, assumptions A1..An, steps (words, then the displayed equation, then what it uses, the
    operation and the check that confirmed it), result (LaTeX and the graded SymPy form), neglected effects,
    questions.

lint() returns format problems. They are reported in PROVENANCE.md and the manifest; they do not change the
verdict, which belongs to the grader. The checks:
  - every step has words (why) and an equation (latex), and the equation has no $ delimiters;
  - text fields have no Unicode Greek or math symbols (write $\\phi$, not the Unicode letter);
  - prose outside $...$ has no ASCII symbol names (phi, rho_v, eps, ...);
  - each step's check names a function defined in workspace/scratch/check.py, or says why there is none;
  - every Greek letter or single Latin letter left in a step equation, after the declared symbols are removed,
    is reported as possibly undefined.
"""
import argparse
import json
import re
import sys
from pathlib import Path

GREEK = ("alpha beta gamma delta epsilon varepsilon zeta eta theta vartheta iota kappa lambda mu nu xi pi rho "
         "varrho sigma tau upsilon phi varphi chi psi omega Gamma Delta Theta Lambda Xi Pi Sigma Upsilon Phi Psi "
         "Omega").split()
ASCII_NAMES = set(GREEK) | {"eps", "rho_v", "rho_vs", "rho_i", "rho_ice", "rho_a", "phi_t", "T_t"}
UNICODE_MATH = re.compile(r"[\u0370-\u03ff\u2202\u2207\u222b\u2211\u2248\u2264\u2265\u00b7\u00d7\u2212\u2032]")
LETTER_OK = {"d"}  # the differential in dV, dX


def units_latex(u):
    """'m^3 J^-1 s^-1' -> \\mathrm{m^{3}\\,J^{-1}\\,s^{-1}}; '1' -> 1 (dimensionless)."""
    u = (u or "").strip()
    if u in ("", "1"):
        return "1"
    parts = []
    for tok in u.split():
        m = re.fullmatch(r"([A-Za-z]+)\^?(-?\d+(?:/\d+)?)?", tok)
        parts.append(f"{m.group(1)}^{{{m.group(2)}}}" if m and m.group(2) else tok)
    return "$\\mathrm{" + "\\,".join(parts) + "}$"


def base_name(name):
    return name.split("(")[0].strip()


def latex_of(sym):
    """The symbol's LaTeX; without one (older work orders), its name, with a Greek word as the letter."""
    if sym.get("latex"):
        return sym["latex"]
    head, _, sub = base_name(sym["name"]).partition("_")
    head = f"\\{head}" if head in GREEK else head
    return f"{head}_{{{sub}}}" if sub else head


def symbol_rows(task, answer):
    rows = []
    for s in task.get("symbols", []):
        rows.append((latex_of(s), s["meaning"], s["units"], "",
                     "work order"))
    for s in answer.get("symbols", []):
        rows.append((latex_of(s), s["meaning"], s["units"],
                     f"{s.get('typical'):.3g}" if isinstance(s.get("typical"), (int, float)) else "", "introduced"))
    return rows


def cell(text):
    return str(text).replace("|", "\\|").replace("\n", " ")


def render(answer, task, run_id=""):
    """The derivation.md text for one answer."""
    t, out = task, []
    title = t.get("task_id", "derivation")
    out += [f"# {title}", "",
            f"Run `{run_id}` · status **{answer.get('status')}** · written by `scripts/agentkit/render_derivation.py` "
            "from `answer.json`; do not edit by hand.", ""]
    out += ["## Work order (verbatim)", "", "> " + cell(t.get("statement", "")), ""]
    out += ["## Symbols", "", "| Symbol | Meaning | Units | Typical value (SI) | Source |", "|---|---|---|---|---|"]
    for lat, meaning, units, typ, src in symbol_rows(t, answer):
        out.append(f"| ${lat}$ | {cell(meaning)} | {units_latex(units)} | {typ} | {src} |")
    out.append("")
    if t.get("given"):
        out += ["## Given", ""] + [f"- `{g}`" for g in t["given"]] + [""]
    out += ["## Assumptions", ""]
    out += [f"- **A{i}.** {cell(a)}" for i, a in enumerate(answer.get("assumptions_used", []), 1)] or ["- none"]
    out += ["", "## Derivation", ""]
    for s in answer.get("steps", []):
        out += [f"### {s.get('id')}. {s.get('title', '')}".rstrip(), "", s.get("why", s.get("statement", "")), ""]
        if s.get("latex"):
            out += ["$$", s["latex"].strip(), "$$", ""]
        chk = s.get("check") or {}
        check = (f"`{chk['function']}` in `scratch/check.py`: {chk.get('printed', '')}" if chk.get("function")
                 else f"none ({chk.get('printed', 'not stated')})")
        out += [f"*Uses:* {', '.join(s.get('from', [])) or 'nothing'}. *Operation:* {s.get('operation', '')}. "
                f"*Checked by:* {check}", ""]
    out += ["## Result", ""]
    if answer.get("result_latex"):
        out += ["$$", answer["result_latex"].strip(), "$$", ""]
    out += ["Graded SymPy form:", "", "```", answer.get("result_expr", ""), "```", ""]
    out += ["## Neglected effects", "", "| Effect | Why it is left out |", "|---|---|"]
    out += [f"| {cell(n['effect'])} | {cell(n['why'])} |" for n in answer.get("neglected", [])] or ["| none | |"]
    out.append("")
    if answer.get("questions"):
        out += ["## Questions", ""] + [f"- {cell(q)}" for q in answer["questions"]] + [""]
    return "\n".join(out)


def strip_math(text):
    return re.sub(r"\$\$.*?\$\$|\$[^$]*\$", " ", text or "", flags=re.S)


def check_functions(workspace):
    f = Path(workspace) / "scratch" / "check.py"
    return set(re.findall(r"^def\s+([A-Za-z_]\w*)\s*\(", f.read_text(), re.M)) if f.exists() else None


ATOM = re.compile(r"(\\[A-Za-z]+|[A-Za-z])(?:_(\{[^{}]*\}|[A-Za-z0-9]))?")
WRAPPERS = re.compile(r"\\(?:mathbf|boldsymbol|bm|hat|dot|ddot|tilde|bar|vec|overline)\s*(?:\{([^{}]*)\}|(\\?[A-Za-z]))")
NOTATION = {"d", "n", "V", "S", "A", "\\partial", "\\nabla", "\\delta", "\\infty"}  # differentials, normal, region


def atoms(latex):
    """The symbols of a LaTeX string as (base, subscript) pairs: \\rho_{vs} -> ('\\rho', 'vs'), h_i -> ('h', 'i').
    Text blocks are dropped, decorations (\\mathbf{q}, \\dot m) are unwrapped, and anything inside the
    subscript braces is not looked at."""
    s = re.sub(r"\\(text|mathrm|operatorname|mathit)\{[^{}]*\}", " ", latex or "")
    s = WRAPPERS.sub(lambda m: m.group(1) or m.group(2), s)
    out = []
    for m in ATOM.finditer(s):
        base, sub = m.group(1), (m.group(2) or "").strip("{}").replace(" ", "")
        out.append((base, sub))
    return out


def declared_atoms(latex):
    """The symbol a declared LaTeX string names: its atoms outside any argument list, h_i(T) -> h_i."""
    head = re.split(r"[(\[]", latex or "", maxsplit=1)[0]
    return {(b, sub) for b, sub in atoms(head)}


def leftover_symbols(latex, declared):
    """Greek letters and single Latin letters in a step equation that no declared symbol accounts for."""
    known = set()
    for d in declared:
        known |= declared_atoms(d)
    found = set()
    for base, sub in atoms(latex):
        is_symbol = (base.startswith("\\") and base[1:] in GREEK) or len(base) == 1
        derivative_or_index = sub and (sub.isdigit() or set(sub) <= set("xyztX"))   # T_t, w_x, phi_0
        if not is_symbol or base in NOTATION or (base, sub) in known or (derivative_or_index and (base, "") in known):
            continue
        found.add(base + (f"_{{{sub}}}" if sub else ""))
    return sorted(found)


def lint(answer, task, workspace):
    problems = []
    declared = [latex_of(s) for s in task.get("symbols", []) + answer.get("symbols", [])]
    declared += [m.group(1) for g in task.get("given", []) for m in [re.match(r"\s*([A-Za-z]\w*)\s*=", g)] if m]
    funcs = check_functions(workspace)
    texts = [("assumption", a) for a in answer.get("assumptions_used", [])]
    texts += [(f"neglected '{n.get('effect', '')[:30]}'", n.get("effect", "") + " " + n.get("why", ""))
              for n in answer.get("neglected", [])]
    texts += [(f"symbol {s.get('name')}", s.get("meaning", "")) for s in answer.get("symbols", [])]
    for s in answer.get("steps", []):
        sid = s.get("id", "?")
        if not (s.get("why") or "").strip():
            problems.append(f"{sid}: no words (why)")
        if not (s.get("latex") or "").strip():
            problems.append(f"{sid}: no equation (latex)")
        elif "$" in s["latex"]:
            problems.append(f"{sid}: latex contains $ delimiters")
        texts += [(f"{sid} title", s.get("title", "")), (f"{sid} why", s.get("why", ""))]
        if UNICODE_MATH.search(s.get("latex", "")):
            problems.append(f"{sid}: Unicode symbol in the equation")
        chk = s.get("check") or {}
        fn = chk.get("function")
        if fn and funcs is not None and fn not in funcs:
            problems.append(f"{sid}: check function {fn!r} is not defined in scratch/check.py")
        if fn and funcs is None:
            problems.append(f"{sid}: names check {fn!r} but there is no scratch/check.py")
        if not fn and not (chk.get("printed") or "").strip():
            problems.append(f"{sid}: no check and no reason given")
        left = leftover_symbols(s.get("latex", ""), declared)
        if left:
            problems.append(f"{sid}: possibly undefined in the equation: {', '.join(left)}")
    for where, text in texts:
        if UNICODE_MATH.search(text or ""):
            problems.append(f"{where}: Unicode Greek or math symbol in text (write LaTeX, e.g. $\\phi$)")
        names = sorted({w for w in re.findall(r"\b[A-Za-z]+(?:_[A-Za-z]+)?\b", strip_math(text))
                        if w in ASCII_NAMES})
        if names:
            problems.append(f"{where}: symbol names in prose outside $...$: {', '.join(names)}")
    if (answer.get("status") == "done") and not (answer.get("result_latex") or "").strip():
        problems.append("result: no result_latex")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--stdout", action="store_true")
    a = ap.parse_args(argv)
    run = Path(a.run_dir)
    answer = json.loads((run / "answer.json").read_text())
    task = json.loads((run / "workspace" / "task.json").read_text())
    text = render(answer, task, run.name)
    if a.stdout:
        print(text)
    else:
        (run / "workspace" / "derivation.md").write_text(text + "\n")
    for p in lint(answer, task, run / "workspace"):
        print(f"format: {p}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
