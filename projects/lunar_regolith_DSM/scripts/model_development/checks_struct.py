#!/usr/bin/env python3
"""Structural done-checks for the model line's stages (WORKFLOW.md §3.2).

These check form and sourcing, not truth: the sections a stage must write, the ledger's typing and sourcing,
that every quoted citation is on its page (library_search.py quote), and the test list's fields. Truth is the
verifier's job. The Stop hook runs the check of the current stage (so an agent cannot end its turn with the
work unfinished), and the station runs it again after the call.

    venv_workflows/bin/python scripts/model_development/checks_struct.py <stage> --run <run folder>
        [--min-candidates N] [--tests <tests.json>] [--doc <path>]
Prints one line per problem and exits 1 if there is any; prints "OK <stage>" and exits 0 otherwise.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = ROOT / "venv_workflows" / "bin" / "python"

M1_SECTIONS = ["Problem", "Mechanisms", "Scales", "Constitutive laws", "Notation", "Out of scope",
               "Choices for the user"]
M4_SECTIONS = ["Weak form", "Discretization", "Residual", "Jacobian", "Auxiliary functions",
               "Time scales and time step", "Solver settings", "Manufactured solutions", "Layers"]
M9_SECTIONS = ["Symbols", "Problem", "Unknowns and assumptions", "Strong form", "Nondimensional form",
               "Sharp-interface limit", "Weak form", "Discretization", "Residual and Jacobian", "Time step",
               "Code map", "Tests and what each proves", "Performance choices", "Limitations"]
STEP_TYPES = {"cited", "derived", "assumption", "choice"}


def headings(text, level=2):
    return [m.group(1).strip() for m in re.finditer(rf"^{'#' * level}\s+(.+?)\s*$", text, re.M)]


def sections_in_order(path, wanted, prefix_ok=True):
    if not path.is_file():
        return [f"missing {path}"]
    hs = headings(path.read_text())
    pos, problems = -1, []
    for w in wanted:
        idx = [i for i, h in enumerate(hs) if h == w or (prefix_ok and h.lower().startswith(w.lower()))]
        if not idx:
            problems.append(f"{path.name}: no '## {w}' section")
            continue
        if idx[0] < pos:
            problems.append(f"{path.name}: '## {w}' is out of order")
        pos = max(pos, idx[0])
    return problems


def quote_ok(quote, paper, page, index=None):
    cmd = [str(PY), str(HERE / "library_search.py")]
    if index:
        cmd += ["--index", str(index)]
    cmd += ["quote", quote, "--paper", str(paper), "--page", str(page)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode == 0


def check_citations(run, index=None):
    path = Path(run) / "derivation" / "citations.json"
    if not path.is_file():
        return [f"missing {path}"]
    try:
        cits = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        return [f"citations.json is not JSON: {e}"]
    if not isinstance(cits, list):
        return ["citations.json is not a list"]
    problems, ids = [], set()
    for i, c in enumerate(cits):
        missing = [k for k in ("id", "paper_id", "page_pdf", "quote") if k not in c]
        if missing:
            problems.append(f"citation {i}: missing {missing}")
            continue
        if c["id"] in ids:
            problems.append(f"citation id {c['id']} repeated")
        ids.add(c["id"])
        if not quote_ok(c["quote"], c["paper_id"], c["page_pdf"], index):
            problems.append(f"citation {c['id']}: quote not found on {c['paper_id']} p. {c['page_pdf']}")
    return problems


def check_ledger(run, index=None):
    d = Path(run) / "derivation"
    path = d / "ledger.json"
    if not path.is_file():
        return [f"missing {path}"]
    try:
        led = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        return [f"ledger.json is not JSON: {e}"]
    problems, seen = [], set()
    symbols = led.get("symbols")
    if not isinstance(symbols, list) or not symbols:
        problems.append("ledger: no symbols list")
    steps = led.get("steps")
    if not isinstance(steps, list) or not steps:
        return problems + ["ledger: no steps"]
    for s in steps:
        sid, typ = s.get("id"), s.get("type")
        if not sid:
            problems.append(f"step without id: {str(s)[:80]}")
            continue
        if sid in seen:
            problems.append(f"step {sid}: id repeated")
        if typ not in STEP_TYPES:
            problems.append(f"step {sid}: type {typ!r}")
        if not s.get("statement"):
            problems.append(f"step {sid}: no statement")
        if typ == "cited":
            if not all(s.get(k) not in (None, "") for k in ("paper_id", "page_pdf", "quote")):
                problems.append(f"step {sid}: cited without paper_id, page_pdf and quote")
            elif not quote_ok(s["quote"], s["paper_id"], s["page_pdf"], index):
                problems.append(f"step {sid}: quote not found on {s['paper_id']} p. {s['page_pdf']}")
        elif typ == "derived":
            frm = s.get("from") or []
            if not frm or not s.get("operation"):
                problems.append(f"step {sid}: derived without 'from' and 'operation'")
            for f in frm:
                if f not in seen:
                    problems.append(f"step {sid}: cites {f}, which is not an earlier step")
        elif typ == "assumption" and not (s.get("label") and s.get("if_wrong")):
            problems.append(f"step {sid}: assumption without label and if_wrong")
        elif typ == "choice" and not (s.get("label") and s.get("options")):
            problems.append(f"step {sid}: choice without label and options")
        seen.add(sid)
    md = d / "derivation.md"
    if not md.is_file():
        problems.append("missing derivation.md")
    else:
        text = md.read_text()
        absent = [sid for sid in seen if not re.search(rf"(?<![\w-]){re.escape(sid)}(?![\w-])", text)]
        if absent:
            problems.append(f"derivation.md does not mention steps {sorted(absent)[:20]}")
    return problems


def check_tests(tests_path):
    path = Path(tests_path)
    if not path.is_file():
        return [f"missing {path}"]
    try:
        tests = json.loads(path.read_text())["tests"]
    except (json.JSONDecodeError, KeyError) as e:
        return [f"tests.json: {e}"]
    problems, ids = [], set()
    for t in tests:
        tid = t.get("id")
        for k in ("id", "claim", "runner", "check", "tolerance", "tolerance_reason", "layer"):
            if not t.get(k):
                problems.append(f"test {tid}: no {k}")
        if tid in ids:
            problems.append(f"test {tid}: id repeated")
        ids.add(tid)
        if t.get("check") and not (path.parent / t["check"]).is_file():
            problems.append(f"test {tid}: check script {t['check']} missing")
        if t.get("runner") == "python" and not (path.parent / str(t.get("script"))).is_file():
            problems.append(f"test {tid}: script {t.get('script')} missing")
    return problems


def check_doc(doc):
    path = Path(doc)
    problems = sections_in_order(path, M9_SECTIONS)
    if problems:
        return problems
    text = path.read_text()
    # Symbols table: | symbol | meaning | units | where defined |
    body_start = text.find("\n## ", text.find("## Symbols") + 3)
    table = text[text.find("## Symbols"):body_start]
    rows = [r for r in table.splitlines() if r.startswith("|") and not re.match(r"^\|\s*-", r)][1:]
    if not rows:
        problems.append("Symbols table is empty")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["M1", "M2", "M3", "M4", "M5", "M9", "citations"])
    ap.add_argument("--run", required=True)
    ap.add_argument("--min-candidates", type=int, default=4)
    ap.add_argument("--tests")
    ap.add_argument("--doc")
    ap.add_argument("--index", help="library index (default: library_search.py's)")
    a = ap.parse_args(argv)
    d = Path(a.run) / "derivation"
    problems = []
    if a.stage == "M1":
        problems = sections_in_order(d / "spec_problem.md", M1_SECTIONS) + check_citations(a.run, a.index)
    elif a.stage == "M2":
        sel = d / "selection.md"
        if not sel.is_file():
            problems = [f"missing {sel}"]
        else:
            hs = headings(sel.read_text())
            n = sum(1 for h in hs if h.lower().startswith("candidate:"))
            if n < a.min_candidates:
                problems.append(f"selection.md: {n} candidate sections, need {a.min_candidates}")
            for w in ("Comparison", "Recommendation"):
                if w not in hs:
                    problems.append(f"selection.md: no '## {w}' section")
        problems += check_citations(a.run, a.index)
    elif a.stage == "M3":
        problems = check_ledger(a.run, a.index)
    elif a.stage == "M4":
        problems = sections_in_order(d / "spec.md", M4_SECTIONS)
        cm = d / "code_map.json"
        try:
            rows = json.loads(cm.read_text())
            if not isinstance(rows, list) or not rows:
                problems.append("code_map.json: empty")
            for r in rows if isinstance(rows, list) else []:
                miss = [k for k in ("term", "row", "where", "option", "layer") if k not in r]
                if miss:
                    problems.append(f"code_map entry {str(r)[:60]}: missing {miss}")
        except (OSError, json.JSONDecodeError) as e:
            problems.append(f"code_map.json: {e}")
    elif a.stage == "M5":
        problems = check_tests(a.tests)
    elif a.stage == "M9":
        problems = check_doc(a.doc)
    elif a.stage == "citations":
        problems = check_citations(a.run, a.index)
    for p in problems:
        print(p)
    if problems:
        return 1
    print(f"OK {a.stage}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
