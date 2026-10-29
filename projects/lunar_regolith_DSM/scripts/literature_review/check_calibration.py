#!/usr/bin/env python3
"""Check the user's calibration files (workflows/literature-review/calibration/).

For each file it reports, per entry: YAML errors with their line, missing or
invalid fields, pages that do not exist in the paper, and duplicates. Empty
slots are ignored, so it can be run while the files are being filled in.

Printed page numbers are mapped to PDF pages from the numbers printed in each
page's header and footer: the offset (printed minus PDF index) that most pages
agree on. Physical Review style article pages (031502-4) count by their suffix.
The script reads the PDFs with mutool and never prints their text.

    venv_workflows/bin/python scripts/literature_review/check_calibration.py [--maps]
"""
import argparse
import collections
import os
import re
import subprocess
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAL = os.path.join(ROOT, "workflows", "literature-review", "calibration")
DOCS = os.path.expanduser("~/Library/Mobile Documents/com~apple~CloudDocs/Documents")

PLANT_KINDS = {"fabricated_quote", "wrong_page", "equation_sign_or_exponent", "equation_symbol_swap",
               "number_wrong_quantity", "wrong_unit", "wrong_table_cell", "hedge_made_certain",
               "wrong_basis", "wrong_kind", "wrong_question", "dropped_condition", "scope_widened",
               "wrong_attribution", "wrong_figure"}
CLAIM_KINDS = {"method", "equation", "result", "assumption", "limitation", "parameter"}
TRET_KINDS = {"exact", "paraphrase", "absent"}
BASES = {"derived", "numerically tested", "experimentally tested", "cited", "asserted"}


def blank(v):
    return v is None or (isinstance(v, (str, list, dict)) and len(v) == 0)


def filled(entry, keys):
    return any(not blank(entry.get(k)) for k in keys)


class PageMap:
    """Printed page -> PDF page, from the numbers in each page's first and last lines."""

    def __init__(self, pdf):
        self.pdf = pdf
        self.n = 0
        self.offset = None
        self.agree = 0
        if not pdf or not os.path.exists(pdf):
            return
        out = subprocess.run(["mutool", "info", pdf], capture_output=True, text=True).stdout
        m = re.search(r"Pages:\s*(\d+)", out)
        self.n = int(m.group(1)) if m else 0
        votes = collections.Counter()
        for i in range(1, self.n + 1):
            txt = subprocess.run(["mutool", "draw", "-q", "-F", "txt", "-o", "-", pdf, str(i)],
                                 capture_output=True, text=True).stdout
            lines = [l.strip() for l in txt.splitlines() if l.strip()]
            edge = lines[:3] + lines[-3:]
            cands = set()
            for l in edge:
                for a in re.findall(r"\b\d{6}-(\d{1,3})\b", l):      # 031502-4
                    cands.add(int(a))
                if len(l) <= 12:
                    for a in re.findall(r"^[A-Z]?(\d{1,4})$", l):     # a line that is just a number (R360)
                        cands.add(int(a))
                for a in re.findall(r"^(\d{1,4})\s{2,}\S|\S\s{2,}(\d{1,4})$", l):
                    cands.update(int(x) for x in a if x)
            for c in cands:
                votes[c - i] += 1
        if votes:
            self.offset, self.agree = votes.most_common(1)[0]
            if self.agree < max(2, self.n // 3):
                self.offset = None

    def to_pdf(self, printed):
        if self.offset is None:
            return None
        p = int(printed) - self.offset
        return p if 1 <= p <= self.n else 0

    def describe(self):
        if not self.n:
            return "no file"
        if self.offset is None:
            return f"{self.n} pages; printed numbers not found, so printed pages cannot be checked"
        return (f"{self.n} pages; printed {1 + self.offset} to {self.n + self.offset} "
                f"(offset {self.offset}, agreed by {self.agree} pages)")


def load_yaml(name, problems):
    path = os.path.join(CAL, name)
    try:
        with open(path) as f:
            return yaml.safe_load(f) or {}
    except yaml.YAMLError as e:
        mark = getattr(e, "problem_mark", None)
        where = f"line {mark.line + 1}" if mark else "unknown line"
        problems.append(f"{name}: not valid YAML at {where}: {getattr(e, 'problem', e)}. "
                        "A colon followed by a space inside text needs single quotes around the text.")
        return None


def papers(problems):
    data = load_yaml("papers.yaml", problems)
    out = {}
    for p in data or []:
        f = p.get("file")
        path = os.path.join(DOCS, f) if f else None
        out[p["id"]] = {"file": path, "map": None}
    return out


def page_check(pid, page, pp, where, problems, printed=True):
    if pid not in pp:
        problems.append(f"{where}: paper '{pid}' is not in papers.yaml")
        return
    info = pp[pid]
    if not info["file"]:
        return  # not acquired yet
    if info["map"] is None:
        info["map"] = PageMap(info["file"])
    pages = page if isinstance(page, list) else [page]
    for p in pages:
        if isinstance(p, str) and re.fullmatch(r"[A-Z]\d{1,4}", p.strip()):
            p = int(p.strip()[1:])                                # review pages such as R360
        if not isinstance(p, int):
            problems.append(f"{where}: page '{p}' is not a whole number")
            continue
        if printed:
            q = info["map"].to_pdf(p)
            if q == 0:
                problems.append(f"{where}: printed page {p} is not in the paper ({info['map'].describe()})")
        elif not 1 <= p <= info["map"].n:
            problems.append(f"{where}: PDF page {p} is outside 1 to {info['map'].n}")


def check_planted(pp, problems, counts):
    d = load_yaml("planted_errors.yaml", problems)
    if d is None:
        return
    seen = set()
    n = 0
    for i, e in enumerate(d.get("errors") or [], 1):
        keys = ("paper", "page", "kind", "correct_claim", "planted_claim", "note")
        if not filled(e, keys):
            continue
        n += 1
        w = f"planted_errors.yaml entry {i}"
        for k in keys:
            if blank(e.get(k)):
                problems.append(f"{w}: '{k}' is empty")
        if e.get("kind") and e["kind"] not in PLANT_KINDS:
            problems.append(f"{w}: kind '{e['kind']}' is not one of the listed kinds")
        if e.get("kind") == "dropped_condition" and blank(e.get("condition_page")):
            problems.append(f"{w}: a dropped_condition needs condition_page")
        if e.get("paper") and e.get("page") is not None:
            page_check(e["paper"], e["page"], pp, w, problems)
        if e.get("condition_page") not in (None, ""):
            page_check(e.get("paper"), e["condition_page"], pp, w + " (condition_page)", problems)
        a, b = str(e.get("correct_claim", "")).strip(), str(e.get("planted_claim", "")).strip()
        if a and a == b:
            problems.append(f"{w}: planted_claim is the same as correct_claim")
        if b in seen and b:
            problems.append(f"{w}: duplicate planted_claim")
        seen.add(b)
    counts["planted errors"] = (n, "10 or more")


def check_gold(pp, problems, counts):
    d = load_yaml("gold_statements.yaml", problems)
    if d is None:
        return
    per = collections.Counter()
    seen = set()
    for i, e in enumerate(d.get("statements") or [], 1):
        if not filled(e, ("page", "location", "kind", "statement", "symbol", "value", "unit")):
            continue
        w = f"gold_statements.yaml entry {i} ({e.get('paper')})"
        per[e.get("paper")] += 1
        for k in ("paper", "page", "kind", "statement"):
            if blank(e.get(k)):
                problems.append(f"{w}: '{k}' is empty")
        if e.get("kind") and e["kind"] not in CLAIM_KINDS:
            problems.append(f"{w}: kind '{e['kind']}' is not one of {sorted(CLAIM_KINDS)}")
        if not blank(e.get("basis")) and str(e["basis"]).strip() not in BASES:
            problems.append(f"{w}: basis '{e['basis']}' is not one of {sorted(BASES)}")
        vals = [not blank(e.get(k)) for k in ("symbol", "value", "unit")]
        if any(vals) and not all(vals[1:]):
            problems.append(f"{w}: a value needs both value and unit (symbol optional)")
        v = e.get("value")
        if not blank(v) and not re.fullmatch(r"-?\d+(\.\d+)?([eE][-+]?\d+)?", str(v).strip()):
            problems.append(f"{w}: value must be one number (1.8e6); put ranges, comparisons and "
                            "equations in the statement and leave symbol, value and unit empty")
        if e.get("page") is not None and e.get("paper"):
            page_check(e["paper"], e["page"], pp, w, problems)
        s = str(e.get("statement", "")).strip().lower()
        if s and s in seen:
            problems.append(f"{w}: duplicate statement")
        seen.add(s)
    for pid, k in sorted(per.items()):
        counts[f"gold statements, {pid}"] = (k, "see README")


def check_tres(pp, problems, counts):
    d = load_yaml("tres_pages.yaml", problems)
    if d is None:
        return
    for i, e in enumerate(d.get("pages") or [], 1):
        w = f"tres_pages.yaml entry {i} ({e.get('paper')})"
        pages = e.get("page_pdf")
        if blank(pages):
            continue
        page_check(e.get("paper"), pages, pp, w, problems, printed=False)
        if e.get("source") == "user":
            items = e.get("items") or []
            if not any(filled(it, ("label", "latex")) for it in items):
                problems.append(f"{w}: pages given but no items transcribed")
            for j, it in enumerate(items, 1):
                if filled(it, ("label", "latex")) and (blank(it.get("label")) or blank(it.get("latex"))):
                    problems.append(f"{w} item {j}: needs both label and latex")
                lt = str(it.get("latex") or "")
                if lt.count("{") != lt.count("}"):
                    problems.append(f"{w} item {j}: unbalanced braces in latex")


def check_tret(pp, problems, counts):
    d = load_yaml("tret_questions.yaml", problems)
    if d is None:
        return
    kinds = collections.Counter()
    seen = set()
    for i, e in enumerate(d.get("questions") or [], 1):
        if not filled(e, ("question", "paper", "page", "answer")):
            continue
        w = f"tret_questions.yaml entry {i}"
        k = e.get("kind")
        kinds[k] += 1
        if k not in TRET_KINDS:
            problems.append(f"{w}: kind '{k}' is not one of {sorted(TRET_KINDS)}")
        if blank(e.get("question")):
            problems.append(f"{w}: 'question' is empty")
        if k == "absent":
            if filled(e, ("paper", "page", "answer")):
                problems.append(f"{w}: an absent question leaves paper, page and answer empty")
        else:
            for f in ("paper", "page", "answer"):
                if blank(e.get(f)):
                    problems.append(f"{w}: '{f}' is empty")
            if e.get("paper") and e.get("page") is not None:
                page_check(e["paper"], e["page"], pp, w, problems)
        q = str(e.get("question", "")).strip().lower()
        if q and q in seen:
            problems.append(f"{w}: duplicate question")
        seen.add(q)
    counts["search questions"] = (sum(kinds.values()),
                                  f"25: exact {kinds['exact']}/10, paraphrase {kinds['paraphrase']}/10, "
                                  f"absent {kinds['absent']}/5")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--maps", action="store_true", help="also print each paper's page map")
    args = ap.parse_args(argv)
    problems, counts = [], {}
    pp = papers(problems)
    check_planted(pp, problems, counts)
    check_tres(pp, problems, counts)
    check_gold(pp, problems, counts)
    check_tret(pp, problems, counts)
    print("Filled in:")
    for k, (n, target) in counts.items():
        print(f"  {k}: {n} (target {target})")
    if args.maps:
        print("Page maps:")
        for pid, info in pp.items():
            if info["file"] and info["map"] is None:
                info["map"] = PageMap(info["file"])
            print(f"  {pid}: {info['map'].describe() if info['map'] else 'not acquired yet'}")
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for p in problems:
            print(" -", p)
        return 1
    print("\nNo problems.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
