#!/usr/bin/env python3
"""X1, second part: compare the independent reads of each equation and value the memo cites (WORKFLOW.md §5).

    compare_reads.py --run <run folder> [--question Q1]

For each item blind_read.py read, the reads are:
- reader:       the claim's equation LaTeX or parameter (the library's latest version);
- checker:      the checker's blind transcription in E4 step 1 (<run>/untrusted/<paper>/blind.json), if it read it;
- blind_reader: morty-blind-reader's read (<run>/untrusted/x1/<paper>.json);
- parser:       for a numbered equation, the PDF parser's LaTeX for that number on that page (P1 equations.jsonl),
                a read made without Claude.
Not compared yet: the authors' LaTeX source (Moure & Fu only); recorded as not compared.
LaTeX is compared after tres.canon (tested equivalences only; fonts kept); a value as value plus unit, braces and
1e-3 notation ignored. An item agrees when every available read is equal and none is "not found". Any other item
is a disagreement: it opens the blocking X1 gate and goes to the review queue; the user decides it with
u3_review.py. Details (with LaTeX) go to <run>/untrusted/x1/compare.json; counts to <run>/summary_x1_<question>.json.
"""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import blind_read as br  # noqa: E402
import revise  # noqa: E402
import run_b1  # noqa: E402
import tres  # noqa: E402
from open_gate import open_gate  # noqa: E402


def vu(q):
    """A value and unit in one canonical string."""
    v = re.sub(r"(\d)[eE]([-+]?\d+)", r"\1\\times10^{\2}", str(q.get("value", "")))
    return re.sub(r"[{}]", "", tres.canon(v + str(q.get("unit_printed", ""))))


def number_key(n):
    m = re.search(r"\d+[a-z]?", str(n or ""))
    return m.group(0) if m else None


def parser_equations(run, paper):
    pdf = next(p["pdf"] for p in run_b1.papers(check_hash=False) if p["id"] == paper)
    h = hashlib.sha256(Path(pdf).read_bytes()).hexdigest()[:16]
    f = Path(run) / "untrusted" / "content" / paper / h / "equations.jsonl"
    if not f.exists():
        return {}
    out = {}
    for line in f.read_text().splitlines():
        e = json.loads(line)
        k = number_key(e.get("number"))
        if k:
            out[(e["page_pdf"], k)] = e["latex"]
    return out


def compare_item(c, checker, blind, parser):
    """{reads: {source: canonical string or 'NOT FOUND'}, agree, kind}."""
    reads = {}
    if c.get("equation"):
        kind = "equation"
        reads["reader"] = tres.canon(c["equation"]["latex"])
        for name, r in (("checker", checker), ("blind_reader", blind)):
            if r is not None:
                reads[name] = "NOT FOUND" if r.get("not_found") else tres.canon(r.get("latex", ""))
        k = number_key(c["equation"].get("number"))
        if k and (c["page_pdf"], k) in parser:
            reads["parser"] = tres.canon(parser[(c["page_pdf"], k)])
    else:
        kind = "value"
        reads["reader"] = vu(c["parameter"])
        for name, r in (("checker", checker), ("blind_reader", blind)):
            if r is not None:
                reads[name] = "NOT FOUND" if r.get("not_found") else vu(r.get("parameter") or {})
    vals = set(reads.values())
    return {"kind": kind, "reads": reads, "agree": len(vals) == 1 and "NOT FOUND" not in vals,
            "independent_reads": len(reads) - 1}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--question", default="Q1")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    version, items = br.cited_items(run, a.question)
    x1 = run / "untrusted" / "x1"
    rows = []
    for paper, c in items:
        bf = x1 / f"{paper}.json"
        if not bf.exists():
            raise SystemExit(f"no blind read for {paper}; run blind_read.py first")
        blind = {b["claim_id"]: b for b in json.loads(bf.read_text())["items"]}
        cf = run / "untrusted" / paper / "blind.json"
        checker = {b["claim_id"]: b for b in json.loads(cf.read_text())["items"]} if cf.exists() else {}
        row = compare_item(c, checker.get(c["id"]), blind.get(c["id"]), parser_equations(run, paper))
        row.update({"claim_id": c["id"], "paper": paper, "author_source": "not compared"})
        rows.append(row)
    (x1 / "compare.json").write_text(json.dumps(rows, indent=1))
    dis = [r for r in rows if not r["agree"]]
    summary = {"question": a.question, "memo_version": version, "items": len(rows),
               "agree": len(rows) - len(dis), "disagree": len(dis),
               "with_parser_read": sum("parser" in r["reads"] for r in rows),
               "with_checker_read": sum("checker" in r["reads"] for r in rows),
               "author_source_compared": 0}
    if dis:
        open_gate(run, "X1", f"{len(dis)} equation or value read(s) disagree: decide each in u3_review.py, then "
                              "approve X1", count=len(dis))
        revise.add_review(run, [{"stage": "X1", "item": r["paper"], "claim_id": r["claim_id"],
                                 "kind": "reads_disagree"} for r in dis])
    (run / f"summary_x1_{a.question}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
