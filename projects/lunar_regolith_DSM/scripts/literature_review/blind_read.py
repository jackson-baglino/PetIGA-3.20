#!/usr/bin/env python3
"""X1, first part: a blind read of every equation and value the latest memo cites (WORKFLOW.md §3 V1, §5).

    blind_read.py --run <run folder> [--question Q1]

For each claim the latest memo cites that has an equation or a parameter (from the library's latest version of
the claim), morty-blind-reader gets only the claim ID, the equation's printed number or the quantity's name, and
the page; never the claimed LaTeX or value. One fresh call per paper, with that paper's cited pages (P1's images,
as the reader saw them). Reads go to <run>/untrusted/x1/<paper>.json; compare_reads.py compares them. Calls go to
the run's ledger (<run>/calls.jsonl); a read already made for the same memo version is reused, and
<run>/summary_x1_reads_<question>.json counts the calls made and the reads reused, by paper.
"""
import argparse
import json
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call as ac  # noqa: E402
import check_citations as cx  # noqa: E402
import read_check as rc  # noqa: E402
import synthesize as sy  # noqa: E402


def cited_items(run, question):
    """[(paper, claim)] for each equation or value claim the latest memo cites."""
    mdir = Path(run) / "untrusted" / "memos"
    k = cx.latest_version(mdir, question)
    cited = set(re.findall(r"\[([a-z0-9_]+-c\d{3})", (mdir / f"{question}.v{k}.md").read_text()))
    claims, _, _ = sy.library_claims(run)
    return k, [(c["paper"], c) for c in claims if c["id"] in cited and (c.get("equation") or c.get("parameter"))]


def item_line(c):
    if c.get("equation"):
        return f"- {c['id']}: the equation numbered {c['equation'].get('number') or '(unnumbered)'} on page " \
               f"{c['page_pdf']}"
    return f"- {c['id']}: the value of {c['parameter']['quantity']} on page {c['page_pdf']}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--question", default="Q1")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    version, items = cited_items(run, a.question)
    out = run / "untrusted" / "x1"
    out.mkdir(parents=True, exist_ok=True)
    by_paper = {}
    for paper, c in items:
        by_paper.setdefault(paper, []).append(c)
    n, reused = 0, []
    for paper, cs in sorted(by_paper.items()):
        f = out / f"{paper}.json"
        if f.exists() and json.loads(f.read_text()).get("memo_version") == version:
            reused.append(paper)                    # an earlier read of this memo version: counted, not hidden
            continue
        pages = {int(p.stem[1:]): p for p in (run / "untrusted" / paper / "pages").glob("p*.png")}
        show = [pages[p] for p in sorted({c["page_pdf"] for c in cs}) if p in pages][:rc.MAX_IMAGES]
        with tempfile.TemporaryDirectory() as tmp:
            files = rc.folder_with(show, tmp)
            res, meta = ac.call(tmp, ac.prompt("blind_reader", files=files, items="\n".join(map(item_line, cs))),
                                "blind", "morty-blind-reader", ledger=run / "calls.jsonl", step="X1", paper=paper)
        if res is None:
            raise SystemExit(f"X1 blind read failed for {paper}: {meta.get('error')}")
        f.write_text(json.dumps({"memo_version": version, **res}, indent=1))
        n += 1
    summary = {"question": a.question, "memo_version": version, "items": len(items), "papers": len(by_paper),
               "calls": n, "reused": len(reused), "reused_papers": reused}
    (run / f"summary_x1_reads_{a.question}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
