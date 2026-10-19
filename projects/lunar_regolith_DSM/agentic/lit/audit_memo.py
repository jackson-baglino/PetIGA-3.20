#!/usr/bin/env python3
"""Audit a gate memo against the papers it cites: every memo line that carries citations [record, p. N] is judged by an independent
`claim-checker` (Opus) that sees ONLY the memo line, the record's own quotes on the cited pages and the cited page texts (never the
synthesizer's reasoning). Also checks that each cited (record, page) exists in the record. Output: $LIT_REVIEW_DIR/memos/<name>_audit.md

  python -m agentic.lit.audit_memo --memo gate_R8_diffuse_domain [--workers 4] [--model opus]
"""
from __future__ import annotations

import argparse
import json
import re
from concurrent.futures import ThreadPoolExecutor

from agentic.lit import common
from agentic.lit.verify_record import split_pages
from agentic.lit.check_claims import SCHEMA

GROUP = re.compile(r"\[([^\[\]]*?\bp\.[^\[\]]*)\]")


def parse_cites(line: str):
    cites = []
    for g in GROUP.findall(line):
        rec = None
        for tok in g.split(";"):
            m = re.search(r"([a-z]+_\d{4}_[a-z0-9]+)", tok)
            if m:
                rec = m.group(1)
            for p in re.findall(r"p\.\s*(\d+)", tok):
                if rec:
                    cites.append((rec, int(p)))
    return sorted(set(cites))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--memo", required=True); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--model", default=None)
    a = ap.parse_args()
    memo = (common.REVIEW / "memos" / f"{a.memo}.md").read_text()
    rd = common.REVIEW / "records"
    units = []
    for n, line in enumerate(memo.splitlines(), 1):
        cs = parse_cites(line)
        if cs:
            units.append((n, line.strip(), cs))
    from agentic.library.desk import Desk
    desk = Desk()

    def one(u):
        n, line, cs = u
        parts, missing = [], []
        for slug, p in cs[:5]:
            try:
                page_text = desk.ask("claim-checker", f"audit memo {a.memo}", "page", paper_id=slug, page=p)["text"]
                qs = [c["quote"] for c in desk.ask("claim-checker", f"audit memo {a.memo}", "claims", paper_id=slug, page=p)["claims"]]
            except (FileNotFoundError, KeyError):
                missing.append(f"{slug} (not in the library)"); continue
            if not page_text:
                missing.append(f"{slug} p.{p} (page not in text)"); continue
            parts.append(f"--- {slug} PAGE {p} ---\nRECORD'S OWN QUOTES ON THIS PAGE:\n" + "\n".join("* " + q for q in qs if q) + f"\nPAGE TEXT:\n{page_text[:9000]}")
        if not parts:
            return n, line, {"verdict": "cannot_tell", "quote_found_on_page": False, "problem": "no cited page could be loaded: " + "; ".join(missing)}, missing
        prompt = (f"CLAIM (a sentence from a literature-review memo; judge only the factual content, ignore INFERENCE-labelled parts and project-context labels): {re.sub(r'INFERENCE.*', '', line)[:1500]}\n\n"
                  f"QUOTE: see the record's own quotes below each page.\n\nCITED PAGES:\n" + "\n\n".join(parts))
        out, _ = common.call_claude("claim-checker", prompt, model=a.model, schema=SCHEMA, label=f"auditmemo:{n}")
        return n, line, out, missing

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        res = list(ex.map(one, units))
    cnt = {}
    for _, _, o, _ in res:
        cnt[o["verdict"]] = cnt.get(o["verdict"], 0) + 1
    lines = [f"# Audit of memo {a.memo}", "", f"{len(res)} cited lines judged. Verdicts: {cnt}. Missing citations: {sum(1 for r in res if r[3])}.", ""]
    for n, line, o, miss in res:
        if o["verdict"] != "supported" or miss:
            lines += [f"- line {n}: **{o['verdict']}**. {o['problem']}", f"  memo: {line[:400]}", f"  missing: {miss}" if miss else ""]
    (common.REVIEW / "memos" / f"{a.memo}_audit.md").write_text("\n".join(lines) + "\n")
    print(len(res), "lines;", cnt)


if __name__ == "__main__":
    main()
