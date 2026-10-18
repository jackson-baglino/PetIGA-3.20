#!/usr/bin/env python3
"""Extract, verify and write records for chosen papers: one fresh worker per paper, full text in its context.

  python -m agentic.lit.run_extract --paper "10.1103/6jb2-wthl=ji_karma_2026:R1,R3" --model sonnet
  python -m agentic.lit.run_extract --paper "<doi-or-registry-id>=<slug>:<topics>" ... [--workers 3] [--tag _opus] [--no-doi]

Per paper: prepare text (mutool) -> worker (`paper-extractor`, schema-enforced JSON) -> deterministic verification
(verify_record) -> record .json + .md and verification .json in $LIT_REVIEW_DIR/records/. Usage goes to the ledger.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from agentic.lit import common
from agentic.lit.verify_record import verify

SCHEMA = json.loads((common.ROOT / "agentic" / "schemas" / "paper_record.schema.json").read_text())
PROBES = json.loads((common.ROOT / "agentic" / "lit" / "probes.json").read_text())
TEXT_DIR = common.STATE / "text"


def registry_row(pid: str):
    con = sqlite3.connect(common.STATE / "registry.sqlite")
    r = con.execute("SELECT id, doi, title, year, pdf_paths FROM papers WHERE id=?", (pid.lower(),)).fetchone()
    if not r:
        r = con.execute("SELECT id, doi, title, year, pdf_paths FROM papers WHERE doi=?", (pid.lower(),)).fetchone()
    return r


def build_prompt(slug, row, topics, text):
    qs = [q for t in topics for q in PROBES.get(t, [])]
    head = (f"PAPER ID: {slug}\nREGISTRY TITLE: {row[2]}\nREGISTRY DOI: {row[1]}\nREGISTRY YEAR: {row[3]}\n\n"
            "PROBE QUESTIONS (answer each in `probes`, copying its id into `qid`):\n" + "\n".join(f"P{i + 1}: {q}" for i, q in enumerate(qs)) +
            "\n\nPAPER TEXT (page-marked; text layer is imperfect):\n\n")
    return head + text


def md_record(slug, rec, ver):
    mark = {"page_exact": "OK", "page_adjacent": "OK (adjacent page)", "elsewhere": "FOUND ON ANOTHER PAGE",
            "not_found": "NOT FOUND", "quote_too_short": "quote too short"}
    cv = {c["id"]: c for c in ver["claims"]}
    L = [f"# {slug}", "", f"- title: {rec['paper']['title']}",
         f"- year: {rec['paper']['year']}  venue: {rec['paper'].get('venue', '')}  doi: {rec['paper'].get('doi', '')}",
         f"- confidence (worker): {rec['confidence']}", f"- verification: {json.dumps(ver['metrics'])}", ""]
    if ver.get("problems"):
        L += ["**Problems:** " + "; ".join(ver["problems"]), ""]
    L += ["## Problem", rec["problem"], "", "## Probes"]
    for p in rec["probes"]:
        L.append(f"- **{p.get('qid', '')} {p['question']}** [{p['status']}] {p['answer']}"
                 + (f" (p. {p['page']}) \"{p.get('quote', '')}\"" if p.get("quote") else ""))
    L += ["", "## Claims (page, verbatim quote, string-check result)"]
    for c in rec["claims"]:
        L.append(f"- {c['id']} [{c['type']}] {c['claim']} (p. {c['page']}) \"{c['quote']}\"  => "
                 f"{mark.get(cv[c['id']]['status'], cv[c['id']]['status'])}")
    L += ["", "## Equations"]
    for e in rec["equations"]:
        st = e.get("image_status")
        tag = f" [IMAGE: {st}]" if st else (" [IMAGE CHECK PENDING]" if not e["text_layer_ok"] else "")
        L.append(f"- {e['label']} (p. {e['page']}, {e['role']}){tag}: {e.get('latex_image') or e.get('latex', '')}")
    L += ["", "## Numbers"]
    for n in rec["numbers"]:
        st = n.get("image_status")
        tag = f" [IMAGE: {st}]" if st else (" [IMAGE CHECK PENDING]" if not n.get("text_layer_ok", True) else "")
        val = n.get("value_image") or f"{n['value']} {n['unit']}"
        L.append(f"- {n['quantity']}: {val} (p. {n['page']}){tag}")
    for k, t in (("Assumptions", "assumptions"), ("Key results", "key_results"), ("Open questions", "open_questions")):
        if rec.get(t):
            L += ["", f"## {k}"] + [f"- {x}" for x in rec[t]]
    L += ["", "## Relevance to us", rec["relevance_to_us"], ""]
    ic = ver.get("image_check")
    if ic:
        L += [f"## Image check ({ic['model']}): {ic['counts']}"]
        L += [f"- CORRECTED p. {c['page']} {c['id']}: was `{c['was'][:80]}` now `{c['now'][:120]}` ({c['note']})" for c in ic["corrections"]]
        L += [f"- UNREADABLE p. {u['page']} {u['id']}: {u['note']}" for u in ic["unreadable"]]
    else:
        L += [f"## Image-check queue ({len(ver['image_queue'])}), not yet checked"] + [f"- p. {q['page']}: {q['what']}" for q in ver["image_queue"]]
    return "\n".join(L) + "\n"


def run_one(spec: str, model: str | None, tag: str, do_doi: bool):
    ident, _, topics = spec.partition(":")
    pid, _, slug = ident.partition("=")
    slug = slug or pid.replace("/", "_")
    topics = [t for t in topics.split(",") if t] or ["R1"]
    row = registry_row(pid)
    if not row:
        raise SystemExit(f"{pid}: not in the registry")
    pdfs = [p for p in json.loads(row[4]) if p and Path(p).exists()]
    if not pdfs:
        raise SystemExit(f"{pid}: no PDF on disk")
    TEXT_DIR.mkdir(parents=True, exist_ok=True)
    if not (TEXT_DIR / f"{slug}.txt").exists():
        subprocess.run([sys.executable, str(common.ROOT / "agentic/lit/prepare_text.py"), "--pdf", pdfs[0], "--id", slug,
                        "--out", str(TEXT_DIR)], check=True)
    text = (TEXT_DIR / f"{slug}.txt").read_text()
    meta = json.loads((TEXT_DIR / f"{slug}.pages.json").read_text())
    rec, usage = common.call_claude("paper-extractor", build_prompt(slug, row, topics, text), model=model, schema=SCHEMA,
                                    label=f"{slug}{tag}")
    ver = verify(rec, text, meta, do_doi=do_doi)
    out = common.REVIEW / "records"
    out.mkdir(parents=True, exist_ok=True)
    stem = f"{slug}{tag}"
    (out / f"{stem}.json").write_text(json.dumps(rec, indent=1))
    (out / f"{stem}.verification.json").write_text(json.dumps(ver, indent=1))
    (out / f"{stem}.md").write_text(md_record(stem, rec, ver))
    m = ver["metrics"]
    print(f"{stem}: {m['n_claims']} claims, quote match {m['quote_match_rate']}, exact page {m['exact_page_rate']}, "
          f"not found {m['not_found']}, coverage {m['coverage']}, image queue {m['image_queue_items']} | {usage['model']} "
          f"in {usage['input'] + usage['cache_read'] + usage['cache_write']} out {usage['output']} "
          f"${usage['cost_usd_list']:.3f} {usage['seconds']}s", flush=True)
    return stem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper", action="append", required=True, help="DOI_OR_ID=slug:R1,R3")
    ap.add_argument("--model", default=None, help="override the agent's model (sonnet|opus|haiku)")
    ap.add_argument("--tag", default="", help="suffix for the record name (e.g. _opus, _rep2)")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--no-doi", action="store_true")
    a = ap.parse_args()
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(lambda s: run_one(s, a.model, a.tag, not a.no_doi), a.paper))


if __name__ == "__main__":
    main()
