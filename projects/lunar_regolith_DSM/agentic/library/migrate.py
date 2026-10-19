"""Build (or refresh) the library from the working files the pipeline has produced so far.  Staging area (workspace) -> library: the curator copies
validated records, texts, memos, search sets and the candidate table into the library; PDFs stay where they are and are recorded by path and hash.
  python -m agentic.library.migrate [--dry-run]"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sqlite3
from pathlib import Path

from agentic.library.store import Curator, ROOT, sha256

REVIEW = ROOT / "research" / "Review"
STATE = ROOT / ".agent_state"


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--dry-run", action="store_true"); a = ap.parse_args()
    cur = Curator(actor="curator:migrate")
    reg = sqlite3.connect(STATE / "registry.sqlite"); reg.row_factory = sqlite3.Row
    pdfs = {}
    for r in reg.execute("SELECT doi,pdf_paths FROM papers"):
        ps = [p for p in json.loads(r["pdf_paths"] or "[]") if p and os.path.exists(p)]
        if r["doi"] and ps: pdfs[r["doi"].lower()] = ps[0]
    waves, topics = {}, {}
    for wf in sorted(f for f in glob.glob(str(REVIEW / "waves" / "wave*.json")) if re.fullmatch(r"wave\d+[a-z]?", os.path.basename(f)[:-5])):
        for r in json.load(open(wf)):
            waves[r["slug"]] = os.path.basename(wf)[:-5]; topics[r["slug"]] = ",".join(r["topics"])
    memo_txt = (REVIEW / "memos" / "gate_R8_diffuse_domain.md").read_text() if (REVIEW / "memos" / "gate_R8_diffuse_domain.md").exists() else ""
    n_papers = n_variants = 0
    seen = {}
    for f in sorted(glob.glob(str(REVIEW / "records" / "*.json"))):
        if f.endswith(".verification.json"): continue
        stem = os.path.basename(f)[:-5]
        base = re.sub(r"(_(opus|sonnet|haiku|rep\d+))+$", "", stem)
        if base != stem:
            seen.setdefault(base, []).append(stem); continue
        rec = json.load(open(f)); ver = json.load(open(f[:-5] + ".verification.json"))
        md = Path(f[:-5] + ".md").read_text() if Path(f[:-5] + ".md").exists() else None
        doi = ((rec.get("paper") or {}).get("doi") or ver.get("doi") or "").lower()
        if a.dry_run: n_papers += 1; continue
        cur.publish_paper(stem, rec, ver, md, STATE / "text" / f"{stem}.txt", pdfs.get(doi), wave=waves.get(stem, "pilot"), topics=topics.get(stem, ""),
                          status="flagged" if (ver.get("checker") or {}).get("flagged") else "extracted")
        n_papers += 1
    for base, variants in seen.items():                       # replicates and model variants are kept with their paper, as extra files
        for v in variants:
            src = REVIEW / "records" / f"{v}.json"
            if not a.dry_run:
                for ext in (".json", ".verification.json"):
                    if (REVIEW / "records" / f"{v}{ext}").exists():
                        cur.publish_file(f"papers/{base}/variants/{v}{ext}", REVIEW / "records" / f"{v}{ext}", "variant")
            n_variants += 1
    n_memos = 0
    kinds = (("gate_", "gate_memo"), ("wave", "wave_report_or_triage"), ("f2_red_team", "red_team"), ("_audit", "audit"))
    for f in sorted(glob.glob(str(REVIEW / "memos" / "*"))):
        name = os.path.basename(f); stem = Path(name).stem
        kind = next((k for pre, k in kinds if pre in name), "memo")
        if not a.dry_run: cur.publish_memo(stem, Path(f), kind)
        n_memos += 1
    for extra, kind in ((REVIEW / "PAPER_LIST.md", "paper_list"), (REVIEW / "INSPECT_R8.md", "inspect_list"), (REVIEW / "REVIEW_SHEET.md", "review_sheet")):
        if extra.exists() and not a.dry_run: cur.publish_memo(extra.stem, extra, kind)
    for wf in sorted(glob.glob(str(REVIEW / "waves" / "*"))):
        if not a.dry_run: cur.publish_file(f"searches/waves/{os.path.basename(wf)}", Path(wf), "wave_manifest")
    for qf in sorted(glob.glob(str(ROOT / "agentic" / "lit" / "queries_*.json"))):
        if not a.dry_run: cur.publish_file(f"searches/queries/{os.path.basename(qf)}", Path(qf), "search_queries")
    if not a.dry_run:
        with cur.rw() as con:
            con.execute("DELETE FROM searches")
            for qf in sorted(glob.glob(str(ROOT / "agentic" / "lit" / "queries_*.json"))):
                for topic, qs in json.load(open(qf)).items():
                    for q in qs: con.execute("INSERT INTO searches VALUES (?,?,?)", (topic, q, os.path.basename(qf)))
            con.execute("DELETE FROM candidates")
            have = {r[0] for r in con.execute("SELECT lower(doi) FROM papers")}
            tri = {r["id"]: r for r in reg.execute("SELECT * FROM triage")} if reg.execute("SELECT name FROM sqlite_master WHERE name='triage'").fetchone() else {}
            meta = {r["id"]: r for r in reg.execute("SELECT * FROM meta")} if reg.execute("SELECT name FROM sqlite_master WHERE name='meta'").fetchone() else {}
            for r in reg.execute("SELECT id,doi,title,year FROM papers"):
                t = tri.get(r["id"]); m = meta.get(r["id"])
                con.execute("INSERT INTO candidates VALUES (?,?,?,?,?,?,?,?,?)", (r["id"], r["doi"], r["title"], r["year"], (m["authors"] if m else ""), ((m["abstract"] or "")[:1500] if m else ""),
                            (t["scores"] if t else ""), (t["best_topic"] if t else ""), int((r["doi"] or "").lower() in have)))
        for doc in ("POLICY.md",):
            src = Path(__file__).with_name(doc)
            if src.exists(): cur.publish_file(doc, src, "policy")
    print(f"{n_papers} papers, {n_variants} replicate/variant records, {n_memos} memos {'(dry run)' if a.dry_run else 'published'}")


if __name__ == "__main__":
    main()
