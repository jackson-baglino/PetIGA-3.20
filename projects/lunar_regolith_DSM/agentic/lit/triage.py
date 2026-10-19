#!/usr/bin/env python3
"""Stage 3: score candidate papers from title and abstract only (Haiku, tiny contexts). Writes scores to the registry
(table `triage`) and a ranked markdown list for human approval. No paper is read in full here.

  python -m agentic.lit.triage --source openalex [--limit 200] [--workers 6]
"""
from __future__ import annotations

import argparse
import json
import sqlite3
from concurrent.futures import ThreadPoolExecutor

from agentic.lit import common

SCHEMA = {"type": "object", "required": ["scores", "best_topic", "reason", "needs_full_read"], "properties": {
    "scores": {"type": "object", "properties": {f"R{i}": {"type": "integer", "minimum": 0, "maximum": 3} for i in range(1, 9)}, "required": [f"R{i}" for i in range(1, 9)]},
    "best_topic": {"type": "string"}, "reason": {"type": "string"}, "needs_full_read": {"type": "boolean"}}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="openalex")
    ap.add_argument("--limit", type=int, default=300)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out-name", default="wave1_triage.md", help="ranking file under Review/memos (use a new name per wave)")
    ap.add_argument("--only-topics", default="", help="comma list, e.g. R8: list only candidates scoring >= 2 on these topics")
    a = ap.parse_args()
    con = sqlite3.connect(common.STATE / "registry.sqlite")
    con.execute("CREATE TABLE IF NOT EXISTS triage (id TEXT PRIMARY KEY, scores TEXT, best_topic TEXT, reason TEXT, needs_full_read INTEGER)")
    rows = con.execute("""SELECT p.id, p.title, p.year, m.authors, m.venue, m.abstract FROM papers p
                          LEFT JOIN meta m ON m.id = p.id
                          WHERE p.sources LIKE ? AND p.id NOT IN (SELECT id FROM triage) AND m.abstract IS NOT NULL
                          LIMIT ?""", (f"%{a.source}%", a.limit)).fetchall()
    print(f"triaging {len(rows)} candidates")

    def one(r):
        pid, title, year, au, venue, ab = r
        out, u = common.call_claude("triage-scorer", f"TITLE: {title}\nAUTHORS: {au}\nYEAR: {year}\nVENUE: {venue}\nABSTRACT: {ab}", schema=SCHEMA, label=f"triage:{pid[:30]}")
        return pid, out

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for pid, out in ex.map(one, rows):
            con.execute("INSERT OR REPLACE INTO triage VALUES (?,?,?,?,?)", (pid, json.dumps(out["scores"]), out["best_topic"], out["reason"], int(out["needs_full_read"])))
    con.commit()
    ranked = con.execute("""SELECT p.id, p.title, p.year, t.scores, t.best_topic, t.reason, m.authors FROM triage t JOIN papers p ON p.id=t.id LEFT JOIN meta m ON m.id=p.id""").fetchall()
    ranked.sort(key=lambda r: -max(json.loads(r[3]).values()))
    out = ["# Triage ranking (abstract-level; Haiku). Score 0-3 per topic; approve before any full read\n"]
    for pid, title, year, sc, bt, why, au in ranked:
        s = json.loads(sc)
        sel = [t for t in a.only_topics.split(",") if t] or list(s)
        if max([s.get(t, 0) for t in sel]) >= 2:
            out.append(f"- [{bt}] {max(s.values())} | {year} {au} | {title[:100]} | {pid} | {why}")
    (common.REVIEW / "memos").mkdir(parents=True, exist_ok=True)
    (common.REVIEW / "memos" / a.out_name).write_text("\n".join(out) + "\n")
    print(f"{sum(1 for l in out if l.startswith('- '))} candidates scored >= 2 -> memos/{a.out_name}")


if __name__ == "__main__":
    main()
