#!/usr/bin/env python3
"""Build a wave manifest (papers to read in full) from the registry for human approval.
Rule: a PDF is on disk AND (triage score >= 2 on some topic OR the title matches a Wave-1 keyword topic) AND no verified record exists.
Also links PDFs found in an incoming folder to their registry rows by DOI.

  python -m agentic.lit.make_manifest --incoming research/Literature/Wave1_incoming --out research/Review/waves/wave1
"""
from __future__ import annotations

import argparse
import json
import re
import sqlite3
from pathlib import Path

from agentic.lit import common

KEYS = {
    "R1": r"thin.?interface|asymptotic|anti.?trapping|antitrapping|quantitative phase.?field|sharp.?interface|karma|echebarria|kaempfer|solidification",
    "R2": r"knudsen|free.?molecular|rarefied|vapou?r (transport|diffusion|migration)|sublimation|cold trap|diffusion",
    "R3": r"grand.?potential|free.?entropy|variational|thermodynamically consistent|multi.?phase.?field",
    "R4": r"contact angle|wetting|nucleation|surface energy|triple|three.?phase",
    "R8": r"diffuse.?domain|smoothed.?boundary|diffuse.?solid|complex geometr|embedded boundary|immersed boundary|ghost phase|pore.?space",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incoming", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    con = sqlite3.connect(common.STATE / "registry.sqlite")
    if a.incoming:
        for f in Path(a.incoming).glob("*.pdf"):
            doi = f.stem.replace("_", "/", 1).lower()
            # file names are doi with '/' and other punctuation replaced by '_'; match on the sanitised form
            for pid, d, pdfs in con.execute("SELECT id, doi, pdf_paths FROM papers WHERE doi IS NOT NULL").fetchall():
                if re.sub(r"[^A-Za-z0-9._-]+", "_", d) == f.stem or re.sub(r"[^A-Za-z0-9._-]+", "_", d).lower() == f.stem.lower():
                    lst = json.loads(pdfs)
                    if str(f) not in lst:
                        con.execute("UPDATE papers SET pdf_paths=? WHERE id=?", (json.dumps(lst + [str(f)]), pid))
        con.commit()
    have = {p.stem.split(".")[0] for p in (common.REVIEW / "records").glob("*.json") if ".verification" not in p.name}
    rows = con.execute("""SELECT p.id, p.doi, p.title, p.year, p.pdf_paths, t.scores, t.best_topic, m.authors FROM papers p
                          LEFT JOIN triage t ON t.id = p.id LEFT JOIN meta m ON m.id = p.id""").fetchall()
    out, seen = [], set()
    for pid, doi, title, year, pdfs, sc, bt, au in rows:
        pdfl = [x for x in json.loads(pdfs or "[]") if x and Path(x).exists()]
        if not pdfl or not title or pid in seen:
            continue
        score = max(json.loads(sc).values()) if sc else 0
        topics = [k for k, rx in KEYS.items() if re.search(rx, title, re.I)]
        if score < 2 and not topics:
            continue
        top = (bt if bt in KEYS or (bt and bt.startswith("R")) else None) or (topics[0] if topics else "R1")
        first = re.sub(r"[^a-z]", "", (au or "x").split(";")[0].lower()) or "x"
        slug = f"{first}_{year}_{re.sub(r'[^a-z0-9]', '', (doi or pid).lower())[-8:]}"
        if slug in have:
            continue
        seen.add(pid)
        out.append({"id": doi or pid, "slug": slug, "topics": sorted(set([top] + topics))[:3], "title": title, "year": year, "triage": score,
                    "pdf": pdfl[0], "pages": None})
    out.sort(key=lambda r: (-r["triage"], r["year"] or ""))
    outp = Path(a.out); outp.parent.mkdir(parents=True, exist_ok=True)
    outp.with_suffix(".json").write_text(json.dumps(out, indent=1))
    md = ["# Wave 1 manifest: papers to read in full (approve, trim or add before running)\n", f"{len(out)} papers with a PDF on disk. Topics decide the probe questions.\n"]
    md += [f"- [{','.join(r['topics'])}] {r['triage']} | {r['year']} | {r['title'][:100]} | `{r['slug']}`" for r in out]
    outp.with_suffix(".md").write_text("\n".join(md) + "\n")
    print(f"{len(out)} papers -> {outp}.md / .json")


if __name__ == "__main__":
    main()
