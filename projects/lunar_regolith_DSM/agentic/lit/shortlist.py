#!/usr/bin/env python3
"""First-pass triage by keyword on titles (deterministic, no LLM). Writes a markdown
shortlist per topic with flags: PDF on disk? existing note (quality)? DOI verified?
The human approves or edits this list before any paper is read by a worker.

Usage: shortlist.py --db .agent_state/registry.sqlite --out shortlist.md
"""
import argparse, json, re, sqlite3

TOPICS = {
    "R1 thin-interface / asymptotics / anti-trapping": r"thin.?interface|asymptotic|anti.?trapping|antitrapping|quantitative phase.?field|sharp.?interface|karma|almgren|echebarria|kaempfer|plapp|unequal|kinetic coefficient|dendritic|solidification",
    "R2 vapor transport: Fickian / Knudsen / free-molecular": r"knudsen|free.?molecular|rarefied|mean free path|vapou?r (transport|diffusion|migration)|water vapou?r|sublimation|hertz|adsorption|cold trap|diffusion in (porous|regolith)|dusty gas|transition regime",
    "R3 variational / grand potential / entropy": r"grand.?potential|free.?entropy|entropy|variational|thermodynamically consistent|penrose|lyapunov|gradient flow|plapp|choudhury|nestler|steinbach|multi.?phase.?field|kim.*kim|suzuki",
    "R4 ice-solid interface / contact angle / wetting": r"contact angle|wetting|wetted|pre.?melt|quasi.?liquid|ice (on|nucleation)|nucleation|surface energy|interfacial energy|heterogeneous|three.?phase|triple",
}
ap = argparse.ArgumentParser()
ap.add_argument("--db", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--pins", help="JSON list of {doi, topic, why}: papers the user requires; listed first")
a = ap.parse_args()
con = sqlite3.connect(a.db)
rows = con.execute("SELECT id, doi, title, year, pdf_paths, record_path, record_quality, status FROM papers").fetchall()
out = ["# Wave 1 shortlist (keyword triage on titles; for human approval)\n",
       "Flags: PDF = file on disk; NOTE = existing note from the earlier review (stub = under 700 bytes, treat as title-level); `?doi` = DOI did not resolve.\n"]
if a.pins:
    out.append("\n## PINNED by the user\n")
    for pin in json.load(open(a.pins)):
        r = con.execute("SELECT title, year, pdf_paths, record_path FROM papers WHERE id=?", (pin["doi"].lower(),)).fetchone()
        have = ("PDF " if r and json.loads(r[2]) else "") + ("RECORD " if r and r[3] else "")
        out.append(f"- [{pin['topic']}] {(r[0] if r and r[0] else pin['doi'])[:90]} ({r[1] if r else '?'}) | {pin['doi']} | {have.strip() or 'not in registry'}\n    why: {pin['why']}")
for topic, rx in TOPICS.items():
    pat = re.compile(rx, re.I)
    hit = [r for r in rows if pat.search((r[2] or "") + " " + (r[0] or ""))]
    out.append(f"\n## {topic} ({len(hit)})\n")
    for pid, doi, title, year, pdfs, rec, rq, st in sorted(hit, key=lambda r: (r[3] or "")):
        flags = ("PDF " if json.loads(pdfs) else "") + (f"NOTE:{rq} " if rec else "") + ("?doi" if st == "doi-unverified" else "")
        out.append(f"- [{year or '----'}] {(title or pid)[:110]} | {doi or '-'} | {flags.strip() or 'candidate only (no PDF)'}")
open(a.out, "w").write("\n".join(out) + "\n")
print("\n".join(l for l in out if l.startswith("\n##")))
