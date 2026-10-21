"""Collector librarian: gathers candidate papers and PDFs and hands them on; it never edits content.
  snowball : add the references and the most-cited citing papers of papers already in the library to the candidate pool (OpenAlex), tagged with their source
  tier     : choose reading tiers from the triage scores (A full pipeline, B extraction with sampled checks, C catalogued from title and abstract)
  kit      : an HTML page of links (publisher, Caltech proxy, Cornell proxy, Docuserve text) for the PDFs that need YOUR library login; downloads are filed by ingest_pdfs
  pass     : give checked workspace records to the curator (the only writer), refusing to change anything already published
Every item keeps its provenance (source, query or seed, route, time).  No credentials are used or stored.
  python -m agentic.library.collector snowball|tier|kit|pass ..."""
from __future__ import annotations

import argparse
import html
import json
import re
import sqlite3
import time
import urllib.parse
import urllib.request
from pathlib import Path

from agentic.library.store import ROOT

STATE = ROOT / ".agent_state"
MAILTO = "jbaglino@caltech.edu"
CORE = ("R4", "R8", "R9", "R10", "R11")


def get(url: str, tries: int = 4):
    for i in range(tries):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": f"lit-registry/0.1 (mailto:{MAILTO})"}), timeout=40))
        except Exception:
            time.sleep(3 * (i + 1))
    return None


def abstract(inv):
    if not inv:
        return None
    pos = {}
    for w, idx in inv.items():
        for i in idx:
            pos[i] = w
    return " ".join(pos[i] for i in sorted(pos))


def add_work(con, w: dict, source: str, seed: str) -> bool:
    doi = (w.get("doi") or "").replace("https://doi.org/", "").lower() or None
    title = w.get("display_name") or ""
    if not title:
        return False
    pid = doi or "title:" + re.sub(r"[^a-z0-9]+", "", title.lower())[:80]
    new = not con.execute("SELECT 1 FROM papers WHERE id=?", (pid,)).fetchone()
    if new:
        con.execute("INSERT INTO papers(id, doi, title, year, sources) VALUES (?,?,?,?,?)", (pid, doi, title, str(w.get("publication_year") or ""), json.dumps(["openalex", source])))
    auth = "; ".join(x["author"]["display_name"].split()[-1] for x in (w.get("authorships") or [])[:4])
    venue = ((w.get("primary_location") or {}).get("source") or {}).get("display_name") or ""
    con.execute("INSERT OR IGNORE INTO meta VALUES (?,?,?,?)", (pid, venue, auth, abstract(w.get("abstract_inverted_index"))))
    con.execute("INSERT OR REPLACE INTO hits VALUES (?,?,?,?)", (pid, "SNOW", seed, w.get("cited_by_count") or 0))
    return new


def snowball(seeds: list[str], refs: int, cites: int) -> None:
    con = sqlite3.connect(STATE / "registry.sqlite")
    new = 0
    for n, doi in enumerate(seeds, 1):
        w = get(f"https://api.openalex.org/works/doi:{urllib.parse.quote(doi)}?mailto={MAILTO}")
        if not w:
            continue
        ids = [x.rsplit("/", 1)[-1] for x in (w.get("referenced_works") or [])][:refs]
        for i in range(0, len(ids), 40):
            r = get("https://api.openalex.org/works?" + urllib.parse.urlencode({"filter": "openalex:" + "|".join(ids[i:i + 40]), "per-page": 40, "mailto": MAILTO}))
            for x in (r or {}).get("results", []):
                new += add_work(con, x, "snowball:reference", doi)
        wid = w["id"].rsplit("/", 1)[-1]
        r = get("https://api.openalex.org/works?" + urllib.parse.urlencode({"filter": f"cites:{wid}", "sort": "cited_by_count:desc", "per-page": cites, "mailto": MAILTO}))
        for x in (r or {}).get("results", []):
            new += add_work(con, x, "snowball:citing", doi)
        if n % 10 == 0:
            con.commit(); print(f"  {n}/{len(seeds)} seeds, {new} new candidates", flush=True)
        time.sleep(0.15)
    con.commit(); print(f"snowball: {new} new candidates; registry {con.execute('SELECT count(*) FROM papers').fetchone()[0]}")


def tier(read_dois: set[str], out: Path, target: int) -> None:
    con = sqlite3.connect(STATE / "registry.sqlite"); con.row_factory = sqlite3.Row
    rows = con.execute("SELECT p.id,p.doi,p.title,p.year,t.scores,m.authors,m.venue FROM triage t JOIN papers p ON p.id=t.id LEFT JOIN meta m ON m.id=p.id").fetchall()
    cites = dict(con.execute("SELECT id, max(cited_by) FROM hits GROUP BY id").fetchall())
    scored = []
    for r in rows:
        if (r["doi"] or "").lower() in read_dois or r["id"].startswith("title:"):
            continue
        s = json.loads(r["scores"]); core = max(s.get(t, 0) for t in CORE); anyv = max(s.values()); n2 = sum(v >= 2 for v in s.values())
        if anyv < 2:
            continue
        key = (core, anyv, n2, cites.get(r["id"]) or 0)
        scored.append((key, r, s))
    scored.sort(key=lambda x: x[0], reverse=True)
    A = [x for x in scored if x[0][0] >= 3]
    B = [x for x in scored if x not in A and (x[0][1] >= 3 or x[0][0] >= 2)]
    C = [x for x in scored if x not in A and x not in B]
    pick = (A + B + C)[:target]
    lines = ["# Wave 4 reading plan (tiers)", "", f"{len(scored)} candidates scored 2 or more on some topic (title and abstract only). Target {target}. A = full pipeline; B = extraction with a quote check and sampled independent checks; C = catalogued from title and abstract, read later if needed.", ""]
    for name, grp in (("A", A), ("B", B), ("C", C)):
        sel = [x for x in grp if x in pick]
        lines += [f"## Tier {name}: {len(sel)} papers", "", "| year | first author | title | topics scoring 2 or 3 | DOI |", "|---|---|---|---|---|"]
        for key, r, s in sel:
            tp = ",".join(f"{k}:{v}" for k, v in s.items() if v >= 2)
            lines.append(f"| {r['year']} | {(r['authors'] or '').split(';')[0][:18]} | {(r['title'] or '')[:100]} | {tp} | {r['doi'] or ''} |")
        lines.append("")
    out.write_text("\n".join(lines) + "\n")
    json.dump({n: [x[1]["doi"] for x in grp if x in pick and x[1]["doi"]] for n, grp in (("A", A), ("B", B), ("C", C))}, open(out.with_suffix(".json"), "w"), indent=1)
    print({n: len([x for x in g if x in pick]) for n, g in (("A", A), ("B", B), ("C", C))}, "->", out)


def kit(dois: list[str], out: Path) -> None:
    prox = json.loads((Path(__file__).with_name("proxies.json")).read_text())
    con = sqlite3.connect(STATE / "registry.sqlite"); con.row_factory = sqlite3.Row
    rows = []
    for d in dois:
        r = con.execute("SELECT p.title,p.year,m.authors,m.venue FROM papers p LEFT JOIN meta m ON m.id=p.id WHERE lower(p.doi)=?", (d.lower(),)).fetchone()
        rows.append((d, r["title"] if r else d, r["year"] if r else "", (r["authors"] if r else "") or "", (r["venue"] if r else "") or ""))
    h = ["<html><head><meta charset='utf-8'><title>Retrieval kit</title><style>body{font-family:Helvetica;font-size:13px} td{padding:3px 8px;border-bottom:1px solid #ddd}</style></head><body>",
         "<h2>PDFs that need your library login</h2><p>Click a proxy link while signed in, save the PDF anywhere (for example Downloads) and run <code>python -m agentic.lit.ingest_pdfs --from ~/Downloads --hours 24 --watch --out research/Literature/Wave4_incoming</code>; it files each PDF by DOI. If neither library has it, use the Docuserve text in the last column. Nothing here logs in for you.</p>",
         "<table><tr><th>#</th><th>paper</th><th>publisher</th><th>Caltech</th><th>Cornell</th><th>Docuserve request text</th></tr>"]
    for i, (d, title, year, au, venue) in enumerate(rows, 1):
        doi_url = "https://doi.org/" + d
        cal = f"<a href='{html.escape(prox['caltech'] + urllib.parse.quote(doi_url, safe=''))}'>Caltech</a>" if prox.get("caltech") else "<i>set prefix</i>"
        cor = f"<a href='{html.escape(prox['cornell'] + urllib.parse.quote(doi_url, safe=''))}'>Cornell</a>" if prox.get("cornell") else "<i>set prefix</i>"
        req = html.escape(f"{au}, {title}. {venue} ({year}). doi:{d}")
        h.append(f"<tr><td>{i}</td><td>{html.escape(au.split(';')[0])} {year}: {html.escape(title)}</td><td><a href='{doi_url}'>DOI</a></td><td>{cal}</td><td>{cor}</td><td><small>{req}</small></td></tr>")
    h.append("</table></body></html>")
    out.write_text("\n".join(h)); print(f"{len(rows)} papers -> {out}")


def route(dois: list[str], tiers: dict, out: Path) -> None:
    """Ask the collecting librarian for a retrieval route per paper (it only proposes; this function writes the plan outside the library)."""
    from agentic.lit import common
    con = sqlite3.connect(STATE / "registry.sqlite"); con.row_factory = sqlite3.Row
    schema = {"type": "object", "required": ["items", "notes_for_human"], "properties": {"notes_for_human": {"type": "array", "items": {"type": "string"}},
              "items": {"type": "array", "items": {"type": "object", "required": ["doi", "route", "priority", "reason"], "properties": {
                  "doi": {"type": "string"}, "route": {"type": "string"}, "priority": {"type": "string"}, "reason": {"type": "string"}, "duplicate_of": {"type": ["string", "null"]}, "preprint_hint": {"type": ["string", "null"]}}}}}}
    plan, notes = [], []
    for i in range(0, len(dois), 25):
        batch = []
        for d in dois[i:i + 25]:
            r = con.execute("SELECT p.id,p.title,p.year,p.pdf_paths,m.authors,m.venue FROM papers p LEFT JOIN meta m ON m.id=p.id WHERE lower(p.doi)=?", (d.lower(),)).fetchone()
            if r:
                have = any(p_ and Path(p_).exists() for p_ in json.loads(r["pdf_paths"] or "[]"))
                batch.append({"doi": d, "title": r["title"], "authors": r["authors"], "year": r["year"], "venue": r["venue"], "open_access_pdf_on_disk": have, "tier": tiers.get(d, "C")})
        res, _ = common.call_claude("librarian-collector", json.dumps({"wanted": batch}, indent=1), schema=schema, label="library:collector")
        plan += res["items"]; notes += res["notes_for_human"]
    out.write_text(json.dumps({"plan": plan, "notes_for_human": notes}, indent=1))
    import collections
    print(collections.Counter(x["route"] for x in plan), "->", out)


def main() -> None:
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("snowball"); s.add_argument("--refs", type=int, default=40); s.add_argument("--cites", type=int, default=30)
    t = sub.add_parser("tier"); t.add_argument("--out", required=True); t.add_argument("--target", type=int, default=350)
    r = sub.add_parser("route"); r.add_argument("--tiers-json", required=True); r.add_argument("--out", required=True); r.add_argument("--tiers", default="A,B")
    k = sub.add_parser("kit"); k.add_argument("--dois-file", required=True); k.add_argument("--out", required=True)
    n = ap.parse_args()
    if n.cmd == "snowball":
        from agentic.library.desk import Desk
        con = sqlite3.connect("file:" + str(Desk().catalog) + "?mode=ro", uri=True)
        seeds = [r[0] for r in con.execute("SELECT doi FROM papers WHERE doi!=''")]
        snowball(seeds, n.refs, n.cites)
    elif n.cmd == "tier":
        con = sqlite3.connect("file:" + str(__import__("agentic.library.store", fromlist=["LIB"]).LIB / "catalog.sqlite") + "?mode=ro", uri=True)
        read = {r[0].lower() for r in con.execute("SELECT doi FROM papers WHERE doi!=''")}
        tier(read, Path(n.out), n.target)
    elif n.cmd == "route":
        t = json.load(open(n.tiers_json)); want = n.tiers.split(","); tiers = {d: k for k in want for d in t.get(k, [])}
        route(list(tiers), tiers, Path(n.out))
    elif n.cmd == "kit":
        kit([l.strip() for l in Path(n.dois_file).read_text().split() if l.strip()], Path(n.out))


if __name__ == "__main__":
    main()
