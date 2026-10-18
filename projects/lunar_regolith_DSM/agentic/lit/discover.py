#!/usr/bin/env python3
"""Discovery: run topic queries against OpenAlex and add hits to the registry as candidates
(with abstract, for later triage). Deterministic, no LLM.

Usage: discover.py --db REG --mailto ADDR --queries queries.json [--per-query 12]
queries.json: {"R1": ["query", ...], "R2": [...]}
"""
import argparse, json, re, sqlite3, time, urllib.parse, urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("--db", required=True)
ap.add_argument("--mailto", required=True)
ap.add_argument("--queries", required=True)
ap.add_argument("--per-query", type=int, default=12)
a = ap.parse_args()
con = sqlite3.connect(a.db)
con.execute("CREATE TABLE IF NOT EXISTS meta (id TEXT PRIMARY KEY, venue TEXT, authors TEXT, abstract TEXT)")
con.execute("CREATE TABLE IF NOT EXISTS hits (id TEXT, topic TEXT, query TEXT, cited_by INTEGER, PRIMARY KEY (id, topic, query))")


def abstract(inv):
    if not inv:
        return None
    pos = {}
    for w, idx in inv.items():
        for i in idx:
            pos[i] = w
    return " ".join(pos[i] for i in sorted(pos))


new = 0
for topic, qs in json.load(open(a.queries)).items():
    for q in qs:
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(
            {"search": q, "per-page": a.per_query, "filter": "from_publication_date:1995-01-01,type:article|preprint",
             "mailto": a.mailto})
        try:
            res = json.load(urllib.request.urlopen(url, timeout=30))["results"]
        except Exception as e:
            print("query failed:", q, e)
            continue
        for w in res:
            doi = (w.get("doi") or "").replace("https://doi.org/", "").lower() or None
            title = w.get("display_name") or ""
            pid = doi or "title:" + re.sub(r"[^a-z0-9]+", "", title.lower())[:80]
            if not con.execute("SELECT 1 FROM papers WHERE id=?", (pid,)).fetchone():
                con.execute("INSERT INTO papers(id, doi, title, year, sources) VALUES (?,?,?,?,?)",
                            (pid, doi, title, str(w.get("publication_year") or ""), json.dumps(["openalex"])))
                new += 1
            auth = "; ".join((x["author"]["display_name"].split()[-1]) for x in (w.get("authorships") or [])[:4])
            venue = ((w.get("primary_location") or {}).get("source") or {}).get("display_name") or ""
            con.execute("INSERT OR IGNORE INTO meta VALUES (?,?,?,?)", (pid, venue, auth, abstract(w.get("abstract_inverted_index"))))
            con.execute("INSERT OR REPLACE INTO hits VALUES (?,?,?,?)", (pid, topic, q, w.get("cited_by_count") or 0))
        time.sleep(0.2)
con.commit()
print("new candidate papers:", new, "| registry total:", con.execute("SELECT count(*) FROM papers").fetchone()[0])
