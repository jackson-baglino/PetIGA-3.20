#!/usr/bin/env python3
"""Fill title/year/venue for registry rows that have a DOI but no title, using Crossref.
Also verifies that each DOI resolves (rows that 404 get status 'doi-unverified').
Deterministic, no LLM. Polite-pool header uses the mailto given on the command line.

Usage: fetch_metadata.py --db .agent_state/registry.sqlite --mailto you@example.edu
"""
import argparse, json, sqlite3, time, urllib.error, urllib.parse, urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("--db", required=True)
ap.add_argument("--mailto", required=True)
a = ap.parse_args()
con = sqlite3.connect(a.db)
con.execute("CREATE TABLE IF NOT EXISTS meta (id TEXT PRIMARY KEY, venue TEXT, authors TEXT, abstract TEXT)")
rows = con.execute("SELECT id, doi FROM papers WHERE doi IS NOT NULL AND (title IS NULL OR title='')").fetchall()
ok = bad = 0
for pid, doi in rows:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": f"lit-registry/0.1 (mailto:{a.mailto})"})
    try:
        m = json.load(urllib.request.urlopen(req, timeout=20))["message"]
        title = (m.get("title") or [""])[0]
        year = str((m.get("issued", {}).get("date-parts") or [[None]])[0][0] or "")
        venue = (m.get("container-title") or [""])[0]
        authors = "; ".join(f"{x.get('family','')}" for x in (m.get("author") or [])[:4])
        con.execute("UPDATE papers SET title=?, year=? WHERE id=?", (title, year, pid))
        con.execute("INSERT OR REPLACE INTO meta VALUES (?,?,?,?)", (pid, venue, authors, None))
        ok += 1
    except urllib.error.HTTPError as e:
        con.execute("UPDATE papers SET status='doi-unverified' WHERE id=?", (pid,))
        bad += 1
        print("unresolved:", doi, e.code)
    except Exception as e:
        print("error:", doi, e)
    time.sleep(0.15)
con.commit()
print(f"filled {ok}, unresolved {bad}, of {len(rows)}")
