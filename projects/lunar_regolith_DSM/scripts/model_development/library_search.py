#!/usr/bin/env python3
"""Keyword search over the paper library, with verbatim results and quote checks.

The model line's math agent (beth-modeler) uses this to find passages in the
papers; it then opens the page itself before relying on one. Everything here is
deterministic: SQLite FTS5 (porter + unicode61 tokenizer), no embeddings, no
model calls. The index is built from the library files and can be rebuilt and
compared at any time.

What is indexed, from library/papers/<paper_id>/:
  - content/<copy>/pages/pNNN.md  split into paragraph chunks (blank-line blocks)
  - content/<copy>/equations.jsonl  one chunk per displayed equation
  - claims/<claim_id>.vN.json       one chunk per published claim (its quote)
Each chunk records paper_id, copy hash, PDF page, kind, path and the verbatim text.

Paper text is untrusted data. This script never interprets it, and a normal
session should not print it (project rule: library/ is read only in review mode
or by the modeler's headless profile). `selfcheck` prints counts only.

Commands (run with venv_workflows/bin/python):
  build     [--library library] [--index index/library_fts.sqlite]
  search    QUERY [--paper ID] [--kind page|equation|claim] [--k 10] [--log PATH]
  quote     TEXT --paper ID --page N     exit 0 if TEXT is on that page, else 1
  selfcheck                              every published claim's quote is found on
                                         its own page, and search ranks that page
                                         first; prints counts only

Output of search and quote is JSON on stdout. Every search is appended to the
query log (default index/library_queries.jsonl) with the time and the hit IDs.
"""
import argparse
import datetime
import hashlib
import json
import re
import sqlite3
import sys
from pathlib import Path

DEFAULT_LIBRARY = Path("library")
DEFAULT_INDEX = Path("index/library_fts.sqlite")
DEFAULT_LOG = Path("index/library_queries.jsonl")
SCHEMA_VERSION = 1

PAGE_RE = re.compile(r"^p(\d+)\.md$")
CLAIM_RE = re.compile(r"^(?P<cid>.+)\.v(?P<ver>\d+)\.json$")


def normalize(text):
    """Whitespace- and case-insensitive form used for verbatim quote matching.

    Only whitespace runs are collapsed and case folded; punctuation, digits and
    symbols must match exactly, so a changed sign or exponent fails.
    """
    return re.sub(r"\s+", " ", text).strip().casefold()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def find_key(obj, key):
    """First value of `key` anywhere in a nested dict or list (claims are wrapped
    differently by different stations)."""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            r = find_key(v, key)
            if r is not None:
                return r
    elif isinstance(obj, list):
        for v in obj:
            r = find_key(v, key)
            if r is not None:
                return r
    return None


def iter_chunks(library):
    """Yield chunk dicts from the library tree, in a fixed (sorted) order."""
    papers = library / "papers"
    for paper_dir in sorted(p for p in papers.iterdir() if p.is_dir()):
        pid = paper_dir.name
        content = paper_dir / "content"
        if content.is_dir():
            for copy_dir in sorted(c for c in content.iterdir() if c.is_dir()):
                copy = copy_dir.name
                pages = copy_dir / "pages"
                if pages.is_dir():
                    for page_file in sorted(pages.iterdir()):
                        m = PAGE_RE.match(page_file.name)
                        if not m:
                            continue
                        page = int(m.group(1))
                        text = page_file.read_text(encoding="utf-8", errors="replace")
                        blocks = [b.strip() for b in re.split(r"\n\s*\n", text) if b.strip()]
                        for i, block in enumerate(blocks, 1):
                            yield {
                                "chunk_id": f"{pid}/{copy[:12]}/p{page:03d}/b{i:03d}",
                                "paper_id": pid, "copy": copy, "page_pdf": page,
                                "kind": "page", "path": str(page_file), "text": block,
                            }
                eq_file = copy_dir / "equations.jsonl"
                if eq_file.is_file():
                    for i, line in enumerate(eq_file.read_text(encoding="utf-8").splitlines(), 1):
                        if not line.strip():
                            continue
                        try:
                            eq = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        latex = eq.get("latex") or eq.get("text") or ""
                        page = eq.get("page_pdf", eq.get("page"))
                        number = eq.get("number") or eq.get("printed_number") or ""
                        yield {
                            "chunk_id": f"{pid}/{copy[:12]}/eq{i:04d}",
                            "paper_id": pid, "copy": copy,
                            "page_pdf": int(page) if isinstance(page, (int, float, str)) and str(page).isdigit() else None,
                            "kind": "equation", "path": str(eq_file),
                            "text": (f"({number}) " if number else "") + latex,
                        }
        claims = paper_dir / "claims"
        if claims.is_dir():
            latest = {}
            for claim_file in sorted(claims.iterdir()):
                m = CLAIM_RE.match(claim_file.name)
                if not m:
                    continue
                cid, ver = m.group("cid"), int(m.group("ver"))
                if cid not in latest or ver > latest[cid][0]:
                    latest[cid] = (ver, claim_file)
            for cid, (ver, claim_file) in sorted(latest.items()):
                try:
                    rec = json.loads(claim_file.read_text(encoding="utf-8"))
                except json.JSONDecodeError:
                    continue
                quote = find_key(rec, "quote")
                page = find_key(rec, "page_pdf")
                if not quote:
                    continue
                yield {
                    "chunk_id": f"{pid}/claim/{cid}.v{ver}",
                    "paper_id": pid, "copy": "",
                    "page_pdf": int(page) if isinstance(page, int) else None,
                    "kind": "claim", "path": str(claim_file), "text": quote,
                }


def cmd_build(args):
    library, index = Path(args.library), Path(args.index)
    if not (library / "papers").is_dir():
        sys.exit(f"no papers folder under {library}")
    index.parent.mkdir(parents=True, exist_ok=True)
    tmp = index.with_suffix(".tmp")
    if tmp.exists():
        tmp.unlink()
    con = sqlite3.connect(tmp)
    con.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
    con.execute("CREATE TABLE chunks (chunk_id TEXT PRIMARY KEY, paper_id TEXT, copy TEXT,"
                " page_pdf INTEGER, kind TEXT, path TEXT, text TEXT, norm TEXT)")
    con.execute("CREATE VIRTUAL TABLE fts USING fts5(text, content='chunks',"
                " content_rowid='rowid', tokenize='porter unicode61')")
    files = {}
    n = 0
    for ch in iter_chunks(library):
        con.execute("INSERT INTO chunks VALUES (?,?,?,?,?,?,?,?)",
                    (ch["chunk_id"], ch["paper_id"], ch["copy"], ch["page_pdf"], ch["kind"],
                     ch["path"], ch["text"], normalize(ch["text"])))
        files.setdefault(ch["path"], None)
        n += 1
    con.execute("INSERT INTO fts(rowid, text) SELECT rowid, text FROM chunks")
    manifest = {p: sha256_file(p) for p in sorted(files)}
    meta = {
        "schema_version": SCHEMA_VERSION,
        "built_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "library": str(library),
        "n_chunks": n,
        "n_files": len(manifest),
        "files_sha256": hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest(),
        "sqlite_version": sqlite3.sqlite_version,
        "tokenizer": "porter unicode61",
    }
    con.executemany("INSERT INTO meta VALUES (?,?)", [(k, json.dumps(v)) for k, v in meta.items()])
    con.execute("CREATE TABLE files (path TEXT PRIMARY KEY, sha256 TEXT)")
    con.executemany("INSERT INTO files VALUES (?,?)", manifest.items())
    con.commit()
    con.close()
    tmp.replace(index)
    counts = summary_counts(index)
    print(json.dumps({"index": str(index), **meta, **counts}, indent=1))


def summary_counts(index):
    con = sqlite3.connect(index)
    rows = con.execute("SELECT kind, COUNT(*), COUNT(DISTINCT paper_id) FROM chunks GROUP BY kind").fetchall()
    papers = con.execute("SELECT COUNT(DISTINCT paper_id) FROM chunks").fetchone()[0]
    con.close()
    return {"papers": papers, "by_kind": {k: {"chunks": c, "papers": p} for k, c, p in rows}}


def fts_query(text):
    """Turn free text into a safe FTS5 query: quoted terms joined by OR, so
    punctuation in the query cannot be read as FTS syntax."""
    terms = re.findall(r"\w+", text, flags=re.UNICODE)
    if not terms:
        return None
    return " OR ".join('"' + t.replace('"', "") + '"' for t in terms)


def search(index, query, paper=None, kind=None, k=10):
    q = fts_query(query)
    if q is None:
        return []
    con = sqlite3.connect(index)
    sql = ("SELECT c.chunk_id, c.paper_id, c.copy, c.page_pdf, c.kind, c.path, c.text, bm25(fts) AS score"
           " FROM fts JOIN chunks c ON c.rowid = fts.rowid WHERE fts MATCH ?")
    params: list = [q]
    if paper:
        sql += " AND c.paper_id = ?"
        params.append(paper)
    if kind:
        sql += " AND c.kind = ?"
        params.append(kind)
    sql += " ORDER BY score, c.chunk_id LIMIT ?"
    params.append(k)
    rows = con.execute(sql, params).fetchall()
    con.close()
    keys = ["chunk_id", "paper_id", "copy", "page_pdf", "kind", "path", "text", "score"]
    return [dict(zip(keys, r)) for r in rows]


def quote_on_page(index, text, paper, page):
    """True if the normalized quote is a substring of the normalized page text
    (all page chunks of that page joined, so a quote may span paragraphs)."""
    con = sqlite3.connect(index)
    rows = con.execute("SELECT copy, text FROM chunks WHERE paper_id = ? AND page_pdf = ? AND kind = 'page'"
                       " ORDER BY chunk_id", (paper, page)).fetchall()
    con.close()
    by_copy = {}
    for copy, t in rows:
        by_copy.setdefault(copy, []).append(t)
    target = normalize(text)
    found_in = [c for c, ts in by_copy.items() if target and target in normalize(" ".join(ts))]
    return {"found": bool(found_in), "copies_with_quote": found_in, "copies_checked": sorted(by_copy)}


def log_query(log, record):
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")


def cmd_search(args):
    hits = search(Path(args.index), args.query, args.paper, args.kind, args.k)
    log_query(Path(args.log), {
        "time": datetime.datetime.now().isoformat(timespec="seconds"),
        "query": args.query, "paper": args.paper, "kind": args.kind, "k": args.k,
        "hits": [h["chunk_id"] for h in hits],
    })
    print(json.dumps({"query": args.query, "abstain": not hits, "hits": hits}, indent=1, ensure_ascii=False))


def cmd_quote(args):
    r = quote_on_page(Path(args.index), args.text, args.paper, args.page)
    print(json.dumps({"paper": args.paper, "page_pdf": args.page, **r}, indent=1))
    sys.exit(0 if r["found"] else 1)


def cmd_selfcheck(args):
    """Counts only. For each published claim with a quote and a page:
    found_on_page  - the quote is on its own page (verbatim, whitespace/case-normalized)
    top1_same_page - searching the quote within its paper ranks a chunk of that page first
    """
    index = Path(args.index)
    con = sqlite3.connect(index)
    claims = con.execute("SELECT paper_id, page_pdf, text FROM chunks WHERE kind = 'claim'"
                         " AND page_pdf IS NOT NULL").fetchall()
    con.close()
    per_paper = {}
    for pid, page, quote in claims:
        s = per_paper.setdefault(pid, {"claims": 0, "found_on_page": 0, "top1_same_page": 0})
        s["claims"] += 1
        if quote_on_page(index, quote, pid, page)["found"]:
            s["found_on_page"] += 1
        hits = search(index, quote, paper=pid, kind="page", k=1)
        if hits and hits[0]["page_pdf"] == page:
            s["top1_same_page"] += 1
    total = {k: sum(s[k] for s in per_paper.values()) for k in ("claims", "found_on_page", "top1_same_page")}
    print(json.dumps({"index": str(index), "total": total, "per_paper": per_paper}, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", default=str(DEFAULT_INDEX))
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--library", default=str(DEFAULT_LIBRARY))
    s = sub.add_parser("search")
    s.add_argument("query")
    s.add_argument("--paper")
    s.add_argument("--kind", choices=["page", "equation", "claim"])
    s.add_argument("--k", type=int, default=10)
    s.add_argument("--log", default=str(DEFAULT_LOG))
    q = sub.add_parser("quote")
    q.add_argument("text")
    q.add_argument("--paper", required=True)
    q.add_argument("--page", type=int, required=True)
    sub.add_parser("selfcheck")
    args = ap.parse_args(argv)
    {"build": cmd_build, "search": cmd_search, "quote": cmd_quote, "selfcheck": cmd_selfcheck}[args.cmd](args)


if __name__ == "__main__":
    main()
