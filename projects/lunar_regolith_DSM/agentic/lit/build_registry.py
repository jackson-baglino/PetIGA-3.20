#!/usr/bin/env python3
"""Build the literature registry: one row per paper, merged across sources.

No LLM is used. Sources (all read-only):
  zotero     Zotero library (zotero.sqlite is copied first; the live file is locked)
  folder     any directory of PDFs, searched recursively (DOI read from pages 1-2)
  notes      per-paper note files from an earlier review (<citekey>.md), kept as records
  doi-list   plain-text files, every DOI in them becomes a candidate row

Identity: lower-case DOI when known, else a slug of the title. Rows with the same
identity are merged. The registry lives outside git (default .agent_state/).

Usage:
  build_registry.py --db .agent_state/registry.sqlite \
      --zotero ~/Zotero/zotero.sqlite --folder DIR [--folder DIR ...] \
      --notes DIR [--doi-list FILE ...]
"""
import argparse, json, os, re, shutil, sqlite3, sys, tempfile

DOI_RE = re.compile(r"\b(10\.\d{4,9}/[^\s\"<>,;()\[\]]+)", re.I)


def slug(text):
    return re.sub(r"[^a-z0-9]+", "", (text or "").lower())[:80]


def clean_doi(d):
    return d.strip().rstrip(".").lower() if d else None


def paper_id(doi, title):
    return clean_doi(doi) or ("title:" + slug(title))


SCHEMA = """
CREATE TABLE IF NOT EXISTS papers (
  id TEXT PRIMARY KEY, doi TEXT, title TEXT, year TEXT,
  pdf_paths TEXT DEFAULT '[]', sources TEXT DEFAULT '[]',
  record_path TEXT, record_quality TEXT, status TEXT DEFAULT 'new',
  topics TEXT DEFAULT '[]', triage_score TEXT);
"""


def upsert(con, pid, **kw):
    row = con.execute("SELECT doi,title,year,pdf_paths,sources,record_path,record_quality FROM papers WHERE id=?", (pid,)).fetchone()
    if row is None:
        con.execute("INSERT INTO papers(id) VALUES (?)", (pid,))
        row = (None, None, None, "[]", "[]", None, None)
    doi, title, year, pdfs, srcs, rec, rq = row
    pdfs, srcs = json.loads(pdfs), json.loads(srcs)
    doi = doi or kw.get("doi")
    title = title or kw.get("title")
    year = year or kw.get("year")
    if kw.get("pdf") and kw["pdf"] not in pdfs:
        pdfs.append(kw["pdf"])
    if kw.get("source") and kw["source"] not in srcs:
        srcs.append(kw["source"])
    rec = kw.get("record_path") or rec
    rq = kw.get("record_quality") or rq
    con.execute("UPDATE papers SET doi=?,title=?,year=?,pdf_paths=?,sources=?,record_path=?,record_quality=? WHERE id=?",
                (doi, title, year, json.dumps(pdfs), json.dumps(srcs), rec, rq, pid))


BIB = {}


def load_bib(path):
    """citekey -> (doi, title, year) from a BibTeX file (Better BibTeX export)."""
    txt = open(path, errors="ignore").read()
    for m in re.finditer(r"@\w+\{([^,\s]+),(.*?)\n\}", txt, re.S):
        body = m.group(2)
        g = lambda f: (re.search(r"\b" + f + r"\s*=\s*[{\"](.*?)[}\"],?\s*\n", body, re.S | re.I) or [None, None])[1]
        BIB[m.group(1)] = (clean_doi(g("doi")), re.sub(r"[{}]", "", g("title") or ""), g("year"))
    return len(BIB)


def from_zotero(con, db):
    tmp = tempfile.mkdtemp()
    cp = os.path.join(tmp, "z.sqlite")
    shutil.copy(db, cp)
    z = sqlite3.connect(cp)
    q = """SELECT i.itemID, f.fieldName, v.value FROM items i
           JOIN itemData d ON d.itemID=i.itemID JOIN fields f ON f.fieldID=d.fieldID
           JOIN itemDataValues v ON v.valueID=d.valueID
           WHERE f.fieldName IN ('title','DOI','date') AND i.itemID NOT IN (SELECT itemID FROM deletedItems)"""
    meta = {}
    for iid, fn, val in z.execute(q):
        meta.setdefault(iid, {})[fn] = val
    storage = os.path.join(os.path.dirname(db), "storage")
    att = {}
    for parent, key, path in z.execute(
        "SELECT a.parentItemID, i.key, a.path FROM itemAttachments a JOIN items i ON i.itemID=a.itemID "
        "WHERE a.contentType='application/pdf' AND a.parentItemID IS NOT NULL"):
        if path and path.startswith("storage:"):
            att.setdefault(parent, []).append(os.path.join(storage, key, path[len("storage:"):]))
    n = 0
    for iid, m in meta.items():
        if not m.get("title"):
            continue
        year = (m.get("date") or "")[:4]
        pdfs = att.get(iid, [None])
        for p in pdfs:
            upsert(con, paper_id(m.get("DOI"), m["title"]), doi=clean_doi(m.get("DOI")), title=m["title"], year=year, pdf=p, source="zotero")
        n += 1
    return n


def first_pages_doi(pdf):
    try:
        import pypdf
        r = pypdf.PdfReader(pdf)
        txt = " ".join((r.pages[i].extract_text() or "") for i in range(min(2, len(r.pages))))
        m = DOI_RE.search(txt)
        return clean_doi(m.group(1)) if m else None
    except Exception:
        return None


def from_folder(con, folder, label):
    n = 0
    for root, _, files in os.walk(folder):
        for f in files:
            if not f.lower().endswith(".pdf"):
                continue
            p = os.path.join(root, f)
            doi = first_pages_doi(p)
            title = os.path.splitext(f)[0]
            if not doi and title in BIB:
                doi, title = BIB[title][0], BIB[title][1] or title
            # a DOI-keyed row merges with Zotero; otherwise key by file stem
            upsert(con, paper_id(doi, title), doi=doi, title=title, pdf=p, source=label)
            n += 1
    return n


def from_notes(con, folder, label, stub_bytes=700):
    n = 0
    for f in sorted(os.listdir(folder)):
        if not f.endswith(".md"):
            continue
        p = os.path.join(folder, f)
        size = os.path.getsize(p)
        key = os.path.splitext(f)[0]
        doi = BIB[key][0] if key in BIB else None
        with open(p, errors="ignore") as fh:
            m = DOI_RE.search(fh.read(4000))
            doi = doi or (clean_doi(m.group(1)) if m else None)
        quality = "stub" if size < stub_bytes else "note"
        prev = con.execute("SELECT record_path FROM papers WHERE id=?", (paper_id(doi, key),)).fetchone()
        if prev and prev[0] and os.path.exists(prev[0]) and os.path.getsize(prev[0]) >= size:
            upsert(con, paper_id(doi, key), source=label)  # keep the larger existing note
            n += 1
            continue
        # a citekey like karmaPhaseFieldFormulationQuantitative2001 has no DOI; key by citekey slug
        upsert(con, paper_id(doi, key), doi=doi, title=key, source=label, record_path=p, record_quality=quality)
        n += 1
    return n


def from_doi_list(con, path, label):
    n = 0
    for doi in sorted(set(clean_doi(d) for d in DOI_RE.findall(open(path, errors="ignore").read()))):
        upsert(con, paper_id(doi, ""), doi=doi, source=label)
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--zotero")
    ap.add_argument("--folder", action="append", default=[], help="DIR or LABEL=DIR")
    ap.add_argument("--notes", action="append", default=[], help="DIR or LABEL=DIR")
    ap.add_argument("--bib", action="append", default=[], help="BibTeX file mapping citekeys to DOIs")
    ap.add_argument("--doi-list", action="append", default=[])
    a = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.db)), exist_ok=True)
    con = sqlite3.connect(a.db)
    con.executescript(SCHEMA)
    for b in a.bib:
        print("bib entries:", load_bib(b))
    if a.zotero:
        print("zotero items:", from_zotero(con, a.zotero))
    for spec in a.folder:
        label, d = spec.split("=", 1) if "=" in spec else ("folder:" + os.path.basename(spec.rstrip("/")), spec)
        print(f"{label}: {from_folder(con, d, label)} pdfs")
    for spec in a.notes:
        label, d = spec.split("=", 1) if "=" in spec else ("notes:" + os.path.basename(spec.rstrip("/")), spec)
        print(f"{label}: {from_notes(con, d, label)} notes")
    for f in a.doi_list:
        print(f"doi-list {os.path.basename(f)}:", from_doi_list(con, f, "doi-list:" + os.path.basename(f)))
    con.commit()
    total, with_pdf, with_rec = con.execute(
        "SELECT count(*), sum(pdf_paths!='[]'), sum(record_path IS NOT NULL) FROM papers").fetchone()
    print(f"registry: {total} papers, {with_pdf} with a PDF, {with_rec} with an existing note/record")


if __name__ == "__main__":
    sys.exit(main())
