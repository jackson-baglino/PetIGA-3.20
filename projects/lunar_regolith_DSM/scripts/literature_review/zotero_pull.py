#!/usr/bin/env python3
"""Copy the PDFs the user saved in Zotero into Literature/incoming/ (U2, WORKFLOW.md §3).

The user saves papers the pipeline could not get with the Zotero Connector. This script
finds them through Zotero's local API on this machine (http://localhost:23119/api/, read
only), matches them to the fetch list by DOI, and copies each stored PDF attachment into a
new folder:

    Literature/incoming/<review>/zotero_<attachment key>/<doi>.pdf

and records it in Literature/incoming/<review>/zotero_pull.json (DOI, Zotero keys, sha256,
size, path, date). A1 and A2 then check each copy like any other.

- The local API needs no key, so the Keychain is not read.
- Only localhost:23119 is contacted; the zotero-local sandbox profile enforces the same.
- Nothing in Zotero is changed.
- PDFs are untrusted: they are copied and hashed, never opened or printed. A file that does
  not start with the PDF signature is skipped and reported.
- A rerun copies nothing already recorded (same attachment, or same hash).

The fetch list is any text file with DOIs in it (fetch_list.md has one DOI link per paper);
DOIs can also be given with --doi.

    venv_workflows/bin/python scripts/literature_review/zotero_pull.py --review <name> \\
        (--fetch-list fetch_list.md | --doi 10.x/y ...) [--collection <name>] [--dry-run]
    venv_workflows/bin/python scripts/literature_review/zotero_pull.py --list-collections
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import sys
import urllib.error
import urllib.parse
import urllib.request

API = "http://localhost:23119/api/users/0"
DOCS = os.path.expanduser("~/Library/Mobile Documents/com~apple~CloudDocs/Documents")
INCOMING = os.path.join(DOCS, "Research", "Projects", "lunar_regolith_DSM", "Literature", "incoming")
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\"'<>\])]+", re.I)
STORED = {"imported_file", "imported_url"}       # attachments Zotero keeps a copy of


def norm_doi(s):
    m = DOI_RE.search(s or "")
    return m.group(0).rstrip(".,;").lower() if m else None


def dois_in(text):
    return sorted({norm_doi(m.group(0)) for m in DOI_RE.finditer(text)})


def slug(doi):
    return re.sub(r"[^a-z0-9._-]+", "_", doi)


class Local:
    """Read-only client for Zotero's local API. Refuses any other host."""

    def __init__(self, base=API):
        u = urllib.parse.urlparse(base)
        if (u.hostname, u.port) not in (("localhost", 23119), ("127.0.0.1", 23119)):
            raise SystemExit(f"zotero_pull.py talks only to localhost:23119, not {u.netloc}")
        self.base = base.rstrip("/")

    def _open(self, path, params=None):
        url = self.base + path + ("?" + urllib.parse.urlencode(params) if params else "")
        req = urllib.request.Request(url, headers={"Zotero-API-Version": "3"})
        return urllib.request.urlopen(req, timeout=30)

    def get_all(self, path, params=None):
        """Every page of a list endpoint."""
        out, start = [], 0
        while True:
            p = dict(params or {}, format="json", limit=100, start=start)
            with self._open(path, p) as r:
                batch = json.load(r)
                total = int(r.headers.get("Total-Results", len(batch)))
            out += batch
            start += len(batch)
            if not batch or start >= total:
                return out

    def file_path(self, att):
        """Local path of a stored attachment, from its enclosure link or the /file redirect."""
        href = (att.get("links", {}).get("enclosure", {}) or {}).get("href", "")
        if href.startswith("file://"):
            return urllib.request.url2pathname(urllib.parse.urlparse(href).path)

        class NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, *a, **k):
                return None
        opener = urllib.request.build_opener(NoRedirect)
        req = urllib.request.Request(f"{self.base}/items/{att['key']}/file",
                                     headers={"Zotero-API-Version": "3"})
        try:
            opener.open(req, timeout=30)
        except urllib.error.HTTPError as e:
            loc = e.headers.get("Location", "")
            if e.code in (301, 302, 303, 307) and loc.startswith("file://"):
                return urllib.request.url2pathname(urllib.parse.urlparse(loc).path)
        return None


def item_doi(data):
    d = norm_doi(data.get("DOI", ""))
    if d:
        return d
    for line in (data.get("extra") or "").splitlines():     # preprints keep the DOI in Extra
        if line.lower().startswith("doi:"):
            return norm_doi(line)
    return None


def collection_paths(z):
    """Full path of every collection ("lunar_regolith_DSM/B1") -> key."""
    cs = {c["key"]: c["data"] for c in z.get_all("/collections")}

    def path(k):
        d = cs[k]
        parent = d.get("parentCollection")
        return (path(parent) + "/" if parent in cs else "") + d["name"]
    return {path(k): k for k in cs}


def find_collection(z, name):
    paths = collection_paths(z)
    if name in paths:
        return paths[name]
    raise SystemExit(f"no Zotero collection '{name}'; give its full path, as --list-collections prints it")


def load_manifest(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"copies": []}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pull(z, wanted, review, out_root, collection=None, dry_run=False, today=None):
    """Copy stored PDFs of items whose DOI is wanted. Returns a summary dict."""
    today = today or datetime.date.today().isoformat()
    top = f"/collections/{find_collection(z, collection)}/items/top" if collection else "/items/top"
    items = [i for i in z.get_all(top) if i["data"].get("itemType") not in ("attachment", "note")]
    rdir = os.path.join(out_root, review)
    mpath = os.path.join(rdir, "zotero_pull.json")
    man = load_manifest(mpath)
    have_keys = {c["attachment_key"] for c in man["copies"]}
    have_hash = {c["sha256"] for c in man["copies"]}
    found = {}
    for it in items:
        d = item_doi(it["data"])
        if d in wanted:
            found.setdefault(d, []).append(it)
    s = {"copied": [], "already": [], "no_pdf": [], "not_in_zotero": sorted(set(wanted) - set(found)),
         "not_a_pdf": []}
    for doi in sorted(found):
        atts = []
        for it in found[doi]:
            for a in z.get_all(f"/items/{it['key']}/children"):
                ad = a["data"]
                if (ad.get("itemType") == "attachment" and ad.get("contentType") == "application/pdf"
                        and ad.get("linkMode") in STORED):
                    atts.append((it["key"], a))
        if not atts:
            s["no_pdf"].append(doi)
            continue
        for item_key, a in atts:
            if a["key"] in have_keys:
                s["already"].append(doi)
                continue
            src = z.file_path(a)
            if not src or not os.path.isfile(src):
                s["no_pdf"].append(doi)
                continue
            with open(src, "rb") as f:
                if f.read(5) != b"%PDF-":
                    s["not_a_pdf"].append(doi)
                    continue
            h = sha256(src)
            if h in have_hash:
                s["already"].append(doi)
                continue
            dest = os.path.join(rdir, f"zotero_{a['key']}", slug(doi) + ".pdf")
            if not dry_run:
                os.makedirs(os.path.dirname(dest))          # a new folder per copy
                shutil.copyfile(src, dest)
                if sha256(dest) != h:
                    raise SystemExit(f"copy of {doi} does not hash the same as Zotero's file")
                man["copies"].append({"doi": doi, "item_key": item_key, "attachment_key": a["key"],
                                      "sha256": h, "size": os.path.getsize(dest),
                                      "path": os.path.relpath(dest, rdir), "pulled": today})
                have_keys.add(a["key"])
                have_hash.add(h)
            s["copied"].append(doi)
    if s["copied"] and not dry_run:
        tmp = mpath + ".tmp"
        with open(tmp, "w") as f:
            json.dump(man, f, indent=1)
        os.replace(tmp, mpath)
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--review", help="review name: the folder under Literature/incoming/")
    ap.add_argument("--fetch-list", help="a text file with the DOIs to look for (fetch_list.md)")
    ap.add_argument("--doi", action="append", default=[], help="a DOI to look for (repeatable)")
    ap.add_argument("--collection", help="look only in this Zotero collection")
    ap.add_argument("--out-root", default=INCOMING)
    ap.add_argument("--dry-run", action="store_true", help="report what would be copied")
    ap.add_argument("--list-collections", action="store_true",
                    help="print collection names and item counts, then stop")
    args = ap.parse_args(argv)
    z = Local()
    if args.list_collections:
        counts = {c["key"]: c.get("meta", {}).get("numItems", "?") for c in z.get_all("/collections")}
        for p, k in sorted(collection_paths(z).items(), key=lambda pk: pk[0].lower()):
            print(f"{p}: {counts[k]} items")
        return 0
    if not args.review or not re.fullmatch(r"[A-Za-z0-9._-]+", args.review):
        ap.error("--review is required: letters, digits, dot, dash or underscore")
    wanted = {norm_doi(d) for d in args.doi if norm_doi(d)}
    if args.fetch_list:
        with open(args.fetch_list) as f:
            wanted |= set(dois_in(f.read()))
    if not wanted:
        ap.error("no DOIs: give --fetch-list or --doi")
    s = pull(z, wanted, args.review, args.out_root, args.collection, args.dry_run)
    verb = "would copy" if args.dry_run else "copied"
    for key, label in (("copied", verb), ("already", "already copied"), ("no_pdf", "in Zotero, no stored PDF"),
                       ("not_in_zotero", "not in Zotero"), ("not_a_pdf", "attachment is not a PDF")):
        if s[key]:
            print(f"{label} ({len(s[key])}):")
            for d in sorted(set(s[key])):
                print("  " + d)
    return 0


if __name__ == "__main__":
    sys.exit(main())
