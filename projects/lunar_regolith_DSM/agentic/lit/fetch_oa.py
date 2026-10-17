#!/usr/bin/env python3
"""Fetch open-access PDFs by DOI through Unpaywall into a NEW folder (downloads are untrusted: only parsed, never executed).
Paywalled papers are listed for you to download through Caltech access (Zotero connector) and drop in the same folder.

  python -m agentic.lit.fetch_oa --dois-file dois.txt --out research/Literature/Wave1_incoming --email you@caltech.edu
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dois-file", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--email", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    got, missing = [], []
    for doi in [l.strip() for l in open(a.dois_file) if l.strip() and not l.startswith("#")]:
        try:
            m = json.load(urllib.request.urlopen(f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}?email={a.email}", timeout=20))
            loc = m.get("best_oa_location") or {}
            url = loc.get("url_for_pdf")
            if not url:
                missing.append((doi, m.get("title", ""), "no open-access PDF")); continue
            dest = out / (re.sub(r"[^A-Za-z0-9._-]+", "_", doi) + ".pdf")
            req = urllib.request.Request(url, headers={"User-Agent": "lit-fetch/0.1"})
            data = urllib.request.urlopen(req, timeout=60).read()
            if not data.startswith(b"%PDF"):
                missing.append((doi, m.get("title", ""), "link did not return a PDF")); continue
            dest.write_bytes(data); got.append((doi, str(dest)))
        except Exception as e:  # noqa: BLE001
            missing.append((doi, "", str(e)[:80]))
        time.sleep(0.2)
    (out / "FETCH_REPORT.md").write_text("# Fetched\n" + "\n".join(f"- {d} -> {p}" for d, p in got) + "\n\n# Needs manual download (Caltech access)\n" + "\n".join(f"- {d} | {t} | {why}" for d, t, why in missing) + "\n")
    print(f"fetched {len(got)}, need manual {len(missing)} -> {out}/FETCH_REPORT.md")


if __name__ == "__main__":
    main()
