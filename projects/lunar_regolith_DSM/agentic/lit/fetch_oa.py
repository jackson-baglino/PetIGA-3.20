#!/usr/bin/env python3
"""Fetch open-access PDFs by DOI, trying every legal open route, into a NEW folder (downloads are untrusted: they are only
parsed, never executed). Routes, in order: Unpaywall best location and every other OA location, OpenAlex locations
(repositories, preprints), Semantic Scholar openAccessPdf. Anything still missing goes to RETRIEVAL_CHECKLIST.md for you
to download through Caltech or Cornell access (see agentic/lit/ingest_pdfs.py for the drop-in step).

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

PUBLISHERS = {"10.1103": "APS (Physical Review)", "10.1021": "ACS", "10.1016": "Elsevier", "10.1002": "Wiley", "10.1088": "IOP",
              "10.1017": "Cambridge", "10.1126": "Science/AAAS", "10.1038": "Nature/Springer", "10.1007": "Springer", "10.1029": "AGU/Wiley",
              "10.1063": "AIP", "10.1137": "SIAM", "10.3847": "AAS/IOP", "10.1098": "Royal Society", "10.1061": "ASCE", "10.1149": "ECS",
              "10.1115": "ASME", "10.1080": "Taylor & Francis", "10.1111": "Wiley", "10.1109": "IEEE"}


def get(url, timeout=25):
    req = urllib.request.Request(url, headers={"User-Agent": "lit-fetch/0.2"})
    return urllib.request.urlopen(req, timeout=timeout).read()


def candidates(doi: str, email: str):
    urls = []
    try:
        m = json.loads(get(f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}?email={email}"))
        for loc in [m.get("best_oa_location")] + (m.get("oa_locations") or []):
            if loc and loc.get("url_for_pdf"):
                urls.append(("unpaywall", loc["url_for_pdf"]))
        title = m.get("title", "")
    except Exception:  # noqa: BLE001
        title = ""
    try:
        w = json.loads(get(f"https://api.openalex.org/works/doi:{urllib.parse.quote(doi)}?mailto={email}"))
        for loc in [w.get("best_oa_location")] + (w.get("locations") or []):
            if loc and loc.get("pdf_url"):
                urls.append(("openalex", loc["pdf_url"]))
        title = title or w.get("display_name", "")
    except Exception:  # noqa: BLE001
        pass
    try:
        s = json.loads(get(f"https://api.semanticscholar.org/graph/v1/paper/DOI:{urllib.parse.quote(doi)}?fields=openAccessPdf,title"))
        if (s.get("openAccessPdf") or {}).get("url"):
            urls.append(("semanticscholar", s["openAccessPdf"]["url"]))
        title = title or s.get("title", "")
    except Exception:  # noqa: BLE001
        pass
    seen, out = set(), []
    for src, u in urls:
        if u not in seen:
            seen.add(u); out.append((src, u))
    return out, title


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dois-file", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--email", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    got, missing = [], []
    for doi in [l.strip().lower() for l in open(a.dois_file) if l.strip() and not l.startswith("#")]:
        dest = out / (re.sub(r"[^A-Za-z0-9._-]+", "_", doi) + ".pdf")
        if dest.exists():
            got.append((doi, "already have", str(dest))); continue
        cands, title = candidates(doi, a.email)
        ok = False
        for src, url in cands:
            try:
                data = get(url, timeout=60)
            except Exception:  # noqa: BLE001
                continue
            if data.startswith(b"%PDF"):
                dest.write_bytes(data); got.append((doi, src, str(dest))); ok = True
                break
        if not ok:
            missing.append((doi, title, "no open copy found" if not cands else "open links did not return a PDF"))
        time.sleep(0.3)
    lines = ["# Fetched (open access)"] + [f"- {d} via {s} -> {p}" for d, s, p in got]
    lines += ["", "# Retrieval checklist: download through Caltech or Cornell access, save anywhere, then run ingest_pdfs", ""]
    for pref in sorted({d.split("/")[0] for d, _, _ in missing}):
        lines.append(f"## {PUBLISHERS.get(pref, pref)}")
        lines += [f"- [ ] {t or d} | https://doi.org/{d}" for d, t, _ in missing if d.startswith(pref + "/")]
        lines.append("")
    (out / "RETRIEVAL_CHECKLIST.md").write_text("\n".join(lines) + "\n")
    (out / "missing_dois.txt").write_text("\n".join(d for d, _, _ in missing) + "\n")
    print(f"open-access copies: {len(got)}; still needed through your access: {len(missing)} -> {out}/RETRIEVAL_CHECKLIST.md")


if __name__ == "__main__":
    main()
