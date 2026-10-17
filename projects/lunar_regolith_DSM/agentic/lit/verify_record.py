#!/usr/bin/env python3
"""Deterministic verification of a paper record against the paper text. No LLM.

Checks (course Lecture 5 'how do you know the agent is right': assert bounds, recompute, external signal):
  1. every quote is verbatim on the cited page (or an adjacent page; flagged otherwise), after normalising case,
     hyphenation and punctuation;  2. quote length <= 25 words;  3. cited pages exist;  4. page coverage;
  5. numbers and equations on pages whose text layer is unreliable are queued for an IMAGE check;
  6. optional DOI check against Crossref (title similarity).
The text layer drops minus signs and exponents, so a quote match proves the words, not the signs: every claim that
contains a digit on a flagged page is queued for the image check as well.
"""
from __future__ import annotations

import difflib
import json
import re
import urllib.parse
import urllib.request

NORM = re.compile(r"[^0-9a-zа-яͰ-Ͽ]+")


def norm(s: str) -> str:
    return NORM.sub("", (s or "").lower())


def split_pages(text: str) -> dict[int, str]:
    parts = re.split(r"=== PAGE (\d+) ===\n", text)
    return {int(parts[i]): parts[i + 1] for i in range(1, len(parts), 2)}


def locate(quote: str, page: int, pages_norm: dict[int, str]) -> tuple[str, int | None]:
    q = norm(quote)
    if len(q) < 12:
        return "quote_too_short", None
    if q in pages_norm.get(page, ""):
        return "page_exact", page
    for p in (page - 1, page + 1):
        if q in pages_norm.get(p, ""):
            return "page_adjacent", p
    for p, t in pages_norm.items():
        if q in t:
            return "elsewhere", p
    return "not_found", None


def check_doi(doi: str, title: str) -> dict:
    try:
        url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
        req = urllib.request.Request(url, headers={"User-Agent": "lit-verify/0.1 (mailto:jbaglino@caltech.edu)"})
        m = json.load(urllib.request.urlopen(req, timeout=20))["message"]
        t = (m.get("title") or [""])[0]
        sim = difflib.SequenceMatcher(None, norm(t), norm(title)).ratio()
        return {"resolves": True, "crossref_title": t, "title_similarity": round(sim, 3), "ok": sim > 0.8}
    except Exception as e:  # noqa: BLE001
        return {"resolves": False, "ok": False, "error": str(e)[:120]}


def verify(record: dict, text: str, pages_meta: dict, do_doi: bool = True) -> dict:
    pages = split_pages(text)
    pn = {p: norm(t) for p, t in pages.items()}
    n_pages = pages_meta["n_pages"]
    flagged = {p["page"]: p["flags"] for p in pages_meta["pages"] if {"garbled_glyphs", "math_dense", "table_like"} & set(p["flags"])}
    res = {"claims": [], "probes": [], "image_queue": [], "problems": []}
    for c in record.get("claims", []):
        status, found = locate(c["quote"], c["page"], pn)
        words = len(c["quote"].split())
        r = {"id": c["id"], "page": c["page"], "status": status, "found_on_page": found, "quote_words": words,
             "over_25_words": words > 25, "has_digit": bool(re.search(r"\d", c["claim"] + c["quote"])),
             "page_in_range": 1 <= c["page"] <= n_pages}
        if (r["has_digit"] or c["type"] in ("number", "equation")) and (found or c["page"]) in flagged:
            r["needs_image"] = True
            res["image_queue"].append({"what": f"claim {c['id']}", "page": found or c["page"]})
        res["claims"].append(r)
    for q in record.get("probes", []):
        if q.get("status") in ("supported", "contradicted") and q.get("quote"):
            status, found = locate(q["quote"], q.get("page", 0), pn)
            res["probes"].append({"question": q["question"][:60], "status": status, "found_on_page": found})
        elif q.get("status") in ("supported", "contradicted"):
            res["probes"].append({"question": q["question"][:60], "status": "no_quote", "found_on_page": None})
    for e in record.get("equations", []):
        if (not e.get("text_layer_ok", True)) or e["page"] in flagged:
            res["image_queue"].append({"what": f"equation {e['label']}", "page": e["page"]})
    for n in record.get("numbers", []):
        if (not n.get("text_layer_ok", True)) or n["page"] in flagged:
            res["image_queue"].append({"what": f"number {n['quantity']}={n['value']} {n['unit']}", "page": n["page"]})
    cl = res["claims"]
    ok = [c for c in cl if c["status"] in ("page_exact", "page_adjacent")]
    res["metrics"] = {
        "n_claims": len(cl), "quote_match_rate": round(len(ok) / len(cl), 3) if cl else None,
        "exact_page_rate": round(sum(c["status"] == "page_exact" for c in cl) / len(cl), 3) if cl else None,
        "not_found": sum(c["status"] == "not_found" for c in cl), "elsewhere": sum(c["status"] == "elsewhere" for c in cl),
        "over_25_words": sum(c["over_25_words"] for c in cl), "pages_out_of_range": sum(not c["page_in_range"] for c in cl),
        "coverage": round(len(set(record.get("pages_read", [])) & set(range(1, n_pages + 1))) / n_pages, 3),
        "probe_quote_match_rate": round(sum(p["status"] in ("page_exact", "page_adjacent", "elsewhere") for p in res["probes"]) / len(res["probes"]), 3) if res["probes"] else None,
        "image_queue_items": len(res["image_queue"])}
    if res["metrics"]["coverage"] < 1.0:
        res["problems"].append(f"pages_read covers {res['metrics']['coverage']:.0%} of the paper")
    if do_doi and record.get("paper", {}).get("doi"):
        res["doi"] = check_doi(record["paper"]["doi"], record["paper"].get("title", ""))
        if not res["doi"]["ok"]:
            res["problems"].append("DOI did not resolve or the title does not match Crossref")
    return res
