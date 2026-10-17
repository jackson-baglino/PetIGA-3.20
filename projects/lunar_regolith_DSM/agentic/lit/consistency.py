#!/usr/bin/env python3
"""Agreement between two independent extractions of the same paper (course Lecture 4: vote only over independent samples of a
comparable quantity, with tolerance, and log the dissent). Reports; it does not decide which record is right."""
from __future__ import annotations

import re

from agentic.lit.verify_record import norm


def _digits(s: str) -> str:
    return "".join(re.findall(r"[0-9.]+", s or ""))


def compare(a: dict, b: dict) -> dict:
    pa = {p["qid"]: p for p in a.get("probes", []) if p.get("qid")}
    pb = {p["qid"]: p for p in b.get("probes", []) if p.get("qid")}
    common = sorted(set(pa) & set(pb), key=lambda x: int(x[1:]))
    agree = [q for q in common if pa[q]["status"] == pb[q]["status"]]
    ea = {_digits(e["label"]) for e in a.get("equations", [])}
    eb = {_digits(e["label"]) for e in b.get("equations", [])}
    na = {_digits(n["value"]) for n in a.get("numbers", []) if _digits(n["value"])}
    nb = {_digits(n["value"]) for n in b.get("numbers", []) if _digits(n["value"])}
    qa = {norm(c["quote"]) for c in a.get("claims", [])}
    qb = {norm(c["quote"]) for c in b.get("claims", [])}
    j = lambda x, y: round(len(x & y) / len(x | y), 3) if (x | y) else None
    return {"probe_status_agreement": round(len(agree) / len(common), 3) if common else None, "n_probes_compared": len(common),
            "probe_dissent": [{"qid": q, "question": pa[q]["question"][:70], "a": pa[q]["status"], "b": pb[q]["status"]} for q in common if pa[q]["status"] != pb[q]["status"]],
            "equation_number_jaccard": j(ea, eb), "number_value_jaccard": j(na, nb), "claim_quote_jaccard": j(qa, qb),
            "n_claims": [len(a.get("claims", [])), len(b.get("claims", []))]}
