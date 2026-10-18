#!/usr/bin/env python3
"""Independent claim check: a separate worker (different model tier, isolated context) judges whether each sampled quote
supports its claim, seeing only the claim, the quote and the cited page text. Sampled: every claim whose quote was not
found on the cited page, every numeric/equation claim, plus a random fraction of the rest (default 25 %).

  python -m agentic.lit.check_claims --record ji_karma_2026 [--model opus] [--fraction 0.25] [--seed 1]
"""
from __future__ import annotations

import argparse
import json
import random
import re

from concurrent.futures import ThreadPoolExecutor

from agentic.lit import common
from agentic.lit.verify_record import split_pages

SCHEMA = {"type": "object", "required": ["verdict", "quote_found_on_page", "problem"], "properties": {
    "verdict": {"type": "string", "enum": ["supported", "partially_supported", "unsupported", "cannot_tell"]},
    "quote_found_on_page": {"type": "boolean"}, "problem": {"type": "string"}}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", required=True, help="record stem in $LIT_REVIEW_DIR/records (slug or slug_tag)")
    ap.add_argument("--model", default=None)
    ap.add_argument("--fraction", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    rd = common.REVIEW / "records"
    rec = json.loads((rd / f"{a.record}.json").read_text())
    ver = json.loads((rd / f"{a.record}.verification.json").read_text())
    slug = re.sub(r"(_(opus|sonnet|haiku|rep\d+))+$", "", a.record)
    pages = split_pages((common.STATE / "text" / f"{slug}.txt").read_text())
    byid = {c["id"]: c for c in rec["claims"]}
    vid = {c["id"]: c for c in ver["claims"]}
    must = [i for i, v in vid.items() if v["status"] not in ("page_exact", "page_adjacent") or byid[i]["type"] in ("number", "equation") or v.get("needs_image")]
    rest = [i for i in vid if i not in must]
    random.Random(a.seed).shuffle(rest)
    chosen = must + rest[: int(round(a.fraction * len(rest)))]

    def one(i):
        c = byid[i]
        page = vid[i]["found_on_page"] or c["page"]
        prompt = f"CLAIM: {c['claim']}\n\nQUOTE: {c['quote']}\n\nCITED PAGE ({page}) TEXT:\n{pages.get(page, '')}"
        out, u = common.call_claude("claim-checker", prompt, model=a.model, schema=SCHEMA, label=f"check:{a.record}:{i}")
        return i, out

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        results = dict(ex.map(one, chosen))
    ver["checker"] = {"model": a.model or common.parse_agent("claim-checker")["model"], "n_checked": len(results), "n_must": len(must),
                      "verdicts": results,
                      "supported_rate": round(sum(v["verdict"] == "supported" for v in results.values()) / len(results), 3) if results else None,
                      "flagged": {i: v for i, v in results.items() if v["verdict"] != "supported"}}
    (rd / f"{a.record}.verification.json").write_text(json.dumps(ver, indent=1))
    k = ver["checker"]
    print(f"{a.record}: checked {k['n_checked']} claims ({k['n_must']} mandatory); supported {k['supported_rate']}; flagged {list(k['flagged'])}")


if __name__ == "__main__":
    main()
