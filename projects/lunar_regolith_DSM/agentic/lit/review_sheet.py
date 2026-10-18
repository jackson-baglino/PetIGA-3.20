#!/usr/bin/env python3
"""Build the human review sheet: ONLY the items an automated check could not settle, one block per record.
You do not re-read papers. For each item decide: keep / reword the claim to what the quote says / drop it.

Items: (1) claims the independent checker rated partially supported (claim says more than its quote) or unsupported;
(2) equations/numbers the image check corrected or could not read; (3) probes on which an independent replicate disagreed.
  python -m agentic.lit.review_sheet  ->  $LIT_REVIEW_DIR/REVIEW_SHEET.md
"""
from __future__ import annotations

import json

from agentic.lit import common
from agentic.lit.consistency import compare

HEAD = """# Review sheet (gate 2): what to check and why

**Purpose.** The pipeline already read every paper in full, string-matched every quote against its page, and had a second model
tier check a sample of claims. What remains are the places where a program or the checker disagreed with the worker. A claim
that goes beyond its quote could later be cited or used in a derivation as if the paper said it; an equation or number
corrected from the page image could change a result. You do not need to re-read the papers.

**For each item, mark one:** `[k]` keep as written, `[r]` reword the claim to what the quote supports, `[d]` drop it, `[p]` look at the paper page.
Items are ordered by how much they could matter: corrections to equations and numbers first, then unsupported or partial claims,
then replicate dissent. A record with no items is ready to use.

"""


def main():
    rd = common.REVIEW / "records"
    out = [HEAD]
    n_items = 0
    for vf in sorted(rd.glob("*.verification.json")):
        stem = vf.name[: -len(".verification.json")]
        if "_rep2" in stem:
            continue
        ver = json.loads(vf.read_text()); rec = json.loads((rd / f"{stem}.json").read_text())
        block = []
        ic = ver.get("image_check")
        for c in (ic or {}).get("corrections", []):
            if c["was"].startswith("(empty"):      # extractor left it empty: a transcription, not a disagreement
                continue
            block.append(f"- [ ] **IMAGE CORRECTION** p. {c['page']} {c['id']}: extractor `{c['was'][:90]}` -> page shows `{c['now'][:140]}` ({c['note']})")
        for u in (ic or {}).get("unreadable", []):
            block.append(f"- [ ] **UNREADABLE** p. {u['page']} {u['id']}: {u['note']}")
        byid = {c["id"]: c for c in rec["claims"]}
        for cid, v in ver.get("checker", {}).get("flagged", {}).items():
            c = byid.get(cid)
            if c:
                block.append(f"- [ ] **{v['verdict'].upper()}** {cid} (p. {c['page']}): claim: {c['claim']}  | quote: \"{c['quote']}\"  | checker: {v['problem']}")
        rep = rd / f"{stem}_rep2.json"
        if rep.exists():
            for d in compare(rec, json.loads(rep.read_text()))["probe_dissent"]:
                block.append(f"- [ ] **REPLICATE DISSENT** {d['qid']} {d['question']}: run 1 `{d['a']}`, run 2 `{d['b']}`")
        out.append(f"## {stem}: {rec['paper']['title'][:90]} ({rec['paper']['year']})")
        out += block if block else ["- nothing to review"]
        out.append("")
        n_items += len(block)
    (common.REVIEW / "REVIEW_SHEET.md").write_text("\n".join(out) + "\n")
    print(f"{n_items} items -> {common.REVIEW}/REVIEW_SHEET.md")


if __name__ == "__main__":
    main()
