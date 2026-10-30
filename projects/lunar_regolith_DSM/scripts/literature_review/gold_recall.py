#!/usr/bin/env python3
"""Calibration diagnostic: how many of the user's gold statements a reading found (WORKFLOW.md §5).

    gold_recall.py --run <run folder> --paper-id <id>

A fresh Opus judge compares the gold statements (calibration/gold_statements.yaml, the
user's own words) with the reader's claims, text only. Writes <run>/untrusted/<id>/recall.json
and prints, per gold statement, its page, kind and result, with no claim text.
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call as ac  # noqa: E402

GOLD = ac.PROJECT / "workflows" / "literature-review" / "calibration" / "gold_statements.yaml"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--paper-id", required=True)
    a = ap.parse_args(argv)
    out = Path(a.run).resolve() / "untrusted" / a.paper_id
    gold = [g for g in yaml.safe_load(open(GOLD))["statements"]
            if g.get("paper") == a.paper_id and g.get("statement")]
    gold_in = [{"index": i, "statement": g["statement"], "value": g.get("value"), "unit": g.get("unit")}
               for i, g in enumerate(gold, 1)]
    claims = json.loads((out / "read.json").read_text())["claims"]
    shown = [{k: c.get(k) for k in ("id", "text", "page_pdf", "equation", "parameter", "conditions", "cited_work")}
             for c in claims]
    with tempfile.TemporaryDirectory() as tmp:
        res, meta = ac.call(tmp, ac.prompt("recall_judge", gold=json.dumps(gold_in, indent=1),
                                           claims=json.dumps(shown)), "recall", "beth-checker")
    if res is None:
        raise SystemExit(f"judge failed: {meta}")
    (out / "recall.json").write_text(json.dumps(res, indent=1))
    by = {r["index"]: r for r in res["results"]}
    n = {"full": 0, "partial": 0, "none": 0}
    for i, g in enumerate(gold, 1):
        r = by.get(i, {"found": "none"})
        n[r["found"]] += 1
        print(f"  {i:2d} p.{g.get('page')} {g.get('kind'):<10} {r['found']}")
    print(f"{a.paper_id}: {len(gold)} gold statements, full {n['full']}, partial {n['partial']}, none {n['none']}, "
          f"claims {len(claims)}")


if __name__ == "__main__":
    main()
