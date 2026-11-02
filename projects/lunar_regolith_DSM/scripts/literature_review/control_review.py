#!/usr/bin/env python3
"""The user's review of the controls (gold statements) the checker flagged in a calibration run (WORKFLOW.md §5).

    control_review.py --run <calibration run folder> --set tuning

Writes <run>/untrusted/control_review_<set>.md for the user to open in an editor of their own: each flagged
control's gold statement and page, the checker's verdict, its corrections and its reason, with a line to mark
"gold wrong" (fix the gold statement) or "checker wrong" (over-flagging). The checker's reasons are model text
about the papers, so the file stays under untrusted/ and a normal session never reads it. Prints counts only.
Only the tuning set: held-out results are never shown item by item (calibrate.py).
"""
import argparse
import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_calibration as cc  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", choices=["tuning"], default="tuning")
    a = ap.parse_args(argv)
    run = Path(a.run)
    gold = yaml.safe_load(open(Path(cc.CAL) / "gold_statements.yaml"))["statements"]
    flagged = [r for r in map(json.loads, (run / f"plants_{a.set}.jsonl").read_text().splitlines())
               if "control_id" in r and r.get("flagged")]
    lines = [f"# Controls the checker flagged ({a.set} set, {run.name})", "",
             "For each: is the gold statement wrong (fix it in gold_statements.yaml), or is the checker wrong",
             "(over-flagging)? Write your ruling after 'Ruling:'.", ""]
    for r in flagged:
        g = gold[int(r["control_id"][1:]) - 1]
        vf = run / "untrusted" / r["paper"] / f"plant_verdicts_{a.set}.json"
        v = json.loads(vf.read_text()).get(r["control_id"], {}) if vf.exists() else {}
        lines += [f"## {r['control_id']} ({r['paper']}, printed p. {g.get('page')})", "",
                  f"**Gold statement:** {g.get('statement')}", "",
                  f"**Checker verdict:** {v.get('text')}; equation {v.get('equation')}; parameter {v.get('parameter')}"]
        for k in ("corrected_text", "corrected_latex", "corrected_parameter"):
            if v.get(k):
                lines.append(f"**{k.replace('_', ' ').capitalize()}:** {v[k]}")
        if v.get("dropped_conditions"):
            lines.append(f"**Dropped conditions:** {json.dumps(v['dropped_conditions'], ensure_ascii=False)}")
        lines += [f"**Reason:** {v.get('reason')}", "", "Ruling: ", ""]
    out = run / "untrusted" / f"control_review_{a.set}.md"
    out.write_text("\n".join(lines))
    print(json.dumps({"flagged_controls": len(flagged), "file": str(out)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
