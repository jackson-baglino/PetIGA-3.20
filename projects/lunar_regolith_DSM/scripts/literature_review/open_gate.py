#!/usr/bin/env python3
"""Open a blocking gate in a B1 run: the line stops there until the user records a decision (run_b1.py approve).

    open_gate.py --run <dir> --stage <id> [--item <paper>] --what <text> [--count N]

Stations call this (or open_gate() from Python) when they need the user's decision, for example A2 for a copy
that is not the version of record, or X1 for reads that disagree. The gate file holds only the stage, the item,
a one-line description and a count, never paper text.
"""
import argparse
import datetime
import json
from pathlib import Path


def open_gate(run, stage, what, item=None, count=None):
    g = Path(run) / "gates"
    g.mkdir(parents=True, exist_ok=True)
    f = g / f"{stage}{'_' + item if item else ''}.json"
    rec = {"stage": stage, "item": item, "status": "waiting", "what": what, "count": count,
           "opened": datetime.datetime.now().astimezone().isoformat(timespec="seconds")}
    f.write_text(json.dumps(rec, indent=1))
    return f


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--item")
    ap.add_argument("--what", required=True)
    ap.add_argument("--count", type=int)
    a = ap.parse_args()
    print(open_gate(a.run, a.stage, a.what, a.item, a.count))
