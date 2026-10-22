#!/usr/bin/env python3
"""Reader-stability check of units-verification (WORKFLOW.md §4): compare two blind reads of the same pages.

    python scripts/units_verification/read_pages.py --run-id <id> --repeat 5 --out-dir reads_repeat --stage S2-repeat
    python scripts/units_verification/reader_stability.py --run-id <id>

For each symbol read in both reads/ and reads_repeat/, compares what S3 relies on: found, printed symbol
(printed_key), unit (dimensions and scale) and printed numbers. Meaning and location are free text and are not
compared. Reports the disagreement count (there is no pass threshold; the design asks for the count).
Writes stability.json and the "stability" entry of the manifest.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import unitparse as up
import uvlib


def unit_of(r):
    try:
        return up.parse_unit(r.get("unit_as_printed")) if r.get("unit_as_printed") is not None else None
    except up.UnitError:
        return "unparsed:" + str(r.get("unit_as_printed"))


def differences(a, b):
    """Fields on which two reads of one symbol differ."""
    diff = []
    if bool(a.get("found")) != bool(b.get("found")):
        return ["found"]
    if not a.get("found"):
        return []
    if up.printed_key(a.get("printed_as")) != up.printed_key(b.get("printed_as")):
        diff.append("printed_as")
    ua, ub = unit_of(a), unit_of(b)
    if isinstance(ua, up.Unit) and isinstance(ub, up.Unit):
        if not up.same_unit(ua, ub):
            diff.append("unit")
    elif ua != ub:
        diff.append("unit")
    na, nb = up.parse_numbers(a.get("value_as_printed")), up.parse_numbers(b.get("value_as_printed"))
    if len(na) != len(nb) or not all(math.isclose(x, y, rel_tol=1e-9) for x, y in zip(na, nb)):
        diff.append("value")
    return diff


def compare_dirs(d1, d2):
    rows = []
    for f2 in sorted(Path(d2).glob("*.json")):
        f1 = Path(d1) / f2.name
        if not f1.exists():
            continue
        r1 = {r["symbol"]: r for r in json.loads(f1.read_text())["reads"]}
        for r in json.loads(f2.read_text())["reads"]:
            if r["symbol"] in r1:
                d = differences(r1[r["symbol"]], r)
                rows.append({"page": f2.stem, "symbol": r["symbol"], "differs_in": d,
                             "first": {k: r1[r["symbol"]].get(k) for k in ("found", "printed_as", "unit_as_printed",
                                                                           "value_as_printed")},
                             "second": {k: r.get(k) for k in ("found", "printed_as", "unit_as_printed",
                                                              "value_as_printed")}})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--repeat-dir", default="reads_repeat")
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    print(f"stability: run={uvlib.rel(rdir)} compare reads/ with {a.repeat_dir}/")
    try:
        uvlib.require_committed(a.allow_dirty)
        if not (rdir / a.repeat_dir).exists():
            raise uvlib.StageError(f"no {a.repeat_dir}/; run read_pages.py --repeat 5 --out-dir {a.repeat_dir}")
        rows = compare_dirs(rdir / "reads", rdir / a.repeat_dir)
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    n_diff = sum(bool(r["differs_in"]) for r in rows)
    pages = sorted({r["page"] for r in rows})
    uvlib.write_json(rdir / "stability.json", {"pages": pages, "symbols": len(rows), "disagreements": n_diff,
                                               "rows": rows})
    uvlib.update_manifest(rdir, "stability", {"pages": pages, "symbols": len(rows), "disagreements": n_diff,
                                              "output_sha256": uvlib.sha256(rdir / "stability.json"),
                                              "verdict": "reported"}, a.allow_dirty)
    print(f"{len(rows)} symbols on {len(pages)} pages read twice: {n_diff} differ")
    for r in rows:
        if r["differs_in"]:
            print(f"  {r['page']} {r['symbol']}: {r['differs_in']}  {r['first']}  vs  {r['second']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
