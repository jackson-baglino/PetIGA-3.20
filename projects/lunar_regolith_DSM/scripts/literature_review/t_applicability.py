#!/usr/bin/env python3
"""T-applicability (WORKFLOW.md §12): does the checker's applicability checklist (E4 step 0) reach a condition
printed on a different page from the result, and is the dropped-condition plant caught?

    t_applicability.py --run <calibration run folder> --set tuning|heldout

Reads the run's plants_<set>.jsonl (caught or not) and each paper's checklist (untrusted/<paper>/checklist.json,
made by check_plants.py) and prints, per dropped-condition plant in the set: caught, and whether the checklist
has an item on the condition's page. IDs and yes/no only; never checklist text. Exit 0 if every plant passes both.
"""
import argparse
import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate as cal  # noqa: E402
import check_calibration as cc  # noqa: E402
import check_plants as cp  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", choices=["tuning", "heldout"], required=True)
    a = ap.parse_args(argv)
    run = Path(a.run)
    ids = set(yaml.safe_load(open(cal.SPLIT))[a.set])
    plants = [p for p in cal.load_plants() if p["id"] in ids and p["kind"] == "dropped_condition"]
    caught = {}
    for line in (run / f"plants_{a.set}.jsonl").read_text().splitlines():
        r = json.loads(line)
        if "plant_id" in r:
            caught[r["plant_id"]] = bool(r.get("caught_by"))
    ok_all = True
    for p in plants:
        offset = cc.PageMap(str(cp.pdf_of(p["paper"]))).offset or 0
        cond_pdf = int(cal.main_page(p["condition_page"])) - offset
        cl = run / "untrusted" / p["paper"] / "checklist.json"
        pages = {i.get("page_pdf") for i in json.loads(cl.read_text())["items"]} if cl.exists() else set()
        on_page = cond_pdf in pages
        ok = caught.get(p["id"], False) and on_page
        ok_all &= ok
        print(f"{'PASS' if ok else 'FAIL'}  {p['id']} ({p['paper']}): caught {caught.get(p['id'])}; "
              f"checklist item on the condition page (PDF p. {cond_pdf}): {on_page}")
    if not plants:
        print(f"no dropped-condition plants in the {a.set} set")
        return 1
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
