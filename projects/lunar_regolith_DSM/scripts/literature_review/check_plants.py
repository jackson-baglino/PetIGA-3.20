#!/usr/bin/env python3
"""Calibration diagnostic: the checker (E4) on planted errors and on correct controls (WORKFLOW.md §5).

    check_plants.py --run <run folder> --set tuning|heldout [--papers ...]

Each plant (calibrate.load_plants, the side named by --set in split.yaml) is given to the
checker as a claim at its planted page, mixed with controls: the user's gold statements,
which are correct. A plant and the correct statement it came from never share a call, so the
checker cannot spot the plant by comparing the two. The checker sees the cited pages as
2,000 px images and its own applicability checklist of the paper (E4 step 0, made here if
missing).

A plant is caught when its verdict is not "supported", or its equation or value is
corrected, or (for a dropped-condition plant) a dropped condition is named. A control is
flagged on the same rule, without the last clause. Results go to
<run>/plants_<set>.jsonl in calibrate.py's format; for the tuning set a report with
per-kind counts is printed; for the held-out set only the totals, and calibrate.py score
is the one that records it.
"""
import argparse
import collections
import json
import sys
import tempfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call as ac  # noqa: E402
import calibrate as cal  # noqa: E402
import check_calibration as cc  # noqa: E402
import read_check as rc  # noqa: E402

DOCS = Path.home() / "Library" / "Mobile Documents" / "com~apple~CloudDocs" / "Documents"
PER_CALL = 8
PAGES_PER_CALL = 6


def pdf_of(paper):
    papers = {p["id"]: p for p in yaml.safe_load(open(Path(cc.CAL) / "papers.yaml"))}
    return DOCS / papers[paper]["file"].split("#")[0].strip()


def pdf_page(printed, offset):
    p = cal.main_page(printed)
    return int(p) - (offset or 0) if p is not None else None


def items_for(paper, ids, plants, offset):
    """Claims to check: plants on this paper in the set, and the paper's gold statements."""
    out = []
    gold = yaml.safe_load(open(Path(cc.CAL) / "gold_statements.yaml"))["statements"]

    def parameter(symbol, value, unit):
        # A value or unit field travels with the claim, as the reader's parameter field does;
        # plants that change only the field are invisible without it.
        if value is None and not unit:
            return None
        return {"quantity": str(symbol or ""), "value": "" if value is None else str(value),
                "unit_printed": str(unit or "")}

    for p in plants:
        if p["paper"] != paper or p["id"] not in ids:
            continue
        fact = ("S", p["from_statement"]) if p.get("from_statement") is not None else ("U", p["id"])
        g = gold[p["from_statement"] - 1] if p.get("from_statement") else {}
        out.append({"id": p["id"], "role": "plant", "kind_planted": p["kind"], "fact": fact,
                    "text": p["planted_claim"], "page_pdf": pdf_page(p.get("planted_page"), offset),
                    "parameter": parameter(g.get("symbol"), p.get("planted_value", p.get("correct_value")),
                                           p.get("planted_unit", p.get("correct_unit")))})
    for i, g in enumerate(gold, 1):
        if g.get("paper") == paper and g.get("statement"):
            out.append({"id": f"G{i:03d}", "role": "control", "fact": ("S", i),
                        "text": g["statement"], "page_pdf": pdf_page(g.get("page"), offset),
                        "parameter": parameter(g.get("symbol"), g.get("value"), g.get("unit"))})
    return [x for x in out if x["page_pdf"]]


def calls(items):
    """Groups of at most PER_CALL items on at most PAGES_PER_CALL pages, one item per fact."""
    left = sorted(items, key=lambda x: (x["page_pdf"], x["id"]))
    groups = []
    while left:
        g, facts, pages, rest = [], set(), set(), []
        for x in left:
            if (len(g) < PER_CALL and x["fact"] not in facts
                    and (x["page_pdf"] in pages or len(pages) < PAGES_PER_CALL)):
                g.append(x)
                facts.add(x["fact"])
                pages.add(x["page_pdf"])
            else:
                rest.append(x)
        groups.append(g)
        left = rest
    return groups


def caught(v, item):
    hit = v["text"] != "supported" or v["equation"] == "corrected" or v["parameter"] == "corrected"
    if item.get("kind_planted") == "dropped_condition" and v["dropped_conditions"]:
        hit = True
    return hit


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", choices=["tuning", "heldout"], required=True)
    ap.add_argument("--papers", nargs="*")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    split = yaml.safe_load(open(cal.SPLIT))
    ids = set(split["tuning"] if a.set == "tuning" else split["heldout"])
    plants = cal.load_plants()
    papers = a.papers or sorted({p["paper"] for p in plants if p["id"] in ids})
    results, by_item = [], {}
    for paper in papers:
        pdf = pdf_of(paper)
        offset = cc.PageMap(str(pdf)).offset
        out = run / "untrusted" / paper
        out.mkdir(parents=True, exist_ok=True)
        pages = rc.render(pdf, out)
        by_num = {int(f.stem[1:]): f for f in pages}
        if not (out / "checklist.json").exists():
            with tempfile.TemporaryDirectory() as tmp:
                files = rc.folder_with(pages, tmp)
                res, meta = ac.call(tmp, ac.prompt("checker_applicability", files=files, questions=rc.questions()),
                                    "checklist", "beth-checker")
            if res is None:
                raise SystemExit(f"checklist failed for {paper}: {meta}")
            (out / "checklist.json").write_text(json.dumps(res, indent=1))
        checklist = json.loads((out / "checklist.json").read_text())
        items = items_for(paper, ids, plants, offset)
        vfile = out / f"plant_verdicts_{a.set}.json"
        verdicts = json.loads(vfile.read_text()) if vfile.exists() else {}
        for k, g in enumerate(calls(items)):
            if all(x["id"] in verdicts for x in g):
                continue
            shown = [{"id": x["id"], "text": x["text"], "kind": "parameter" if x.get("parameter") else "result",
                      "quote": "", "page_pdf": x["page_pdf"], "location": "", "conditions": [],
                      "quote_check": "not given", **({"parameter": x["parameter"]} if x.get("parameter") else {})}
                     for x in g]
            show = [by_num[p] for p in sorted({x["page_pdf"] for x in g}) if p in by_num]
            with tempfile.TemporaryDirectory() as tmp:
                files = rc.folder_with(show, tmp)
                res, meta = ac.call(tmp, ac.prompt("checker_verdicts", files=files,
                                                   checklist=json.dumps(checklist["items"]),
                                                   claims=json.dumps(shown)), "verdicts", "beth-checker")
            if res is None:
                print(f"{paper} call {k + 1}: failed ({meta.get('exit')})", flush=True)
                continue
            for v in res["verdicts"]:
                verdicts[v["claim_id"]] = v
            vfile.write_text(json.dumps(verdicts, indent=1))
            print(f"{paper} call {k + 1}: {len(g)} items, {meta['seconds']} s", flush=True)
        for x in items:
            v = verdicts.get(x["id"])
            if x["role"] == "plant":
                results.append({"plant_id": x["id"], "caught_by": ("E4" if v and caught(v, x) else None)})
                by_item[x["id"]] = (x, v)
            else:
                results.append({"control_id": x["id"], "flagged": bool(v and caught(v, x)), "paper": paper,
                                "verdict": v["text"] if v else None})
    with open(run / f"plants_{a.set}.jsonl", "w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    rep = cal.score(plants, sorted(i for i in ids if any(r.get("plant_id") == i for r in results)), results)
    ctl = [r for r in results if "control_id" in r]
    if a.set == "tuning":
        print(json.dumps({"per_kind": rep["per_kind"], "by_source": rep["by_source"],
                          "controls_flagged": sum(r["flagged"] for r in ctl), "controls": len(ctl),
                          "control_verdicts": dict(collections.Counter(r["verdict"] for r in ctl))}, indent=1))
    else:
        e4 = [r for r in results if "plant_id" in r]
        print(f"held-out: {sum(bool(r['caught_by']) for r in e4)} of {len(e4)} caught; "
              f"controls flagged {sum(r['flagged'] for r in ctl)} of {len(ctl)}")


if __name__ == "__main__":
    main()
