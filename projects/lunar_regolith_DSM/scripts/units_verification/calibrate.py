#!/usr/bin/env python3
"""C0 of units-verification: measure the pipeline on planted errors before any real run.

    python scripts/units_verification/calibrate.py --run-id <id>

Builds <run>/calibration/units_planted.md: 14 rows of the real table, 10 with a planted error (PLANTS) and 4
untouched controls (CONTROLS), then runs S2 -> S3 -> re-reads -> S3 -> S4 -> S4m on them, re-using the run's page
images and Marker reading (S0, S1 and S1b must have passed). Scores (WORKFLOW.md §4):
  S3 flags 8/8 unit and page plants; S4 rejects >= 7/8 of them (says the table is wrong and gives the original
  unit or page); S4m calls 2/2 meaning plants table_wrong; every false alarm on the controls is reported.
Writes <run>/calibration.json and the C0 entry of the manifest.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import agents
import check_rows
import compare_units
import marker_extract
import read_pages
import render_pages
import unitparse as up
import uvlib

# (symbol, kind, field, planted value). Fixed so the calibration is repeatable; kinds from WORKFLOW.md §4.
PLANTS = [
    ("beta", "flipped exponent sign", "unit", "s m^1"),
    ("L_sg", "flipped exponent sign", "unit", "J m^3"),
    ("a", "wrong prefix", "unit", "mm"),
    ("gamma_ia", "wrong prefix", "unit", "mN m^-1"),
    ("kappa_i", "swapped unit", "unit", "J m^-1 K^-1"),
    ("c_p_i", "swapped unit", "unit", "J m^-3 K^-1"),
    ("k", "wrong page", "page", 13),         # Marker finds k only on K&P p4
    ("K_w", "wrong page", "page", 7),        # Marker finds K_w only on SI p9
    ("rho_i", "swapped meaning", "meaning", "density of dry air (printed as ρi)"),
    ("kappa_a", "swapped meaning", "meaning", "heat conductivity of ice (printed as κa)"),
]
CONTROLS = ["gamma", "h", "u_j", "Lambda"]


def planted_table(real_rows, header):
    by = {r["symbol"]: r for r in real_rows}
    rows, plants = [], []
    for sym, kind, field, val in PLANTS:
        r = dict(by[sym])
        plants.append({"symbol": sym, "kind": kind, "field": field, "original": r[field], "planted": val})
        r[field] = val
        rows.append(r)
    rows += [dict(by[s]) for s in CONTROLS]
    lines = [header]
    for r in rows:
        lines.append("| " + " | ".join([f"`{r['symbol']}`", r["meaning"], r["unit"], r["source"],
                                        str(r["page"]), r["quote"], r["verified"]]) + " |")
    return "\n".join(lines) + "\n", plants


def score(plants, comp, s4, s4m):
    by_c = {c["symbol"]: c for c in comp["rows"]}
    unit_page = [p for p in plants if p["field"] in ("unit", "page")]
    meaning = [p for p in plants if p["field"] == "meaning"]
    detail = []
    s3_hits = s4_hits = s4m_hits = 0
    for p in unit_page:
        c, v = by_c[p["symbol"]], s4.get(p["symbol"], {})
        flagged = c["verdict"] in ("DISAGREE", "NOT_FOUND")
        rejected = v.get("verdict") in ("read_right", "both_wrong")
        if rejected and p["field"] == "unit":
            try:
                rejected = v.get("unit") is not None and up.same_unit(up.parse_unit(v["unit"]),
                                                                       up.parse_unit(p["original"]))
            except up.UnitError:
                rejected = False
        if rejected and p["field"] == "page":
            rejected = v.get("page") == p["original"]
        s3_hits += flagged
        s4_hits += rejected
        detail.append(dict(p, s3=c["verdict"], s3_reasons=c["reasons"], s4=v.get("verdict"), s4_unit=v.get("unit"),
                           s4_page=v.get("page"), s4_consistency=v.get("consistency"), caught_s3=flagged,
                           caught_s4=rejected))
    for p in meaning:
        v = s4m.get(p["symbol"], {})
        hit = v.get("meaning_verdict") == "table_wrong"
        s4m_hits += hit
        detail.append(dict(p, s3=by_c[p["symbol"]]["verdict"], s3_reasons=by_c[p["symbol"]]["reasons"],
                           s4m=v.get("meaning_verdict"), s4m_wording=v.get("page_wording"), caught_s4m=hit))
    controls = []
    for sym in CONTROLS:
        c, v, m = by_c[sym], s4.get(sym), s4m.get(sym, {})
        alarms = []
        if c["verdict"] != "AGREE":
            alarms.append(f"S3 {c['verdict']}: {'; '.join(c['reasons'])}")
            if v and v.get("verdict") != "table_right":
                alarms.append(f"S4 {v.get('verdict')}: {v.get('reason')}")
        if m.get("meaning_verdict") != "matches":
            alarms.append(f"S4m {m.get('meaning_verdict')}: {m.get('reason')}")
        controls.append({"symbol": sym, "s3": c["verdict"], "s4": v.get("verdict") if v else None,
                         "s4m": m.get("meaning_verdict"), "false_alarms": alarms})
    crit = {"S3_unit_page": f"{s3_hits}/{len(unit_page)}", "S4_unit_page": f"{s4_hits}/{len(unit_page)}",
            "S4m_meaning": f"{s4m_hits}/{len(meaning)}"}
    passed = s3_hits == len(unit_page) and s4_hits >= len(unit_page) - 1 and s4m_hits == len(meaning)
    return {"criteria": crit, "pass": passed,
            "control_false_alarms": sum(bool(c["false_alarms"]) for c in controls),
            "plants": detail, "controls": controls}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    cdir = rdir / "calibration"
    print(f"C0 calibrate: run={uvlib.rel(rdir)} table={uvlib.rel(a.table)} plants={len(PLANTS)} "
          f"controls={len(CONTROLS)}")
    try:
        uvlib.require_committed(a.allow_dirty)
        stages = uvlib.load_manifest(rdir).get("stages", {})
        for s in ("S0", "S1", "S1b"):
            if stages.get(s, {}).get("verdict") != "pass":
                raise uvlib.StageError(f"{s} has not passed for this run")
        n_pages = {s: v["pages"] for s, v in stages["S0"]["sources"].items()}
        real = uvlib.literature_rows(uvlib.read_table(a.table))
        header = "| symbol | meaning | unit | source | page | quote (text layer) | verified |\n|---|---|---|---|---|---|---|"
        text, plants = planted_table(real, header)
        cdir.mkdir(exist_ok=True)
        (cdir / "units_planted.md").write_text("# Calibration copy with planted errors (C0); not real data\n\n" + text)
        uvlib.write_json(cdir / "plants.json", {"plants": plants, "controls": CONTROLS})
        for name in ("pages", "marker"):
            if not (cdir / name).exists():
                os.symlink(Path("..") / name, cdir / name)
        rows = uvlib.read_table(cdir / "units_planted.md")
        index = marker_extract.build_index(rows, cdir / "marker", n_pages)
        uvlib.write_json(cdir / "marker_index.json", index)
        need = sorted({(r["source"], p) for r in rows for p in [r["page"]] + index[r["symbol"]]["pages"]})
        missing = [f"{s}:{p}" for s, p in need if not (rdir / "pages" / f"{s}_p{p}.png").exists()]
        if missing:
            render_pages.main(["--run-id", a.run_id, *sum((["--add", m] for m in missing), [])]
                              + (["--allow-dirty"] if a.allow_dirty else []))

        print("C0-S2: blind reads", flush=True)
        _, failed = read_pages.run_groups(cdir, "C0-S2", read_pages.cited_groups(rows))

        def s3():
            reads = compare_units.load_reads(cdir)
            res = [compare_units.compare_row(r, reads, cdir / "marker", index) for r in rows]
            uvlib.write_json(cdir / "compare.json", {"table": "calibration/units_planted.md", "rows": res})
            return {"rows": res}
        comp = s3()
        print("C0-S2: re-reads", flush=True)
        _, failed2 = read_pages.run_groups(cdir, "C0-S2", read_pages.reread_groups(cdir, rows))
        comp = s3()
        print("C0-S3: " + ", ".join(f"{c['symbol']}={c['verdict']}" for c in comp["rows"]), flush=True)
        print("C0-S4: adjudication", flush=True)
        s4 = check_rows.run_s4(cdir, rows, comp, index, stage="C0-S4")
        print("C0-S4m: meaning review", flush=True)
        s4m = check_rows.run_s4m(cdir, rows, comp, stage="C0-S4m")
        uvlib.write_json(cdir / "verdicts.json", {"S4": s4, "S4m": s4m})
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    result = score(plants, comp, s4, s4m)
    result.update({"read_failures": failed + failed2})
    uvlib.write_json(rdir / "calibration.json", result)
    uvlib.update_manifest(rdir, "C0", {
        "agents": {n: agents.load_agent(n)["sha256"] for n in ("morty-unit-reader", "beth-unit-checker")},
        "table_sha256": uvlib.sha256(a.table), "planted_table_sha256": uvlib.sha256(cdir / "units_planted.md"),
        "criteria": result["criteria"], "control_false_alarms": result["control_false_alarms"],
        "verdict": "pass" if result["pass"] else "below threshold"}, a.allow_dirty)
    print(f"C0: {result['criteria']} control false alarms {result['control_false_alarms']}/{len(CONTROLS)} "
          f"-> {'PASS' if result['pass'] else 'BELOW THRESHOLD'}")
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
