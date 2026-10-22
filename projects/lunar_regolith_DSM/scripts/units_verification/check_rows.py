#!/usr/bin/env python3
"""S4 and S4m of units-verification: beth-unit-checker on disputed rows and on every row's meaning.

    python scripts/units_verification/check_rows.py --run-id <id> --task S4    # each DISAGREE / NOT_FOUND row
    python scripts/units_verification/check_rows.py --run-id <id> --task S4m   # meanings, one call per cited page

S4 gives the checker, for one row, the page images (cited page, the page the reader found it on, and up to two
pages where Marker's index finds it), the table's claim, the blind read, Marker's lines around the symbol and
S3's reasons. Its answer is then checked by script: the unit it reports must parse, and must equal the table's
unit for "table_right" and the read's for "read_right"; otherwise the verdict is marked inconsistent.
S4m gives the checker one page and, for each row cited on it, the table's and the read's printed symbol and
meaning. Anything but "matches" makes the row disputed.
Writes verdicts.json (keys S4 and S4m, one entry per symbol) and the manifest entry.
"""
import argparse
import json
import sys
from pathlib import Path

import agents
import compare_units
import render_pages
import unitparse as up
import uvlib

AGENT = "beth-unit-checker"
NS = {"type": ["string", "null"]}
S4_SCHEMA = {
    "type": "object", "additionalProperties": False,
    "required": ["symbol", "verdict", "printed_as", "unit", "value", "meaning", "page", "location", "reason"],
    "properties": {"symbol": {"type": "string"},
                   "verdict": {"type": "string", "enum": ["table_right", "read_right", "both_wrong", "cannot_tell"]},
                   "printed_as": NS, "unit": NS, "value": NS, "meaning": NS, "page": {"type": ["integer", "null"]},
                   "location": NS, "reason": {"type": "string"}}}
S4M_SCHEMA = {
    "type": "object", "additionalProperties": False, "required": ["rows"],
    "properties": {"rows": {"type": "array", "items": {
        "type": "object", "additionalProperties": False,
        "required": ["symbol", "meaning_verdict", "printed_as_on_page", "page_wording", "location", "reason"],
        "properties": {"symbol": {"type": "string"},
                       "meaning_verdict": {"type": "string", "enum": ["matches", "table_wrong", "cannot_tell"]},
                       "printed_as_on_page": NS, "page_wording": NS, "location": NS,
                       "reason": {"type": "string"}}}}}}
MARKER_EXCERPT = 700


def marker_excerpt(mdir, source, page, forms):
    """The Marker lines on a page that contain the symbol, trimmed; '' if none."""
    f = Path(mdir) / f"{source}_p{page}.md"
    if not f.exists():
        return ""
    flat = up.flatten_markdown(f.read_text())
    lines = []
    for form in forms:
        for start, _ in up.find_symbol(flat, form):
            a = flat.text.rfind("\n", 0, start) + 1
            b = flat.text.find("\n", start)
            line = flat.text[a:b if b >= 0 else None].strip()
            if line not in lines:
                lines.append(line[:400])
    return "\n".join(lines)[:MARKER_EXCERPT]


def s4_pages(c, index):
    pages = [c["page"]]
    if c.get("read", {}).get("page") not in (None, c["page"]):
        pages.append(c["read"]["page"])
    near = sorted((p for p in index.get(c["symbol"], {}).get("pages", []) if p not in pages),
                  key=lambda p: abs(p - c["page"]))
    return pages + near[:2]


def s4_consistency(c, table_unit):
    """Script check of the checker's answer (the 'correction must re-pass S3' rule)."""
    def check(out):
        if out["symbol"] != c["symbol"]:
            return [f"answered for {out['symbol']!r}, asked {c['symbol']!r}"]
        return []
    return check


def judge_s4(out, c, table_unit):
    """'consistent' or the reason the checker's verdict does not hold up against its own reported unit."""
    if out["verdict"] == "cannot_tell":
        return "consistent"
    try:
        u = up.parse_unit(out["unit"]) if out["unit"] is not None else None
    except up.UnitError as e:
        return f"inconsistent: checker's unit does not parse ({e})"
    if out["verdict"] == "table_right":
        if u is None or table_unit is None or not up.same_unit(u, table_unit):
            return "inconsistent: says table_right but its unit differs from the table's"
    if out["verdict"] == "read_right":
        try:
            r = up.parse_unit(c.get("read", {}).get("unit_as_printed"))
        except up.UnitError:
            r = None
        if u is None or r is None or not up.same_unit(u, r):
            return "inconsistent: says read_right but its unit differs from the read's"
    return "consistent"


def run_s4(rdir, rows, comp, index, stage="S4", run_id=None):
    by_sym = {r["symbol"]: r for r in rows}
    out = {}
    disputed = [c for c in comp["rows"] if c["verdict"] in ("DISAGREE", "NOT_FOUND")]
    for c in disputed:
        row = by_sym[c["symbol"]]
        pages = s4_pages(c, index)
        missing = [p for p in pages if not (rdir / "pages" / f"{c['source']}_p{p}.png").exists()]
        if missing and run_id:
            render_pages.main(["--run-id", run_id, *sum((["--add", f"{c['source']}:{p}"] for p in missing), [])])
        pages = [p for p in pages if (rdir / "pages" / f"{c['source']}_p{p}.png").exists()]
        images = {f"page_{p}.png": rdir / "pages" / f"{c['source']}_p{p}.png" for p in pages}
        forms = compare_units.printed_forms(row, c.get("read"))
        excerpts = {p: marker_excerpt(rdir / "marker", c["source"], p, forms) for p in pages}
        task = ("Task S4: adjudicate one disputed row.\n"
                f"Page images: {', '.join(f'page_{p}.png (page {p})' for p in pages)}. Page {c['page']} is the "
                "page the table cites.\n\n"
                + json.dumps({"table": {"symbol": c["symbol"], "printed_as": row["printed"], "meaning": row["meaning"],
                                        "unit": row["unit"], "page": row["page"]},
                              "blind_read": c.get("read") or "the reader did not find the symbol on the cited page",
                              "marker_lines_by_page": {str(p): e for p, e in excerpts.items() if e},
                              "script_reasons": c["reasons"]}, ensure_ascii=False, indent=1))
        t_unit, _ = compare_units.parse_or_none(row["unit"])
        try:
            res, info = agents.call_agent(rdir, stage, AGENT, f"{c['symbol']}", task, images, S4_SCHEMA,
                                          check=s4_consistency(c, t_unit))
            res["consistency"] = judge_s4(res, c, t_unit)
            res["s3_verdict"], res["s3_reasons"], res["images"] = c["verdict"], c["reasons"], pages
            res["call"] = info
        except agents.AgentError as e:
            res = {"symbol": c["symbol"], "verdict": "call_failed", "error": json.loads(str(e))}
        out[c["symbol"]] = res
        print(f"  S4 {c['symbol']:<11} {res['verdict']:<12} {res.get('consistency', '')}", flush=True)
    return out


def run_s4m(rdir, rows, comp, stage="S4m"):
    by_sym = {c["symbol"]: c for c in comp["rows"]}
    groups = {}
    for r in rows:
        groups.setdefault((r["source"], r["page"]), []).append(r)
    out = {}
    for (src, pg), grp in sorted(groups.items()):
        img = rdir / "pages" / f"{src}_p{pg}.png"
        items = []
        for r in grp:
            rd = by_sym.get(r["symbol"], {}).get("read") or {}
            items.append({"symbol": r["symbol"], "table_printed_as": r["printed"], "table_meaning": r["meaning"],
                          "read_printed_as": rd.get("printed_as"), "read_meaning": rd.get("meaning")})
        task = ("Task S4m: review meanings on one page.\nPage image: page.png.\nRows cited on this page:\n"
                + json.dumps(items, ensure_ascii=False, indent=1))
        want = sorted(r["symbol"] for r in grp)

        def check(o, want=want):
            got = sorted(x["symbol"] for x in o["rows"])
            return [] if got == want else [f"answered {got}, asked {want}"]
        try:
            res, info = agents.call_agent(rdir, stage, AGENT, f"{src}_p{pg}", task, {"page.png": img}, S4M_SCHEMA,
                                          check=check)
            for x in res["rows"]:
                out[x["symbol"]] = dict(x, page=pg, call=info)
        except agents.AgentError as e:
            for r in grp:
                out[r["symbol"]] = {"symbol": r["symbol"], "meaning_verdict": "call_failed",
                                    "error": json.loads(str(e))}
        print(f"  S4m {src} p{pg}: " + ", ".join(f"{r['symbol']}={out[r['symbol']]['meaning_verdict']}"
                                                for r in grp), flush=True)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--task", required=True, choices=["S4", "S4m"])
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    print(f"{a.task} check_rows: run={uvlib.rel(rdir)} table={uvlib.rel(a.table)}")
    try:
        uvlib.require_committed(a.allow_dirty)
        rows = uvlib.literature_rows(uvlib.read_table(a.table))
        comp = json.loads((rdir / "compare.json").read_text())
        pending = [c["symbol"] for c in comp["rows"] if c["verdict"] == "NOT_FOUND" and c["reread_pages"]]
        if a.task == "S4" and pending:
            raise uvlib.StageError(f"re-reads still pending for {pending}; run read_pages.py --rereads, then S3")
        index = json.loads((rdir / "marker_index.json").read_text())
        res = run_s4(rdir, rows, comp, index, run_id=a.run_id) if a.task == "S4" else run_s4m(rdir, rows, comp)
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    vfile = rdir / "verdicts.json"
    verdicts = json.loads(vfile.read_text()) if vfile.exists() else {}
    verdicts[a.task] = res
    uvlib.write_json(vfile, verdicts)
    agent = agents.load_agent(AGENT)
    key = "verdict" if a.task == "S4" else "meaning_verdict"
    counts = {}
    for v in res.values():
        counts[v[key]] = counts.get(v[key], 0) + 1
    uvlib.update_manifest(rdir, a.task, {"agent": AGENT, "model": agent["model"], "agent_sha256": agent["sha256"],
                                         "compare_sha256": uvlib.sha256(rdir / "compare.json"), "counts": counts,
                                         "verdict": "pass"}, a.allow_dirty)
    print(f"{a.task}: {counts} -> {uvlib.rel(vfile)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
