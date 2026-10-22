#!/usr/bin/env python3
"""S3 of units-verification: compare each table row with its blind read and with Marker's reading.

    python scripts/units_verification/compare_units.py --run-id <id> [--marker-only]

Per literature row of the units table:
  AGREE      the table's unit equals the blind read's unit (same dimensions and scale), the table's printed
             symbol equals the read's (ρi == ρ_i; ρv != ρ_vs), and Marker's text on that page shows the symbol
             followed by the same unit (and the read's printed value, if any)
  DISAGREE   anything else; the reasons say which part failed
  NOT_FOUND  the reader did not find the symbol on the cited page; reread_pages lists the pages where
             Marker's index finds it (S2 re-reads them, then S3 runs again)
  MISSING    no read exists for the row (a pipeline error, not a verdict on the row)
Units are compared as dimension vectors and scale (unitparse.py): J m^-3 K^-1 == J/(m3 K), °C difference == K.
Meanings are not compared here (no deterministic check exists); S4's meaning review sees every row.

Reads are runs/units-verification/<id>/reads/<source>_p<N>.json:
  {"page": "<source>_p<N>", "reads": [{symbol, found, printed_as, meaning, value_as_printed, unit_as_printed,
  location}, ...]}
--marker-only skips the reads and reports only whether Marker shows the table's unit next to the symbol on
the cited page (a check of the Marker side before any agent runs).
Writes compare.json (or compare_marker_only.json) and the S3 entry of manifest.json.
"""
import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

import unitparse as up
import uvlib

WINDOW = 300  # characters after a symbol that may hold its unit, within the same line (table row or paragraph)
# A dimensionless unit is printed as "(-)", "(1)", "(dimensionless)" or "kg/kg"; it counts if no other unit
# comes before it after the symbol (in a table it may sit a cell or two later).
DIMENSIONLESS_MARK = re.compile(r"\(\s*(?:-|1|dimensionless)\s*\)|\bkg\s*/\s*kg\b", re.I)


def printed_forms(row, read=None):
    """Glyph strings to look for in Marker's text: the table's printed form and the reader's, if short."""
    forms = [row["printed"]] if row["printed"] else []
    p = (read or {}).get("printed_as")
    if p:
        p = up.flatten(p).replace("_", "").replace(" ", "")
        if 0 < len(p) <= 8 and p not in forms:
            forms.append(p)
    return forms


def marker_check(flat, forms, table_unit, read_numbers):
    """What Marker's page shows after the symbol: located?, first unit after it, numbers nearby."""
    occ = []
    for f in forms:
        occ += up.find_symbol(flat, f)
    out = {"located": bool(occ), "occurrences": len(occ), "units_seen": [], "unit_match": False,
           "numbers": [], "value_match": None}
    nums = []
    for _, end in sorted(set(occ)):
        nl = flat.text.find("\n", end)
        win = flat.text[end:(nl if nl >= 0 else len(flat.text))][:WINDOW]
        runs = up.unit_runs(win)
        nums += up.parse_numbers(win)
        mark = DIMENSIONLESS_MARK.search(win)
        if (table_unit is not None and table_unit.dims == up.ONE.dims and mark
                and (not runs or runs[0][3] is None or runs[0][3].dims == up.ONE.dims or mark.start() <= runs[0][0])):
            out["units_seen"].append("(dimensionless mark)")
            out["unit_match"] = True
            continue
        if runs:
            raw, u = runs[0][2], runs[0][3]
            out["units_seen"].append(raw)
            if u is not None and table_unit is not None and up.same_unit(u, table_unit):
                out["unit_match"] = True
    out["numbers"] = sorted(set(nums))
    if read_numbers:
        out["value_match"] = up.numbers_contained(read_numbers, nums)
    return out


def parse_or_none(text):
    try:
        return up.parse_unit(text), None
    except up.UnitError as e:
        return None, str(e)


def load_reads(rdir):
    """{(source, page): {symbol: read}} from reads/*.json."""
    out = {}
    for f in sorted((rdir / "reads").glob("*_p*.json")):
        m = re.fullmatch(r"(.+)_p(\d+)\.json", f.name)
        data = json.loads(f.read_text())
        out[(m.group(1), int(m.group(2)))] = {r["symbol"]: r for r in data["reads"]}
    return out


def marker_flat(mdir, source, page):
    f = mdir / f"{source}_p{page}.md"
    return up.flatten_markdown(f.read_text()) if f.exists() else None


def compare_row(row, reads, mdir, index, marker_only=False):
    sym, src, page = row["symbol"], row["source"], row["page"]
    res = {"symbol": sym, "source": src, "page": page, "table_unit": row["unit"], "reasons": [],
           "reread_pages": []}
    t_unit, t_err = parse_or_none(row["unit"])
    res["table_unit_si"] = t_unit.canonical() if t_unit else None
    if t_unit is None:
        res["reasons"].append(f"table_unit_unparsed: {t_err}")

    read, read_page = None, page
    if not marker_only:
        read = reads.get((src, page), {}).get(sym)
        if read is None:
            res["verdict"] = "MISSING"
            res["reasons"].append("no read for this row on its cited page")
            return res
        if not read.get("found"):
            others = [(p, reads[(src, p)][sym]) for p in sorted(index.get(sym, {}).get("pages", []))
                      if p != page and sym in reads.get((src, p), {})]
            hit = next(((p, r) for p, r in others if r.get("found")), None)
            if hit is None:
                pending = [p for p in index.get(sym, {}).get("pages", []) if p != page
                           and sym not in reads.get((src, p), {})]
                res["verdict"] = "NOT_FOUND"
                res["reread_pages"] = pending
                res["reasons"].append(f"reader did not find it on p{page}"
                                      + (f"; Marker finds it on p{pending}" if pending else
                                         "; Marker finds it on no other page" if not others else
                                         f"; re-reads on p{[p for p, _ in others]} did not find it either"))
                return res
            read_page, read = hit
            res["reasons"].append(f"page: found on p{read_page}, not on the cited p{page}")
        r_unit, r_err = parse_or_none(read.get("unit_as_printed"))
        res["read"] = {k: read.get(k) for k in ("printed_as", "meaning", "value_as_printed", "unit_as_printed",
                                                "location")}
        res["read"]["page"] = read_page
        res["read_unit_si"] = r_unit.canonical() if r_unit else None
        if r_unit is None:
            res["reasons"].append(f"read_unit_unparsed: {r_err}")
        elif t_unit is not None and not up.same_unit(t_unit, r_unit):
            res["reasons"].append(f"unit: table {t_unit.canonical()} vs read {r_unit.canonical()}")
        t_key, r_key = up.printed_key(row["printed"]), up.printed_key(read.get("printed_as"))
        if t_key != r_key:
            res["reasons"].append(f"printed: table {row['printed']!r} vs read {read.get('printed_as')!r}")

    flat = marker_flat(mdir, src, read_page)
    read_numbers = up.parse_numbers(read.get("value_as_printed")) if read else []
    if flat is None:
        res["marker"] = None
        res["reasons"].append(f"marker: no Marker reading of p{read_page}")
    else:
        mk = marker_check(flat, printed_forms(row, read), t_unit, read_numbers)
        res["marker"] = mk
        if not mk["located"]:
            res["reasons"].append(f"marker: symbol not found on p{read_page}")
        elif not mk["unit_match"]:
            res["reasons"].append(f"marker: unit after the symbol is {mk['units_seen'] or 'none'}")
        if mk["value_match"] is False:
            res["reasons"].append("marker: printed value not found near the symbol")
    if marker_only:
        res["verdict"] = "MARKER_AGREES" if not res["reasons"] else "MARKER_DIFFERS"
    else:
        res["verdict"] = "AGREE" if not res["reasons"] else "DISAGREE"
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--marker-only", action="store_true")
    ap.add_argument("--stage", default="S3", help="manifest key (C0 uses its own)")
    ap.add_argument("--out", default=None, help="output file name in the run folder")
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    out_name = a.out or ("compare_marker_only.json" if a.marker_only else "compare.json")
    print(f"{a.stage} compare_units: run={uvlib.rel(rdir)} table={uvlib.rel(a.table)} marker_only={a.marker_only}")
    try:
        uvlib.require_committed(a.allow_dirty)
        index_file = rdir / "marker_index.json"
        if not index_file.exists():
            raise uvlib.StageError("no marker_index.json; run marker_extract.py first")
        index = json.loads(index_file.read_text())
        reads = {} if a.marker_only else load_reads(rdir)
        rows = uvlib.literature_rows(uvlib.read_table(a.table))
        results = [compare_row(r, reads, rdir / "marker", index, a.marker_only) for r in rows]
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    counts = dict(Counter(r["verdict"] for r in results))
    uvlib.write_json(rdir / out_name, {"table": uvlib.rel(a.table), "table_sha256": uvlib.sha256(a.table),
                                       "counts": counts, "rows": results})
    read_files = {f.name: uvlib.sha256(f) for f in sorted((rdir / "reads").glob("*.json"))} \
        if (rdir / "reads").exists() and not a.marker_only else {}
    uvlib.update_manifest(rdir, a.stage, {
        "inputs": {"table": uvlib.rel(a.table), "table_sha256": uvlib.sha256(a.table),
                   "marker_index_sha256": uvlib.sha256(index_file), "reads": read_files},
        "output": out_name, "output_sha256": uvlib.sha256(rdir / out_name), "marker_only": a.marker_only,
        "counts": counts, "verdict": "pass"}, a.allow_dirty)
    print(f"{len(results)} rows: {counts} -> {uvlib.rel(rdir / out_name)}")
    for r in results:
        if r["verdict"] not in ("AGREE", "MARKER_AGREES"):
            print(f"  {r['verdict']:<15} {r['symbol']:<11} {r['source']} p{r['page']}: {'; '.join(r['reasons'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
