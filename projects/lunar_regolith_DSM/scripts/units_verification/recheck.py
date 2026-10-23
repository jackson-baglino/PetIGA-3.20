#!/usr/bin/env python3
"""S5 re-check of units-verification: rows the user corrected go through S2 -> S3 -> S4 -> S4m again.

    python scripts/units_verification/apply_verdicts.py --run-id <id>   # writes the corrected rows, pending
    python scripts/units_verification/recheck.py --run-id <id>          # re-checks them
    python scripts/units_verification/apply_verdicts.py --run-id <id>   # writes their final status

Takes the rows of units.md whose status is "re-check pending" (or --symbols), and runs the same pipeline as
the real run on them alone in <run>/recheck/, re-using the run's page images and Marker reading (pages a row
now cites are rendered first). Writes recheck/compare.json, recheck/verdicts.json and the manifest entry.
"""
import argparse
import os
import sys
from pathlib import Path

import check_rows
import compare_units
import marker_extract
import read_pages
import render_pages
import uvlib


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--symbols", nargs="*")
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    cdir = rdir / "recheck"
    print(f"S5-recheck: run={uvlib.rel(rdir)} table={uvlib.rel(a.table)}")
    try:
        uvlib.require_committed(a.allow_dirty)
        stages = uvlib.load_manifest(rdir).get("stages", {})
        n_pages = {s: v["pages"] for s, v in stages["S0"]["sources"].items()}
        rows = [r for r in uvlib.literature_rows(uvlib.read_table(a.table))
                if (r["symbol"] in a.symbols if a.symbols else "re-check pending" in r["verified"])]
        if not rows:
            raise uvlib.StageError("no rows to re-check")
        print("  rows: " + ", ".join(f"{r['symbol']} ({r['source']} p{r['page']})" for r in rows))
        cdir.mkdir(exist_ok=True)
        for name in ("pages", "marker"):
            if not (cdir / name).exists():
                os.symlink(Path("..") / name, cdir / name)
        index = marker_extract.build_index(rows, cdir / "marker", n_pages, marker_extract.related_sources())
        uvlib.write_json(cdir / "marker_index.json", index)
        missing = sorted({f"{r['source']}:{r['page']}" for r in rows
                          if not (rdir / "pages" / f"{r['source']}_p{r['page']}.png").exists()})
        if missing:
            render_pages.main(["--run-id", a.run_id, *sum((["--add", m] for m in missing), [])]
                              + (["--allow-dirty"] if a.allow_dirty else []))

        def s3():
            reads = compare_units.load_reads(cdir)
            res = [compare_units.compare_row(r, reads, cdir / "marker", index) for r in rows]
            uvlib.write_json(cdir / "compare.json", {"table": uvlib.rel(a.table), "rows": res})
            return {"rows": res}
        _, failed = read_pages.run_groups(cdir, "S5-S2", read_pages.cited_groups(rows))
        comp = s3()
        regroups = read_pages.reread_groups(cdir, rows)
        need = sorted({f"{s}:{p}" for s, p in regroups if not (rdir / "pages" / f"{s}_p{p}.png").exists()})
        if need:
            render_pages.main(["--run-id", a.run_id, *sum((["--add", m] for m in need), [])]
                              + (["--allow-dirty"] if a.allow_dirty else []))
        _, failed2 = read_pages.run_groups(cdir, "S5-S2", regroups)
        comp = s3()
        for c in comp["rows"]:
            print(f"  S3 {c['symbol']}: {c['verdict']} {'; '.join(c['reasons'])}")
        s4 = check_rows.run_s4(cdir, rows, comp, index, stage="S5-S4", run_id=a.run_id)
        s4m = check_rows.run_s4m(cdir, rows, comp, index, stage="S5-S4m", run_id=a.run_id)
        uvlib.write_json(cdir / "verdicts.json", {"S4": s4, "S4m": s4m})
    except uvlib.StageError as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    uvlib.update_manifest(rdir, "S5-recheck", {
        "rows": [r["symbol"] for r in rows], "table_sha256": uvlib.sha256(a.table),
        "s3": {c["symbol"]: c["verdict"] for c in comp["rows"]},
        "s4": {k: v.get("verdict") for k, v in s4.items()},
        "s4m": {k: v.get("meaning_verdict") for k, v in s4m.items()},
        "read_failures": failed + failed2, "verdict": "pass"}, a.allow_dirty)
    return 0


if __name__ == "__main__":
    sys.exit(main())
