#!/usr/bin/env python3
"""S5 of units-verification: write each row's final status into memory/semantic/units.md.

    python scripts/units_verification/apply_verdicts.py --run-id <id>

Deterministic: it rebuilds the table from the version the run started with (git, at the S0 commit; its sha256
must equal the one S0 recorded), then applies, row by row,
  - the user's U1 decision (inspection.json user_decision): approve, drop, or edit (fields, possibly a new id);
  - for a row not in the inspection set: it must be S3 AGREE and S4m "matches";
  - for an edited row: the re-check (recheck/compare.json, recheck/verdicts.json; recheck.py). Until the re-check
    has run, or if it disputes the row, the row is unresolved.
The `quote (text layer)` column is replaced by `value (printed)`, the value the blind read transcribed from the
page image (corroborated by Marker for agreed rows). The UNRELIABLE header is removed only if no literature row
is unresolved. Writes units.md, applied.json (what happened to each row) and the S5 entry of the manifest.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import uvlib

COLUMNS = ["symbol", "meaning", "unit", "value", "source", "page", "verified"]
HEADERS = ["symbol", "meaning", "unit", "value (printed)", "source", "page", "verified"]
INFERRED = " — meaning inferred from the subscript"
NOTES = uvlib.WORKFLOW_DIR / "table_notes.md"   # user-approved notes printed above the table


def original_table(rdir):
    s0 = uvlib.load_manifest(rdir)["stages"]["S0"]
    path = Path(s0["inputs"]["table"])
    text = subprocess.run(["git", "show", f"{s0['git']['commit']}:./{path}"], cwd=uvlib.ROOT, capture_output=True,
                          text=True, check=True).stdout
    if uvlib.sha256_bytes(text.encode()) != s0["inputs"]["table_sha256"]:
        raise uvlib.StageError("the table at the S0 commit does not match the sha256 S0 recorded")
    tmp = rdir / ".table_at_S0.md"
    tmp.write_text(text)
    try:
        return uvlib.read_table(tmp)
    finally:
        tmp.unlink()


def load(path, default=None):
    return json.loads(Path(path).read_text()) if Path(path).exists() else default


def decide(row, comp, s4m, decision, recheck, run_id, today):
    """(new row or None if dropped, status, record). status: verified | unresolved | dropped | unchanged."""
    sym = row["symbol"]
    if row["source"] not in uvlib.LITERATURE_SOURCES:
        return dict(row, value=""), "unchanged", {"basis": "out of scope (not a literature row)"}
    c = comp.get(sym, {})
    value = (c.get("read") or {}).get("value_as_printed") or ""
    new = dict(row, value=value)
    if s4m.get(sym, {}).get("meaning_verdict") == "matches":
        new["meaning"] = new["meaning"].replace(INFERRED, "")
    stamp = f"verified {today}, run {run_id}"
    if decision is None:
        if c.get("verdict") == "AGREE" and s4m.get(sym, {}).get("meaning_verdict") == "matches":
            new["verified"] = (f"{stamp}: unit, symbol and page agree in table, blind read and Marker; "
                               f"meaning checked against the page")
            return new, "verified", {"basis": "S3 AGREE and S4m matches"}
        new["verified"] = "unresolved: disputed and not decided by the user"
        return new, "unresolved", {"basis": "disputed without a user decision"}
    act, note = decision["action"], decision.get("note", "")
    if act == "drop":
        return None, "dropped", {"basis": f"user (U1): {note}"}
    if act == "approve":
        new["verified"] = f"{stamp}: approved by the user (U1): {note}"
        return new, "verified", {"basis": "user approved"}
    if act == "edit":
        new.update(decision["edit"])
        new["meaning"] = new["meaning"].replace(INFERRED, "")
        m = re.search(r"\(printed as ([^)]+)\)", new["meaning"])
        new["printed"] = m.group(1).strip() if m else None
        rc, rv = (recheck or {}).get("compare", {}), (recheck or {}).get("S4m", {})
        r_c, r_m, r_s4 = rc.get(new["symbol"]), rv.get(new["symbol"], {}), (recheck or {}).get("S4", {})
        if r_c is None:
            new["value"] = ""
            new["verified"] = f"unresolved: corrected by the user (U1), re-check pending: {note}"
            return new, "unresolved", {"basis": "edited; re-check pending", "edit": decision["edit"]}
        new["value"] = (r_c.get("read") or {}).get("value_as_printed") or ""
        meaning_ok = r_m.get("meaning_verdict") == "matches"
        unit_ok = r_c["verdict"] == "AGREE" or r_s4.get(new["symbol"], {}).get("verdict") == "table_right"
        if unit_ok and meaning_ok:
            how = "re-check agrees" if r_c["verdict"] == "AGREE" else "re-check: checker finds the table right"
            new["verified"] = f"{stamp}: corrected by the user (U1), then {how}: {note}"
            return new, "verified", {"basis": how, "edit": decision["edit"]}
        why = "; ".join(r_c["reasons"] + ([f"S4m {r_m.get('meaning_verdict')}"] if not meaning_ok else []))
        new["verified"] = f"unresolved: corrected by the user (U1), re-check disputes it: {why}"
        return new, "unresolved", {"basis": "re-check disputes", "edit": decision["edit"], "reasons": why}
    raise uvlib.StageError(f"{sym}: unknown decision {act!r}")


def header(run_id, today, n_verified, unresolved, dropped):
    lines = ["# Units of symbols (semantic memory)", ""]
    if unresolved:
        lines += [f"> **UNRELIABLE — {len(unresolved)} rows unresolved** ({', '.join(unresolved)}). Do not cite "
                  "a unit from an unresolved row.", ""]
    lines += [f"Checked by the units-verification workflow, run `{run_id}` ({today}; record in",
              f"`runs/units-verification/{run_id}/`). Each literature row's unit, printed symbol,",
              "page and value were read blind from the page image (Sonnet), compared by script with",
              "the table and with an OCR reading (Marker), disputed rows were judged against the",
              "image (Opus), every meaning was reviewed against the page (Opus), and the user",
              "inspected every disputed row and a random 10 % of the rest.", "",
              f"- {n_verified} literature rows verified; {len(unresolved)} unresolved.",
              f"- Removed at the user's decision: {', '.join(dropped) or 'none'}; earlier, `j_T` and `j_v`.",
              "- `value (printed)` is the value as printed on the cited page, transcribed from the",
              "  image; empty where the page prints none (model unknowns such as `rho_v`).",
              "- `M0` and `alpha` are the solver's own parameters, checked against the code, not",
              "  part of this run.", "",
              "Notes (from `workflows/units-verification/table_notes.md`):", ""]
    lines += NOTES.read_text().rstrip("\n").splitlines() + [""]
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--out", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    today = uvlib.now()[:10]
    print(f"S5 apply_verdicts: run={uvlib.rel(rdir)} out={uvlib.rel(a.out)}")
    try:
        uvlib.require_committed(a.allow_dirty)
        rows = original_table(rdir)
        comp = {c["symbol"]: c for c in load(rdir / "compare.json")["rows"]}
        s4m = load(rdir / "verdicts.json", {}).get("S4m", {})
        insp = load(rdir / "inspection.json")
        decisions = {r["symbol"]: r["user_decision"] for r in insp["rows"]}
        undecided = [s for s, d in decisions.items() if d is None]
        if undecided:
            raise uvlib.StageError(f"no user decision yet for {undecided}")
        recheck = None
        if (rdir / "recheck" / "compare.json").exists():
            v = load(rdir / "recheck" / "verdicts.json", {})
            recheck = {"compare": {c["symbol"]: c for c in load(rdir / "recheck" / "compare.json")["rows"]},
                       "S4": v.get("S4", {}), "S4m": v.get("S4m", {})}
        out_rows, record = [], []
        for r in rows:
            new, status, rec = decide(r, comp, s4m, decisions.get(r["symbol"]), recheck, a.run_id, today)
            record.append(dict(rec, symbol=r["symbol"], status=status, new_symbol=new["symbol"] if new else None))
            if new is not None:
                out_rows.append(new)
    except (uvlib.StageError, subprocess.CalledProcessError) as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    lit = [x for x in record if x["status"] not in ("unchanged", "dropped")]
    unresolved = [x["new_symbol"] for x in record if x["status"] == "unresolved"]
    dropped = [x["symbol"] for x in record if x["status"] == "dropped"]
    n_ver = sum(x["status"] == "verified" for x in record)
    text = "\n".join(header(a.run_id, today, n_ver, unresolved, dropped)
                     + uvlib.table_lines(out_rows, COLUMNS, HEADERS)) + "\n"
    a.out.write_text(text)
    uvlib.read_table(a.out)  # the written table must read back
    uvlib.write_json(rdir / "applied.json", {"rows": record, "verified": n_ver, "unresolved": unresolved,
                                             "dropped": dropped})
    uvlib.update_manifest(rdir, "S5", {"out": uvlib.rel(a.out), "out_sha256": uvlib.sha256(a.out),
                                       "verified": n_ver, "unresolved": unresolved, "dropped": dropped,
                                       "recheck_used": recheck is not None,
                                       "verdict": "complete" if not unresolved else "unresolved rows"}, a.allow_dirty)
    print(f"{len(lit)} literature rows: {n_ver} verified, {len(unresolved)} unresolved {unresolved}; "
          f"dropped {dropped} -> {uvlib.rel(a.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
