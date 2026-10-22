#!/usr/bin/env python3
"""U1 of units-verification: build the user's inspection set.

    python scripts/units_verification/inspect_set.py --run-id <id> [--fraction 0.10]

The set is every disputed row (S3 not AGREE, or S4m not "matches", or a failed or inconsistent check) plus a
random sample of the agreed rows: ceil(fraction x AGREE rows), drawn with a seed derived from the run ID so the
draw is repeatable and recorded. Writes
  inspection.json         rows, why each is in the set, all evidence except page text, and an empty
                          user_decision per row (committed)
  inspection/index.html   the same, one section per row with its page images and Marker lines (local only:
                          page images and page text are not committed)
"""
import argparse
import hashlib
import html
import json
import math
import random
import sys
from pathlib import Path

import check_rows
import compare_units
import uvlib


def select(comp_rows, s4, s4m, fraction, seed):
    """(disputed [(symbol, reasons)], sampled agreed [symbol]) in table order."""
    disputed, agreed = [], []
    for c in comp_rows:
        sym, why = c["symbol"], []
        if c["verdict"] != "AGREE":
            why.append(f"S3 {c['verdict']}")
        v = s4.get(sym)
        if v and (v.get("verdict") in ("cannot_tell", "call_failed") or
                  str(v.get("consistency", "")).startswith("inconsistent")):
            why.append(f"S4 {v.get('verdict')} {v.get('consistency', '')}".strip())
        m = s4m.get(sym, {})
        if m.get("meaning_verdict") != "matches":
            why.append(f"S4m {m.get('meaning_verdict')}")
        (disputed.append((sym, why)) if why else agreed.append(sym))
    k = math.ceil(fraction * len(agreed))
    sample = sorted(random.Random(seed).sample(agreed, k), key=agreed.index) if agreed else []
    return disputed, sample


def page_list(c, v, m):
    pages = [(c["source"], c["page"])]
    for p in (v or {}).get("images", []):
        pages.append((c["source"], p))
    for name in (m or {}).get("pages_shown", []):
        src, _, p = name.rpartition("_p")
        pages.append((src, int(p)))
    out = []
    for x in pages:
        if x not in out:
            out.append(x)
    return out


def render_html(entries, rdir, run_id):
    css = ("body{font-family:-apple-system,Helvetica,sans-serif;margin:24px;max-width:1500px}"
           "section{border-top:2px solid #888;padding:12px 0;margin-top:18px}"
           "table{border-collapse:collapse;margin:6px 0}td,th{border:1px solid #ccc;padding:4px 8px;"
           "vertical-align:top;text-align:left;font-size:13px}th{background:#f3f3f3}"
           ".why{color:#b3261e;font-weight:bold}.imgs img{max-width:740px;border:1px solid #aaa;margin:4px}"
           "pre{white-space:pre-wrap;font-size:12px;background:#f7f7f7;padding:6px}")
    e = lambda x: html.escape("" if x is None else str(x))  # noqa: E731
    out = [f"<!doctype html><meta charset='utf-8'><title>units-verification U1 {e(run_id)}</title>"
           f"<style>{css}</style><h1>Units verification: inspection set (run {e(run_id)})</h1>"
           f"<p>{len(entries)} rows. For each row, look at the page image at the given location and decide: "
           "approve the table row as is, or correct it. Record decisions in <code>inspection.json</code> "
           "(<code>user_decision</code>) or tell Claude.</p>"]
    for i, en in enumerate(entries, 1):
        t, rd, v, m = en["table"], en.get("blind_read") or {}, en.get("S4") or {}, en.get("S4m") or {}
        out.append(f"<section id='{e(en['symbol'])}'><h2>{i}. {e(en['symbol'])} "
                   f"<small>({e(t['source'])} p{e(t['page'])})</small></h2>"
                   f"<p class='why'>{e(en['in_set_because'])}</p><table><tr><th></th><th>printed as</th>"
                   "<th>meaning</th><th>unit</th><th>value</th><th>page / location</th></tr>"
                   f"<tr><th>table</th><td>{e(t['printed'])}</td><td>{e(t['meaning'])}</td><td>{e(t['unit'])}</td>"
                   f"<td>—</td><td>p{e(t['page'])}</td></tr>"
                   f"<tr><th>blind read (Sonnet)</th><td>{e(rd.get('printed_as'))}</td><td>{e(rd.get('meaning'))}"
                   f"</td><td>{e(rd.get('unit_as_printed'))}</td><td>{e(rd.get('value_as_printed'))}</td>"
                   f"<td>p{e(rd.get('page'))}: {e(rd.get('location'))}</td></tr>")
        if v:
            out.append(f"<tr><th>S4 checker (Opus): {e(v.get('verdict'))}</th><td>{e(v.get('printed_as'))}</td>"
                       f"<td>{e(v.get('meaning'))}</td><td>{e(v.get('unit'))}</td><td>{e(v.get('value'))}</td>"
                       f"<td>p{e(v.get('page'))}: {e(v.get('location'))}</td></tr>")
        if m:
            out.append(f"<tr><th>S4m meaning (Opus): {e(m.get('meaning_verdict'))}</th>"
                       f"<td>{e(m.get('printed_as_on_page'))}</td><td>{e(m.get('page_wording'))}</td><td></td><td></td>"
                       f"<td>{e(m.get('location'))}</td></tr>")
        out.append("</table>")
        out.append(f"<p><b>S3:</b> {e(en['S3']['verdict'])}: {e('; '.join(en['S3']['reasons']) or 'all agree')}</p>")
        if v.get("reason"):
            out.append(f"<p><b>S4 reason:</b> {e(v['reason'])} <i>({e(v.get('consistency'))})</i></p>")
        if m.get("reason"):
            out.append(f"<p><b>S4m reason:</b> {e(m['reason'])}</p>")
        if en.get("marker_lines"):
            out.append(f"<p><b>Marker lines (OCR, may be wrong):</b></p><pre>{e(en['marker_lines'])}</pre>")
        out.append("<div class='imgs'>" + "".join(
            f"<figure><figcaption>{e(s)} page {p}</figcaption><img src='../pages/{e(s)}_p{p}.png'></figure>"
            for s, p in en["pages"]) + "</div></section>")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--table", type=Path, default=uvlib.UNITS_MD)
    ap.add_argument("--fraction", type=float, default=0.10)
    ap.add_argument("--allow-dirty", action="store_true", help="development only: run with uncommitted definitions")
    a = ap.parse_args(argv)
    rdir = uvlib.run_dir(a.run_id)
    seed = int(hashlib.sha256(a.run_id.encode()).hexdigest()[:12], 16)
    print(f"U1 inspect_set: run={uvlib.rel(rdir)} fraction={a.fraction} seed={seed}")
    try:
        uvlib.require_committed(a.allow_dirty)
        comp = json.loads((rdir / "compare.json").read_text())
        verdicts = json.loads((rdir / "verdicts.json").read_text())
        s4, s4m = verdicts.get("S4", {}), verdicts.get("S4m", {})
        if "S4m" not in verdicts or not s4m:
            raise uvlib.StageError("no S4m verdicts; run check_rows.py --task S4m first")
        rows = {r["symbol"]: r for r in uvlib.literature_rows(uvlib.read_table(a.table))}
    except (uvlib.StageError, FileNotFoundError) as e:
        print(f"STOP: {e}", file=sys.stderr)
        return 2
    disputed, sample = select(comp["rows"], s4, s4m, a.fraction, seed)
    by_c = {c["symbol"]: c for c in comp["rows"]}
    entries, local = [], []
    for sym, why in [(s, w) for s, w in disputed] + [(s, ["random sample of agreed rows"]) for s in sample]:
        c, v, m, r = by_c[sym], s4.get(sym), s4m.get(sym), rows[sym]
        strip = lambda d: {k: x for k, x in (d or {}).items() if k not in ("call", "pass1", "error")}  # noqa: E731
        en = {"symbol": sym, "in_set_because": "; ".join(why),
              "table": {k: r[k] for k in ("printed", "meaning", "unit", "source", "page")},
              "blind_read": c.get("read"), "S3": {"verdict": c["verdict"], "reasons": c["reasons"]},
              "S4": strip(v) or None, "S4m": strip(m) or None, "pages": page_list(c, v, m),
              "user_decision": None}
        entries.append(en)
        lines = check_rows.marker_excerpt(rdir / "marker", c["source"], (c.get("read") or {}).get("page", c["page"]),
                                          compare_units.printed_forms(r, c.get("read")))
        local.append(dict(en, marker_lines=lines))
    uvlib.write_json(rdir / "inspection.json", {
        "run_id": a.run_id, "fraction": a.fraction, "seed": seed, "disputed": len(disputed),
        "agreed": len(comp["rows"]) - len(disputed), "sampled": len(sample), "rows": entries})
    (rdir / "inspection").mkdir(exist_ok=True)
    (rdir / "inspection" / "index.html").write_text(render_html(local, rdir, a.run_id))
    uvlib.update_manifest(rdir, "U1", {"seed": seed, "fraction": a.fraction, "disputed": len(disputed),
                                       "sampled": sample, "output_sha256": uvlib.sha256(rdir / "inspection.json"),
                                       "verdict": "awaiting user"}, a.allow_dirty)
    print(f"{len(disputed)} disputed + {len(sample)} sampled of {len(comp['rows']) - len(disputed)} agreed "
          f"-> {uvlib.rel(rdir / 'inspection.json')}, {uvlib.rel(rdir / 'inspection' / 'index.html')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
