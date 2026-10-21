"""Candidate units for the standard, from one paper: extractor -> quote match -> page-image check.  Writes a candidate file only; the curator adds entries afterwards.

  python -m agentic.library.units_extract --id kaempfer_plapp_2009 --pdf <pdf> --text <txt> --out <dir>
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from agentic.lit.common import call_claude  # noqa: E402
from agentic.lit.verify_record import norm  # noqa: E402

ITEMS = {"type": "object", "required": ["items"], "properties": {"items": {"type": "array", "items": {"type": "object",
         "required": ["symbol", "meaning", "unit", "kind", "page", "quote", "stated_unit_as_printed", "basis"],
         "properties": {k: ({"type": "integer"} if k == "page" else {"type": "string"}) for k in ("symbol", "meaning", "unit", "kind", "page", "quote", "stated_unit_as_printed", "basis")}}}}}
IMG = {"type": "object", "required": ["page_legible", "items"], "properties": {"page_legible": {"type": "boolean"}, "items": {"type": "array", "items": {
       "type": "object", "required": ["id", "status", "image_value", "note"], "properties": {k: {"type": "string"} for k in ("id", "status", "image_value", "note")}}}}}


def pages_of(text: str) -> dict[int, str]:
    parts = re.split(r"=== [Pp][Aa][Gg][Ee] (\d+) ===\n", text)
    return {int(parts[i]): parts[i + 1] for i in range(1, len(parts), 2)}


def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("--id", required=True); ap.add_argument("--pdf", required=True); ap.add_argument("--text", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    pages = pages_of(Path(a.text).read_text())
    tagged = "\n".join(f"=== PAGE {n} ===\n{t}" for n, t in pages.items())
    prev = out / f"candidates_{a.id}.json"
    if prev.exists():       # rerun of the image step only
        items = json.loads(prev.read_text())["items"]
    else:
        ex, _ = call_claude("units-extractor", tagged, schema=ITEMS, label=f"units:{a.id}:extract")
        items = ex["items"]
    for i, it in enumerate(items):
        it["id"] = f"u{i}"
        it["quote_found"] = norm(it["quote"]) in norm(pages.get(it["page"], ""))
    by_page: dict[int, list] = {}
    for it in items:
        by_page.setdefault(it["page"], []).append(it)
    imgdir = (out / "pages").resolve(); imgdir.mkdir(exist_ok=True)
    for pg, its in sorted(by_page.items()):
        png = imgdir / f"p{pg}.png"
        if not png.exists():
            subprocess.run(["mutool", "draw", "-r", "130", "-o", str(png), a.pdf, str(pg)], check=True, capture_output=True)
        lst = "\n".join(f"{it['id']}: {it['symbol']} ({it['meaning']}); extractor believes the unit as printed is: {it['stated_unit_as_printed']} (SI: {it['unit']})" for it in its)
        res, _ = call_claude("image-checker", f"Page image: {png}\n\nItems to check on this page:\n{lst}", schema=IMG, label=f"units:{a.id}:img{pg}", tools="Read", add_dirs=[imgdir])
        got = {x["id"]: x for x in res["items"]}
        for it in its:
            it["image"] = got.get(it["id"], {"status": "missing", "image_value": "", "note": ""})
    (out / f"candidates_{a.id}.json").write_text(json.dumps({"paper": a.id, "items": items}, indent=1, ensure_ascii=False))
    n = {k: sum(1 for it in items if it.get("image", {}).get("status") == k) for k in ("confirmed", "corrected", "unreadable", "not_on_page", "missing")}
    print(a.id, len(items), "items;", n, "; quote not found:", sum(not it["quote_found"] for it in items))
    return 0


if __name__ == "__main__":
    sys.exit(main())
