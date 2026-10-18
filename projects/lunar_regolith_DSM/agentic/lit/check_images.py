#!/usr/bin/env python3
"""Image check: for every page in a record's image queue, a vision-capable worker (Read tool on the rendered page PNG) confirms,
corrects or flags the equations, numbers and numeric claims that came from an unreliable text layer. Closes the gap left by
verify_record (a quote match proves the words, not the signs). Results are stored in the record and its verification file.

  python -m agentic.lit.check_images --record ji_karma_2026_sonnet [--model opus|sonnet] [--workers 4] [--dpi 110]
  python -m agentic.lit.check_images --all            # every record in $LIT_REVIEW_DIR/records that has a non-empty queue
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from agentic.lit import common

SCHEMA = {"type": "object", "required": ["page_legible", "items"], "properties": {
    "page_legible": {"type": "boolean"},
    "items": {"type": "array", "items": {"type": "object", "required": ["id", "status", "image_value", "note"], "properties": {
        "id": {"type": "string"}, "status": {"type": "string", "enum": ["confirmed", "corrected", "unreadable", "not_on_page"]},
        "image_value": {"type": "string"}, "note": {"type": "string"}}}}}}
TEXT = common.STATE / "text"


def png_for(slug: str, page: int, dpi: int) -> Path:
    png = TEXT / f"{slug}_p{page}.png"
    if not png.exists():
        meta = json.loads((TEXT / f"{slug}.pages.json").read_text())
        subprocess.run(["mutool", "draw", "-q", "-r", str(dpi), "-o", str(png), meta["pdf"], str(page)], capture_output=True, check=True)
    return png


def items_by_page(rec: dict, queue_pages: set[int]) -> dict[int, list[dict]]:
    out: dict[int, list[dict]] = {}
    for i, e in enumerate(rec.get("equations", [])):
        if e["page"] in queue_pages or not e.get("text_layer_ok", True):
            out.setdefault(e["page"], []).append({"id": f"E{i}", "kind": "equation", "label": e["label"], "extractor": e.get("latex", "") or "(empty: transcribe it)"})
    for i, n in enumerate(rec.get("numbers", [])):
        if n["page"] in queue_pages or not n.get("text_layer_ok", True):
            out.setdefault(n["page"], []).append({"id": f"N{i}", "kind": "number", "label": n["quantity"], "extractor": f"{n['value']} {n['unit']}"})
    for c in rec.get("claims", []):
        if c["page"] in queue_pages and re.search(r"\d", c["claim"] + c["quote"]):
            out.setdefault(c["page"], []).append({"id": f"C{c['id']}", "kind": "claim", "label": c["claim"][:120], "extractor": f"quote: {c['quote']}"})
    return out


def check_record(stem: str, model: str | None, workers: int, dpi: int) -> str:
    rd = common.REVIEW / "records"
    rec = json.loads((rd / f"{stem}.json").read_text())
    ver = json.loads((rd / f"{stem}.verification.json").read_text())
    slug = re.sub(r"(_(opus|sonnet|haiku|rep\d+))+$", "", stem)
    pages = {q["page"] for q in ver["image_queue"]}
    work = items_by_page(rec, pages)
    if not work:
        return f"{stem}: empty image queue"

    def one(pg):
        png = png_for(slug, pg, dpi)
        lst = "\n".join(f"{it['id']} [{it['kind']}] {it['label']} :: extractor reads: {it['extractor']}" for it in work[pg])
        prompt = f"Page image: {png}\n\nItems to check on this page:\n{lst}"
        out, _ = common.call_claude("image-checker", prompt, model=model, schema=SCHEMA, label=f"img:{stem}:p{pg}", tools="Read", add_dirs=[png.parent])
        return pg, out

    with ThreadPoolExecutor(max_workers=workers) as ex:
        results = dict(ex.map(one, sorted(work)))
    by_id = {it["id"]: dict(it, page=pg) for pg, r in results.items() for it in r["items"]}
    cnt = {"confirmed": 0, "corrected": 0, "unreadable": 0, "not_on_page": 0}
    for it in by_id.values():
        cnt[it["status"]] += 1
    for i, e in enumerate(rec.get("equations", [])):
        r = by_id.get(f"E{i}")
        if r:
            e["image_status"] = r["status"]
            if r["status"] in ("corrected", "confirmed") and r["image_value"]:
                e["latex_image"] = r["image_value"]
    for i, n in enumerate(rec.get("numbers", [])):
        r = by_id.get(f"N{i}")
        if r:
            n["image_status"] = r["status"]
            if r["status"] == "corrected":
                n["value_image"] = r["image_value"]
    rec["image_checked"] = True
    ver["image_check"] = {"model": model or common.parse_agent("image-checker")["model"], "pages": sorted(work), "counts": cnt,
                          "corrections": [dict(id=k, page=v["page"], was=next((x["extractor"] for x in work[v["page"]] if x["id"] == k), ""), now=v["image_value"], note=v["note"])
                                          for k, v in by_id.items() if v["status"] == "corrected"],
                          "unreadable": [dict(id=k, page=v["page"], note=v["note"]) for k, v in by_id.items() if v["status"] == "unreadable"],
                          "illegible_pages": [pg for pg, r in results.items() if not r["page_legible"]]}
    (rd / f"{stem}.json").write_text(json.dumps(rec, indent=1))
    (rd / f"{stem}.verification.json").write_text(json.dumps(ver, indent=1))
    from agentic.lit.run_extract import md_record
    (rd / f"{stem}.md").write_text(md_record(stem, rec, ver))
    return f"{stem}: {len(work)} pages, {sum(cnt.values())} items: {cnt}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record"); ap.add_argument("--all", action="store_true")
    ap.add_argument("--model", default=None); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--dpi", type=int, default=110)
    a = ap.parse_args()
    stems = ([a.record] if a.record else
             sorted(p.name[:-5] for p in (common.REVIEW / "records").glob("*.json") if ".verification" not in p.name and "_rep2" not in p.name and not json.loads(p.read_text()).get("image_checked")))
    for s in stems:
        print(check_record(s, a.model, a.workers, a.dpi), flush=True)


if __name__ == "__main__":
    main()
