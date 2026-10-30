#!/usr/bin/env python3
"""E1, E3 and E4 for one paper: read it, check the quotes, check the claims (WORKFLOW.md §3).

    read_check.py --run <run folder> --paper-id <id> --pdf <file> [--stages E1 E3 E4]

Pages are given to the agents as whole-page images scaled to 2,000 px on the long edge
(T-res, §9 item 45). Steps, each a fresh headless call:
  E1   morty-reader: the whole paper and the brief's questions -> claims, symbols
  E3   script: each quote found on its cited page in the PDF's text layer, unique, >= 8 words
  E4.0 beth-checker: the whole paper, no claims -> applicability checklist
  E4.1 beth-checker: equation and value claims named only by number or quantity and page
       -> blind transcription, compared with the reader's by script
  E4.2 beth-checker, per group of up to 4 cited pages: the claims, the checklist, E3 flags
       -> verdicts with meaning fields
Everything is written under <run>/untrusted/<paper_id>/; a summary with counts only goes to
<run>/summary_read_<paper_id>.json. Nothing here prints paper text.
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call as ac  # noqa: E402
import tres  # noqa: E402

CONTRACT = ac.PROJECT / "workflows" / "literature-review" / "B1_CONTRACT.md"
EDGE = 2000
GROUP = 4                                   # cited pages per E4.2 call
MAX_IMAGES = 20


def questions():
    """Q1 to Q3 from the B1 contract, as plain text for the prompts."""
    text = CONTRACT.read_text()
    out = []
    for q in ("Q1", "Q2", "Q3"):
        m = re.search(r"\*\*" + q + r"[^*]*\*\*(.*?)(?=\n\s*\n|\n- \*\*Q)", text, flags=re.S)
        out.append(f"{q}: " + " ".join(m.group(1).split()).lstrip("(").strip())
    return "\n".join(out)


def render(pdf, outdir):
    """Pages at 300 dpi, scaled to EDGE px on the long edge; and the text layer per page."""
    img, txt = outdir / "pages", outdir / "text"
    if img.exists():
        return sorted(img.glob("p*.png"))
    img.mkdir(parents=True)
    txt.mkdir(parents=True)
    subprocess.run(["mutool", "draw", "-q", "-r", "300", "-c", "gray", "-o", str(img / "p%03d.png"), str(pdf)],
                   check=True)
    for f in sorted(img.glob("p*.png")):
        im = Image.open(f)
        s = EDGE / max(im.size)
        if s < 1:
            im.resize((round(im.size[0] * s), round(im.size[1] * s)), Image.LANCZOS).save(f)
        p = int(f.stem[1:])
        subprocess.run(["mutool", "draw", "-q", "-F", "txt", "-o", str(txt / f"p{p:03d}.txt"), str(pdf), str(p)],
                       check=True)
    return sorted(img.glob("p*.png"))


def folder_with(pages, tmp):
    for f in pages:
        shutil.copy(f, Path(tmp) / f.name)
    return ", ".join(f.name for f in pages)


# ---- E3 ----

def norm(s):
    s = s.replace("ﬁ", "fi").replace("ﬂ", "fl").replace("ﬀ", "ff").replace("ﬃ", "ffi")
    s = re.sub(r"[‘’“”\"']", "", s)
    s = re.sub(r"[‐-―-]\s*\n\s*", "", s)          # hyphenated line breaks
    s = re.sub(r"[‐-―]", "-", s)
    return re.sub(r"\s+", " ", s).strip().lower()


def check_quotes(claims, textdir):
    pages = {int(f.stem[1:]): norm(f.read_text()) for f in textdir.glob("p*.txt")}
    out = []
    for c in claims:
        q = norm(c["quote"])
        where = [p for p, t in pages.items() if q and q in t]
        out.append({"claim_id": c["id"], "words": len(q.split()), "on_cited_page": c["page_pdf"] in where,
                    "found_on": where, "ok": len(q.split()) >= 8 and where == [c["page_pdf"]]})
    return out


# ---- driver ----

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--paper-id", required=True)
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--stages", nargs="+", default=["E1", "E3", "E4"])
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    out = run / "untrusted" / a.paper_id
    out.mkdir(parents=True, exist_ok=True)
    pages = render(Path(a.pdf), out)
    by_num = {int(f.stem[1:]): f for f in pages}
    qs = questions()
    log = {"paper_id": a.paper_id, "pages": len(pages), "calls": []}
    t0 = time.time()

    def save(name, obj):
        (out / f"{name}.json").write_text(json.dumps(obj, indent=1))

    def load(name):
        f = out / f"{name}.json"
        return json.loads(f.read_text()) if f.exists() else None

    if "E1" in a.stages and load("read") is None:
        with tempfile.TemporaryDirectory() as tmp:
            files = folder_with(pages, tmp)
            res, meta = ac.call(tmp, ac.prompt("reader", files=files, questions=qs), "read", "morty-reader")
        log["calls"].append({"step": "E1", **meta})
        if res is None:
            raise SystemExit(f"E1 failed: {meta}")
        for i, c in enumerate(res["claims"], 1):
            c["id"] = f"{a.paper_id}-c{i:03d}"
        save("read", res)
    read = load("read")
    claims = read["claims"]

    if "E3" in a.stages:
        save("e3", check_quotes(claims, out / "text"))
    e3 = {r["claim_id"]: r for r in (load("e3") or [])}

    if "E4" in a.stages:
        if load("checklist") is None:
            with tempfile.TemporaryDirectory() as tmp:
                files = folder_with(pages, tmp)
                res, meta = ac.call(tmp, ac.prompt("checker_applicability", files=files, questions=qs),
                                    "checklist", "beth-checker")
            log["calls"].append({"step": "E4.0", **meta})
            if res is None:
                raise SystemExit(f"E4.0 failed: {meta}")
            save("checklist", res)
        checklist = load("checklist")

        if load("blind") is None:
            items, need = [], set()
            for c in claims:
                if c.get("equation"):
                    items.append(f"- {c['id']}: the equation numbered {c['equation']['number'] or '(unnumbered)'}"
                                 f" on page {c['page_pdf']} ({c['location']})")
                elif c.get("parameter"):
                    items.append(f"- {c['id']}: the value of {c['parameter']['quantity']} on page {c['page_pdf']}")
                else:
                    continue
                need.add(c["page_pdf"])
            blind = {"items": []}
            if items:
                with tempfile.TemporaryDirectory() as tmp:
                    files = folder_with([by_num[p] for p in sorted(need) if p in by_num][:MAX_IMAGES], tmp)
                    blind, meta = ac.call(tmp, ac.prompt("checker_blind", files=files, items="\n".join(items)),
                                          "blind", "beth-checker")
                log["calls"].append({"step": "E4.1", **meta})
                if blind is None:
                    raise SystemExit(f"E4.1 failed: {meta}")
            save("blind", blind)
        blind = {b["claim_id"]: b for b in load("blind")["items"]}

        if load("verdicts") is None:
            cited = sorted({c["page_pdf"] for c in claims})
            groups = [cited[i:i + GROUP] for i in range(0, len(cited), GROUP)]
            verdicts = []
            for g in groups:
                cs = [c for c in claims if c["page_pdf"] in g]
                cond = sorted({k["page_pdf"] for c in cs for k in c["conditions"]} - set(g))
                show = [by_num[p] for p in (g + cond)[:MAX_IMAGES] if p in by_num]
                shown = [{k: c[k] for k in c if k != "symbols"} for c in cs]
                for c in shown:
                    c["quote_check"] = "ok" if e3.get(c["id"], {}).get("ok") else "failed"
                with tempfile.TemporaryDirectory() as tmp:
                    files = folder_with(show, tmp)
                    res, meta = ac.call(tmp, ac.prompt("checker_verdicts", files=files,
                                                       checklist=json.dumps(checklist["items"]),
                                                       claims=json.dumps(shown)), "verdicts", "beth-checker")
                log["calls"].append({"step": "E4.2", "pages": g, **meta})
                if res is None:
                    raise SystemExit(f"E4.2 failed on pages {g}: {meta}")
                verdicts += res["verdicts"]
            save("verdicts", verdicts)
        verdicts = load("verdicts")

        # Compare the reader's equations and values with the blind transcription.
        comp = []
        for c in claims:
            b = blind.get(c["id"])
            if not b:
                continue
            if c.get("equation"):
                same = (not b.get("not_found")) and tres.canon(c["equation"]["latex"]) == tres.canon(b.get("latex", ""))
            else:
                bp = b.get("parameter") or {}
                same = (not b.get("not_found")) and str(c["parameter"]["value"]).replace(" ", "") == \
                    str(bp.get("value", "")).replace(" ", "")
            comp.append({"claim_id": c["id"], "agree": bool(same), "not_found": bool(b.get("not_found"))})
        save("blind_compare", comp)

    # summary: counts only
    s = {"paper_id": a.paper_id, "pages": len(pages), "claims": len(claims),
         "claims_by_kind": {}, "conditions": sum(len(c["conditions"]) for c in claims),
         "symbols": len(read["symbols"]), "flags": len(read["flags"]), "seconds": round(time.time() - t0, 1),
         "calls": [{k: m.get(k) for k in ("step", "agent", "seconds", "exit", "num_turns")} for m in log["calls"]]}
    for c in claims:
        s["claims_by_kind"][c["kind"]] = s["claims_by_kind"].get(c["kind"], 0) + 1
    if e3:
        s["e3_quotes_ok"] = sum(r["ok"] for r in e3.values())
        s["e3_on_cited_page"] = sum(r["on_cited_page"] for r in e3.values())
    if (out / "checklist.json").exists():
        s["checklist_items"] = len(load("checklist")["items"])
    if (out / "verdicts.json").exists():
        v = load("verdicts")
        s["verdicts"] = {k: sum(1 for x in v if x["text"] == k)
                         for k in ("supported", "overstated", "unsupported", "not_on_cited_page")}
        s["verdicts_missing"] = len({c["id"] for c in claims} - {x["claim_id"] for x in v})
        s["dropped_conditions"] = sum(len(x["dropped_conditions"]) for x in v)
        s["equation_corrected"] = sum(1 for x in v if x["equation"] == "corrected")
        s["parameter_corrected"] = sum(1 for x in v if x["parameter"] == "corrected")
    if (out / "blind_compare.json").exists():
        bc = load("blind_compare")
        s["blind_items"] = len(bc)
        s["blind_agree"] = sum(b["agree"] for b in bc)
    (run / f"summary_read_{a.paper_id}.json").write_text(json.dumps(s, indent=1))
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
