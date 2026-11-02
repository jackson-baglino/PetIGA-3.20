#!/usr/bin/env python3
"""C0 baseline: the calibration plants and controls judged by scripts alone, no model (WORKFLOW.md §5).

    baseline.py --run <run folder> --set tuning|heldout [--papers ...] [--content-run <run with P1 content>]
                [--prepare]

The same items check_plants.py gives the checker (check_plants.items_for), on the same cited PDF page. A claim is
flagged when any of these deterministic checks fails on that page:
- numbers: every number in the claim's text and parameter field (decimal, e-notation, x10^n, Unicode minus and
  superscripts normalized) appears on the page;
- units: a number the claim writes with a unit appears on the page with the same unit (normalized);
- equations: each displayed equation in the claim ($...$ with =, <, > or \\le, \\ge) equals, after tres.canon, an
  equation P1 found on that page; checked only when P1 content for the paper exists (else "not checkable").
The page is P1's page text when it exists, else the PDF text layer (mutool), the text E3 uses to place quotes.
Prose without numbers or equations is not flagged: scripts cannot judge it. That is the point of the comparison.

Output: <run>/baseline_<set>.jsonl in calibrate.py's results format (a plant {"plant_id", "caught_by": "script" |
null}, a control {"control_id", "flagged"}) plus, per item, which checks ran; prints counts only. Score it with
`calibrate.py score --baseline`, which reports and never touches the held-out ledger.
"""
import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate as cal  # noqa: E402
import check_calibration as cc  # noqa: E402
import check_plants as cp  # noqa: E402
import tres  # noqa: E402

SUPER = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
MINUS = str.maketrans({"−": "-", "–": "-", "—": "-"})
# mantissa, then an optional power of ten: e-3, ×10^{-3}, \times 10^{-3}, x 10−3, ·10⁻³
NUM = re.compile(r"(?<![\w.])(-?\d+(?:\.\d+)?)(?:\s*(?:[eE]\s*([-+]?\d+)|\s*(?:×|\\times|\\cdot|·|x)\s*10\s*\^?\s*\{?\s*"
                 r"([-+]?\d+)\s*\}?))?(?![\d])")
UNITS = {"◦c": "degc", "◦": "deg", "µm": "um", "μm": "um", "microns": "um", "micron": "um", "um": "um", "nm": "nm", "mm": "mm", "cm": "cm",
         "m": "m", "km": "km", "°c": "degc", "∘c": "degc", "^\\circ c": "degc", "\\circ c": "degc", "k": "k",
         "s": "s", "ms": "ms", "pa": "pa", "kpa": "kpa", "mpa": "mpa", "j": "j", "kj": "kj", "w": "w", "kg": "kg",
         "g": "g", "mol": "mol", "%": "%", "degrees": "deg", "°": "deg"}
UNIT_RE = re.compile(r"\s*(\^\\circ\s*c|\\circ\s*c|◦\s?c|°\s?c|∘c|µm|μm|microns?|um|nm|mm|cm|km|kpa|mpa|pa|kj|kg|mol|ms|"
                     r"degrees|°|%|m|k|s|j|w|g)(?![a-z])", re.I)
EQ = re.compile(r"\$([^$]*(?:=|<|>|\\le|\\ge)[^$]*)\$")


def norm_text(t):
    return (t or "").translate(MINUS).translate(SUPER).replace(" ", " ").replace("\xa0", " ")


def numbers(text):
    """[(value, unit or None)] for every number in the text; a value is (mantissa, exponent) as a float."""
    t = norm_text(text)
    out = []
    for m in NUM.finditer(t):
        if m.start() > 0 and t[m.start() - 1] == "^":         # a unit's exponent (m^-1, K^-1), not a number
            continue
        exp = m.group(2) or m.group(3)
        try:
            v = float(m.group(1)) * (10.0 ** int(exp) if exp else 1.0)
        except (ValueError, OverflowError):
            continue
        u = UNIT_RE.match(t, m.end())
        out.append((v, UNITS.get(re.sub(r"\s+", "", u.group(1).lower()).replace("\\circc", "\\circ c")) if u else None))
    return out


def same(a, b):
    return a == b or (a != 0 and abs(a - b) <= 1e-9 * abs(a))


def page_text(pdf, page, content):
    """P1's page text when available, else the PDF text layer of that page."""
    if content:
        for name in (f"p{page:03d}.md", f"p{page}.md", f"page_{page:03d}.md"):
            f = content / "pages" / name
            if f.exists():
                return f.read_text()
    return subprocess.run(["mutool", "draw", "-q", "-F", "txt", "-o", "-", str(pdf), str(page)],
                          capture_output=True, text=True).stdout


def page_equations(content, page):
    if not content or not (content / "equations.jsonl").exists():
        return None
    return [tres.canon(e["latex"]) for e in map(json.loads, (content / "equations.jsonl").read_text().splitlines())
            if e.get("page_pdf") == page]


def plain(text, drop_equations=True):
    """Text for the number check, inline math such as $-6^\\circ$C turned into -6°C. For a claim, displayed
    equations are dropped (the equation check covers them); for a page (P1's Markdown) they are kept."""
    if drop_equations:
        text = EQ.sub(" ", text)

    def inline(m):
        t = m.group(1).replace("^\\circ", "°").replace("\\circ", "°").replace("\\times", "×")
        t = t.replace("\\mu", "µ").replace("\\,", " ")
        return re.sub(r"\\[a-zA-Z]+", "", t).replace("{", "").replace("}", "")
    return re.sub(r"\$([^$]*)\$", inline, text)


def judge(item, text, equations):
    """{flag, checks: {numbers, units, equations}} for one claim on its page text."""
    claim = plain(item["text"]) + " " + " ".join(str(v) for v in (item.get("parameter") or {}).values())
    on_page = numbers(plain(text, drop_equations=False))
    checks = {}
    nums = numbers(claim)
    if nums:
        checks["numbers"] = all(any(same(v, w) for w, _ in on_page) for v, _ in nums)
        with_unit = [(v, u) for v, u in nums if u]
        if with_unit:
            checks["units"] = all(any(same(v, w) and u == uu for w, uu in on_page) for v, u in with_unit)
    eqs = [tres.canon(m.group(1)) for m in EQ.finditer(item["text"])]
    if eqs:
        checks["equations"] = "not checkable" if equations is None else all(e in equations for e in eqs)
    flag = any(v is False for v in checks.values())
    return {"flag": flag, "checks": checks}


def content_dir(run, paper, pdf):
    if run is None:
        return None
    h = hashlib.sha256(Path(pdf).read_bytes()).hexdigest()[:16]
    d = Path(run) / "untrusted" / "content" / paper / h
    return d if d.exists() else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", choices=["tuning", "heldout"], required=True)
    ap.add_argument("--papers", nargs="*")
    ap.add_argument("--content-run", help="a run folder holding P1 content for these papers (default: --run)")
    ap.add_argument("--prepare", action="store_true",
                    help="run P1 (prepare.py, in this process and sandbox) for any paper without P1 content")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    split = yaml.safe_load(open(cal.SPLIT))
    ids = set(split["tuning"] if a.set == "tuning" else split["heldout"])
    plants = cal.load_plants()
    papers = a.papers or sorted({p["paper"] for p in plants if p["id"] in ids})
    results, ran = [], {"numbers": 0, "units": 0, "equations": 0, "equations_not_checkable": 0, "none": 0}
    for paper in papers:
        pdf = cp.pdf_of(paper)
        offset = cc.PageMap(str(pdf)).offset
        content = content_dir(a.content_run or run, paper, pdf)
        if content is None and a.prepare:
            import prepare
            if prepare.main(["--run", str(run), "--paper-id", paper, "--pdf", str(pdf)]):
                raise SystemExit(f"P1 failed for {paper}")
            content = content_dir(run, paper, pdf)
        for x in cp.items_for(paper, ids, plants, offset):
            j = judge(x, page_text(pdf, x["page_pdf"], content), page_equations(content, x["page_pdf"]))
            for k, v in j["checks"].items():
                ran[k] += 1
                if v == "not checkable":
                    ran["equations_not_checkable"] += 1
            ran["none"] += not j["checks"]
            if x["role"] == "plant":
                results.append({"plant_id": x["id"], "caught_by": "script" if j["flag"] else None,
                                "checks": j["checks"]})
            else:
                results.append({"control_id": x["id"], "flagged": j["flag"], "checks": j["checks"]})
    with open(run / f"baseline_{a.set}.jsonl", "w") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    pl = [r for r in results if "plant_id" in r]
    ct = [r for r in results if "control_id" in r]
    print(json.dumps({"set": a.set, "plants": len(pl), "caught": sum(bool(r["caught_by"]) for r in pl),
                      "controls": len(ct), "flagged": sum(r["flagged"] for r in ct), "checks_ran": ran}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
