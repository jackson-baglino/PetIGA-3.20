#!/usr/bin/env python3
"""Y2: check a draft memo's citations against the claims (WORKFLOW.md §3, §5).

    check_citations.py --memo <memo.md> --claims <claims.json>

Per sentence outside the exempt sections (Model options, Library silent, Gaps): at least one
[claim_id] citation, every cited id a known checked claim. Every number in a cited sentence
must be followed by a tag [claim_id quantity] whose claim contains that number in its text,
equation or parameter. Prints counts and the failing sentence numbers, never memo text.
"""
import argparse
import json
import re

EXEMPT = re.compile(r"^#+\s*(\d+\.\s*)?\**\s*(model options|library silent|gaps)", re.I)
CITE = re.compile(r"\[([a-z0-9_]+-c\d{3})(?:@v\d+)?(?:\s+([^\]]*))?\]")
LABEL = re.compile(r"^\s*\*\*[^*]+\*\*\s*$")             # a bold lead-in line, a subheading in effect
INTERVAL = re.compile(r"\[\s*[-−]?\d+(?:\.\d+)?\s*,\s*[-−]?\d+(?:\.\d+)?\s*\]")   # math such as [0,1]
NUMBER = re.compile(r"(?<![\w.\-])[-−]?\d+(?:\.\d+)?(?:\s*[×x]\s*10\^?\{?[-−]?\d+\}?)?(?![\w])")


def sentences(md):
    """(index, section exempt?, text) for each sentence or bullet."""
    out, exempt, i = [], False, 0
    for line in md.splitlines():
        if line.startswith("#"):
            exempt = bool(EXEMPT.match(line))
            continue
        if LABEL.match(line):
            continue
        protected = re.sub(r"\b(Eqs?|Figs?|Sec|Ref|e\.g|i\.e|al|cf|vs|approx)\.", lambda m: m.group(0)[:-1] + "\u2024", line)
        for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(\[])", protected.strip(" -*")):
            s = s.replace("\u2024", ".")
            if len(s.split()) >= 4:
                i += 1
                out.append((i, exempt, s))
    return out


def claim_text(c):
    return " ".join([c.get("text", ""), c.get("corrected_text", "") or "", json.dumps(c.get("equation") or {}),
                     json.dumps(c.get("parameter") or {}), json.dumps(c.get("conditions") or [])])


def check(md, claims):
    by = {c["id"]: c for c in claims}
    res = {"sentences": 0, "exempt": 0, "uncited": [], "unknown_ids": [], "numbers": 0, "untagged_numbers": [],
           "tag_mismatch": []}
    for i, ex, s in sentences(md):
        res["sentences"] += 1
        cites = CITE.findall(s)
        for cid, _ in cites:
            if cid not in by:
                res["unknown_ids"].append(i)
        if ex:
            res["exempt"] += 1
            continue
        if not cites:
            res["uncited"].append(i)
        # Citations and intervals blanked to the same length, so positions in `bare` are positions in `s`.
        bare = INTERVAL.sub(lambda m: " " * len(m.group(0)), CITE.sub(lambda m: " " * len(m.group(0)), s))
        for m in NUMBER.finditer(bare):
            n = m.group(0).strip()
            if re.fullmatch(r"\(?\d+(\.\d+)?\)?", n) and re.search(r"(Eq|Section|Fig|Table|page|\()\.?\s*$", bare[:m.start()]):
                continue                                    # an equation, section or figure number
            res["numbers"] += 1
            tag = CITE.search(s[m.end():m.end() + 80])      # the first tag after this number
            if not tag:
                res["untagged_numbers"].append(i)
            elif tag.group(1) in by and n.replace("−", "-").lstrip("-") not in claim_text(by[tag.group(1)]):
                res["tag_mismatch"].append(i)
    res["pass"] = not (res["uncited"] or res["unknown_ids"] or res["untagged_numbers"] or res["tag_mismatch"])
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--memo", required=True)
    ap.add_argument("--claims", required=True)
    a = ap.parse_args(argv)
    r = check(open(a.memo).read(), json.load(open(a.claims)))
    print(json.dumps({k: (len(v) if isinstance(v, list) else v) for k, v in r.items()}))
    print("failing sentence numbers:", {k: v for k, v in r.items() if isinstance(v, list) and v})


if __name__ == "__main__":
    main()
