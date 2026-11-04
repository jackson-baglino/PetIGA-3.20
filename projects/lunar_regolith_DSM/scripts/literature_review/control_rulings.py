#!/usr/bin/env python3
"""Read the user's rulings from a control review file (control_review.py) and print them as IDs and categories.

    control_rulings.py --run <calibration run folder> [--set tuning] [--normalize]

The review file is under untrusted/ (it quotes the checker), so a session never opens it. This script prints, per
flagged control, only its ID and the category of the user's ruling: "gold wrong", "checker wrong", "blank" or
"other" (free text the script does not interpret), whether the user left a note under 'Notes:', and the file's
modification time. No text from the file: a note is reported as present or absent, never read out.

A ruling is classified by what it says, not by its first word: "gold wrong", "gold statement is incorrect" and
"checker is right" are "gold wrong"; "checker wrong", "over-flagged" and "the gold statement is right" are
"checker wrong". A ruling that says both, or neither, is "other". --normalize rewrites each classified 'Ruling:' line
to start with its category ("Ruling: gold wrong. <your text>"), keeping the text, after copying the file to
control_review_<set>.md.bak; it still prints only IDs and categories.
"""
import argparse
import datetime
import json
import re
import sys
from pathlib import Path


_CLEAN = re.compile(r"[*_`]")
_GOLD = r"gold(?:\s+statement)?"
_IS = r"(?:\s+(?:is|was))?\s+"
_SAYS = (  # (pattern, category it implies)
    (rf"\b{_GOLD}{_IS}(?:right|correct|fine|ok|accurate|true)\b", "checker wrong"),
    (rf"\bchecker{_IS}(?:wrong|incorrect|mistaken)\b|\bover-?flag", "checker wrong"),
    (rf"\b{_GOLD}{_IS}(?:wrong|incorrect|inaccurate)\b", "gold wrong"),
    (rf"\bchecker{_IS}(?:right|correct)\b", "gold wrong"),
)


def classify(text):
    """'gold wrong', 'checker wrong', 'blank' or 'other' for the text after 'Ruling:'."""
    s = _CLEAN.sub("", text).strip().lower()
    if not s:
        return "blank"
    found = {cat for pat, cat in _SAYS if re.search(pat, s)}
    return found.pop() if len(found) == 1 else "other"


def rulings(md):
    """{control_id: category} from the review file's '## G### (...)' sections and their 'Ruling:' lines."""
    out, cur = {}, None
    for line in md.splitlines():
        m = re.match(r"^## (G\d{3}) ", line)
        if m:
            cur = m.group(1)
            out[cur] = "blank"
            continue
        if cur and line.startswith("Ruling:"):
            out[cur] = classify(line[len("Ruling:"):])
    return out


def normalize(md):
    """The file with each classified 'Ruling:' line starting with its category; the user's text is kept after it."""
    lines, cur = md.split("\n"), None
    for i, line in enumerate(lines):
        m = re.match(r"^## (G\d{3}) ", line)
        if m:
            cur = m.group(1)
            continue
        if cur and line.startswith("Ruling:"):
            rest = line[len("Ruling:"):].strip()
            cat = classify(rest)
            if cat in ("gold wrong", "checker wrong") and not _CLEAN.sub("", rest).lower().startswith(cat):
                lines[i] = f"Ruling: {cat}. {rest}"
    return "\n".join(lines)


# A notes header: 'Notes:', 'Note:', '**Notes:**', '**Notes**:', '### Notes', any case, optionally indented. A colon
# or a heading is required, so a sentence such as "Note that ..." in the checker's reason is not taken for one.
NOTES_HEAD = re.compile(r"^\s*(?:(?:#+\s*)?\**\s*notes?\s*\**\s*:\s*\**|#+\s*notes?\s*$)", re.I)


def noted(md):
    """IDs of the controls whose section has text after a notes header (same line or the lines below it)."""
    out, cur, in_notes = set(), None, False
    for line in md.splitlines():
        m = re.match(r"^## (G\d{3}) ", line)
        if m:
            cur, in_notes = m.group(1), False
            continue
        if not cur:
            continue
        h = NOTES_HEAD.match(line)
        if h:
            in_notes = True
            line = line[h.end():]
        elif line.startswith("Ruling:"):
            in_notes = False
            continue
        if in_notes and line.strip():
            out.add(cur)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", choices=["tuning"], default="tuning")
    ap.add_argument("--normalize", action="store_true", help="start each classified ruling with its category")
    a = ap.parse_args(argv)
    f = Path(a.run) / "untrusted" / f"control_review_{a.set}.md"
    md = f.read_text()
    if a.normalize:
        f.with_name(f.name + ".bak").write_text(md)
        md = normalize(md)
        f.write_text(md)
    r = rulings(md)
    n = noted(md)
    made = datetime.datetime.fromtimestamp(f.stat().st_mtime).isoformat(timespec="seconds")
    counts = {k: sum(1 for v in r.values() if v == k) for k in ("gold wrong", "checker wrong", "other", "blank")}
    counts["with notes"] = len(n)
    print(json.dumps({"file_modified": made, "counts": counts, "rulings": r, "noted": sorted(n)}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
