#!/usr/bin/env python3
"""Read the user's rulings from a control review file (control_review.py) and print them as IDs and categories.

    control_rulings.py --run <calibration run folder> [--set tuning]

The review file is under untrusted/ (it quotes the checker), so a session never opens it. This script prints, per
flagged control, only its ID and the category of the user's ruling: "gold wrong", "checker wrong", "blank" or
"other" (free text the script does not interpret), plus the file's modification time. No text from the file.
"""
import argparse
import datetime
import json
import re
import sys
from pathlib import Path


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
            r = line[len("Ruling:"):].strip().lower()
            if not r:
                out[cur] = "blank"
            elif re.match(r"^\**\s*gold", r):
                out[cur] = "gold wrong"
            elif re.match(r"^\**\s*checker", r):
                out[cur] = "checker wrong"
            else:
                out[cur] = "other"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--set", choices=["tuning"], default="tuning")
    a = ap.parse_args(argv)
    f = Path(a.run) / "untrusted" / f"control_review_{a.set}.md"
    r = rulings(f.read_text())
    made = datetime.datetime.fromtimestamp(f.stat().st_mtime).isoformat(timespec="seconds")
    counts = {k: sum(1 for v in r.values() if v == k) for k in ("gold wrong", "checker wrong", "other", "blank")}
    print(json.dumps({"file_modified": made, "counts": counts, "rulings": r}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
