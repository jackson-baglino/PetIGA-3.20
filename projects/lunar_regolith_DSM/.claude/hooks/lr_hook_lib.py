"""Shared helpers for the literature-review hooks (H1 to H6). Standard library only; never prints paper text."""
import json
import os
import re
import sys
from pathlib import Path

PROJECT = Path(os.environ.get("CLAUDE_PROJECT_DIR") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(PROJECT / "scripts" / "literature_review"))


def stdin_json():
    try:
        return json.load(sys.stdin)
    except ValueError:
        return {}


def run_b1():
    import run_b1 as rb
    return rb


def frozen_patterns():
    return run_b1().FROZEN


def is_frozen(path):
    p = Path(path)
    try:
        rel = p.resolve().relative_to(PROJECT.resolve())
    except ValueError:
        return False
    return any(rel.match(pat) or str(rel) == pat for pat in frozen_patterns())


def stations():
    rb = run_b1()
    names = set()
    for s in rb.stages():
        for step in s.get("steps", []):
            names.add(step[0])
    return names - {"open_gate.py"}


def calls_station(command):
    names = stations()
    return sorted(n for n in names if re.search(r"(^|[\s/])" + re.escape(n) + r"(\s|$)", command))
