#!/usr/bin/env python3
"""Where the B1 line stands, computed from the stage graph and the run folder, never written by hand.

    lr_status.py --brief            one line: phase, freeze, run, next stage, gates and review items waiting
    lr_status.py --write            refresh the generated table in B1_CONTRACT.md ("Line status")
    lr_status.py --check            exit 1 if that table is stale (used by the Stop hook)
    [--run <dir>]                   default: $LR_RUN, else the newest run under runs/literature-review/b1/

Reads only manifest.json, runs.db, gate files and the done files' names and top-level keys; prints no paper text.
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_b1  # noqa: E402

BEGIN = "<!-- BEGIN GENERATED: lr_status.py -->"
END = "<!-- END GENERATED -->"


def table(run):
    if run is None:
        return "No B1 run yet (`run_b1.py init`)."
    rows = run_b1.status_rows(run)
    out = [f"Run `{Path(run).name}`. {run_b1.brief(run)}", "",
           "| Stage | Station | Done | Gate | Waiting |", "|---|---|---|---|---|"]
    for r in rows:
        out.append(f"| {r['stage']} | {'built' if r['built'] else 'not built'} | {r['done']}/{r['total']} | "
                   f"{r['gate']} | {r['waiting'] or ''} |")
    return "\n".join(out)


def block(run):
    return f"{BEGIN}\n{table(run)}\n{END}"


def current_text():
    t = run_b1.CONTRACT.read_text()
    if BEGIN not in t:
        return None, t
    i, j = t.index(BEGIN), t.index(END, t.index(BEGIN)) + len(END)
    return t[i:j], t


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--brief", action="store_true")
    g.add_argument("--write", action="store_true")
    g.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    run = a.run or os.environ.get("LR_RUN") or run_b1.current_run()
    if a.brief:
        print(run_b1.brief(run))
        return 0
    new = block(run)
    old, text = current_text()
    if a.check:
        stale = old is None or old != new
        if stale:
            print("B1_CONTRACT.md line status is stale: run lr_status.py --write")
        return 1 if stale else 0
    if old is None:
        anchor = "## Status\n"
        if anchor not in text:
            raise SystemExit("B1_CONTRACT.md has no '## Status' section")
        text = text.replace(anchor, anchor + "\n### Line status (generated)\n\n" + new + "\n", 1)
    else:
        text = text.replace(old, new)
    run_b1.CONTRACT.write_text(text)
    print("B1_CONTRACT.md line status written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
