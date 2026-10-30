#!/usr/bin/env python3
"""Y3: two fresh skeptical calls on a draft memo, support and omissions (WORKFLOW.md §3).

    skeptic.py --run <run folder> [--question Q1] [--version 1]

Call 1 sees the memo and the claims it cites; call 2 sees the memo, every checked claim on
the question's papers (all questions) and the checklists. Results go next to the memo;
counts are printed.
"""
import argparse
import collections
import json
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call as ac  # noqa: E402
import read_check as rc  # noqa: E402
import synthesize as sy  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--question", default="Q1")
    ap.add_argument("--version", type=int, default=1)
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    mdir = run / "untrusted" / "memos"
    memo = (mdir / f"{a.question}.v{a.version}.md").read_text()
    used = json.loads((mdir / f"{a.question}.v{a.version}.claims.json").read_text())
    cited = set(re.findall(r"\[([a-z0-9_]+-c\d{3})", memo))
    allc = []
    for q in ("Q1", "Q2", "Q3"):
        allc += [c for c in sy.checked_claims(run, q)[0] if c["id"] not in {x["id"] for x in allc}]
    _, lists, _ = sy.checked_claims(run, a.question)
    q = next(l for l in rc.questions().splitlines() if l.startswith(a.question))
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        out["support"], m1 = ac.call(tmp, ac.prompt("skeptic_support", memo=memo,
                                                    claims=json.dumps([c for c in used if c["id"] in cited])),
                                     "skeptic_support", "beth-skeptic")
    with tempfile.TemporaryDirectory() as tmp:
        out["omissions"], m2 = ac.call(tmp, ac.prompt("skeptic_omissions", question=q, memo=memo,
                                                      claims=json.dumps(allc), checklists=json.dumps(lists)),
                                       "skeptic_omissions", "beth-skeptic")
    if out["support"] is None or out["omissions"] is None:
        raise SystemExit(f"Y3 failed: {m1} {m2}")
    (mdir / f"{a.question}.v{a.version}.y3.json").write_text(json.dumps(out, indent=1))
    s = {"statements": dict(collections.Counter(x["verdict"] for x in out["support"]["statements"])),
         "omissions": dict(collections.Counter((x["kind"], x["importance"]) for x in out["omissions"]["omissions"]))}
    s["omissions"] = {f"{k}/{i}": n for (k, i), n in s["omissions"].items()}
    (run / f"summary_y3_{a.question}.v{a.version}.json").write_text(json.dumps(s, indent=1))
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
