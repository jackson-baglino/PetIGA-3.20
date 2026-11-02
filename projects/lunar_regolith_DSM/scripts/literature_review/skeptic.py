#!/usr/bin/env python3
"""Y3: two fresh skeptical calls on a draft memo, support and omissions (WORKFLOW.md §3).

    skeptic.py --run <run folder> [--question Q1] [--version K] [--source library|run]

Call 1 sees the memo and the claims it cites; call 2 sees the memo, every checked claim on
the papers (all questions) and the checklists. In the line (--source library) the claims are
the published ones (synthesize.library_claims); the memo is the latest version unless
--version is given. Results go next to the memo; counts go to summary_y3_<question>.json
(and .v<K>.json); statements found overstated or unsupported, and high omissions, are added
to the review queue for the user at U4. Not built: the one Y1 revision on an unsupported
statement (WORKFLOW.md §3, Y3); until then those statements go to the user.
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
    ap.add_argument("--version", type=int)
    ap.add_argument("--source", choices=["library", "run"], default="library")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    mdir = run / "untrusted" / "memos"
    if a.version is None:
        import check_citations
        a.version = check_citations.latest_version(mdir, a.question)
    memo = (mdir / f"{a.question}.v{a.version}.md").read_text()
    used = json.loads((mdir / f"{a.question}.v{a.version}.claims.json").read_text())
    cited = set(re.findall(r"\[([a-z0-9_]+-c\d{3})", memo))
    if a.source == "library":
        allc, lists, _ = sy.library_claims(run)
        q = rc.questions()
    else:
        allc = []
        for qq in ("Q1", "Q2", "Q3"):
            allc += [c for c in sy.checked_claims(run, qq)[0] if c["id"] not in {x["id"] for x in allc}]
        _, lists, _ = sy.checked_claims(run, a.question)
        q = next(line for line in rc.questions().splitlines() if line.startswith(a.question))
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        out["support"], m1 = ac.call(tmp, ac.prompt("skeptic_support", memo=memo,
                                                    claims=json.dumps([c for c in used if c["id"] in cited])),
                                     "skeptic_support", "beth-skeptic",
                                     ledger=run / "calls.jsonl", step="Y3.support", question=a.question,
                                     version=a.version)
    with tempfile.TemporaryDirectory() as tmp:
        out["omissions"], m2 = ac.call(tmp, ac.prompt("skeptic_omissions", question=q, memo=memo,
                                                      n_papers=str(len({c["paper"] for c in allc})),
                                                      claims=json.dumps(allc), checklists=json.dumps(lists)),
                                       "skeptic_omissions", "beth-skeptic",
                                       ledger=run / "calls.jsonl", step="Y3.omissions", question=a.question,
                                       version=a.version)
    if out["support"] is None or out["omissions"] is None:
        raise SystemExit(f"Y3 failed: {m1} {m2}")
    (mdir / f"{a.question}.v{a.version}.y3.json").write_text(json.dumps(out, indent=1))
    s = {"statements": dict(collections.Counter(x["verdict"] for x in out["support"]["statements"])),
         "omissions": dict(collections.Counter((x["kind"], x["importance"]) for x in out["omissions"]["omissions"]))}
    s["omissions"] = {f"{k}/{i}": n for (k, i), n in s["omissions"].items()}
    s["usage"] = {"support": m1["usage"], "omissions": m2["usage"]}
    s["version"] = a.version
    for name in (f"summary_y3_{a.question}.v{a.version}.json", f"summary_y3_{a.question}.json"):
        (run / name).write_text(json.dumps(s, indent=1))
    import revise
    revise.add_review(run, [{"stage": "Y3", "item": f"{a.question}.v{a.version}", "claim_id": f"statement-{i}",
                             "kind": x["verdict"]} for i, x in enumerate(out["support"]["statements"], 1)
                            if x["verdict"] != "supported"] +
                      [{"stage": "Y3", "item": f"{a.question}.v{a.version}", "claim_id": f"omission-{i}",
                        "kind": x["kind"]} for i, x in enumerate(out["omissions"]["omissions"], 1)
                       if x["importance"] == "high"])
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
