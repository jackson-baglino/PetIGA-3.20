#!/usr/bin/env python3
"""Y1 and Y2: draft the memo from the checked claims and check its citations (WORKFLOW.md §3).

    synthesize.py --run <run folder> [--question Q1] [--version 1] [--revise]

With --revise, version k is a revision of version k-1: the synthesizer gets its own draft,
the sentences Y2 failed and the same claims, and fixes only those (WORKFLOW.md §3, Y2
failure -> Y1).

Claims enter if the checker (E4 step 2) called them supported, or overstated with a
corrected text, which replaces the reader's; unsupported and not-on-cited-page claims are
left out. Only claims tagged with the question go in. The memo goes to
<run>/untrusted/memos/<question>.v<k>.md; Y2's result (counts only) is printed and saved.
"""
import argparse
import json
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call as ac  # noqa: E402
import check_citations as cx  # noqa: E402
import read_check as rc  # noqa: E402


def checked_claims(run, question):
    out, lists, counts = [], {}, {}
    for d in sorted((run / "untrusted").iterdir()):
        if not (d / "read.json").exists() or not (d / "verdicts.json").exists():
            continue
        claims = json.loads((d / "read.json").read_text())["claims"]
        verdicts = {v["claim_id"]: v for v in json.loads((d / "verdicts.json").read_text())}
        kept = 0
        for c in claims:
            v = verdicts.get(c["id"])
            if not v or question not in c.get("questions", []):
                continue
            if v["text"] == "supported" or (v["text"] == "overstated" and v.get("corrected_text")):
                kept += 1
                out.append({"id": c["id"], "paper": d.name, "verdict": v["text"],
                            "text": v.get("corrected_text") if v["text"] == "overstated" else c["text"],
                            "page_pdf": c["page_pdf"], "kind": c["kind"], "basis": c["basis"],
                            "equation": c.get("equation"), "parameter": c.get("parameter"),
                            "conditions": c["conditions"] + v.get("dropped_conditions", []),
                            "meaning": v.get("meaning")})
        counts[d.name] = {"claims": len(claims), "for_question": kept}
        if (d / "checklist.json").exists():
            lists[d.name] = json.loads((d / "checklist.json").read_text())["items"]
    return out, lists, counts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--question", default="Q1")
    ap.add_argument("--version", type=int, default=1)
    ap.add_argument("--revise", action="store_true")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    claims, lists, counts = checked_claims(run, a.question)
    q = next(l for l in rc.questions().splitlines() if l.startswith(a.question))
    outdir = run / "untrusted" / "memos"
    outdir.mkdir(parents=True, exist_ok=True)
    if a.revise:
        prev = f"{a.question}.v{a.version - 1}"
        claims = json.loads((outdir / f"{prev}.claims.json").read_text())   # the same claims as the draft
        draft = (outdir / f"{prev}.md").read_text()
        y2 = cx.check(draft, claims)
        text = {i: t for i, _, t in cx.sentences(draft)}
        failures = [f"- {k}: {text[i]}" for k in ("uncited", "unknown_ids", "untagged_numbers", "tag_mismatch")
                    for i in dict.fromkeys(y2[k])]
        if not failures:
            raise SystemExit(f"v{a.version - 1} already passes Y2; nothing to revise")
        text = ac.prompt("synthesizer_revise", question=q, memo=draft, failures="\n".join(failures),
                         claims=json.dumps(claims), checklists=json.dumps(lists))
    else:
        text = ac.prompt("synthesizer", question=q, claims=json.dumps(claims), checklists=json.dumps(lists))
    (outdir / f"{a.question}.v{a.version}.claims.json").write_text(json.dumps(claims, indent=1))
    with tempfile.TemporaryDirectory() as tmp:
        res, meta = ac.call(tmp, text, "memo", "beth-synthesizer", ledger=run / "calls.jsonl",
                            step="Y1.revise" if a.revise else "Y1", question=a.question, version=a.version)
    if res is None:
        raise SystemExit(f"Y1 failed: {meta}")
    memo = outdir / f"{a.question}.v{a.version}.md"
    memo.write_text(res["memo"])
    y2 = cx.check(res["memo"], claims)
    cited = sorted(set(re.findall(r"\[([a-z0-9_]+-c\d{3})", res["memo"])))
    summary = {"question": a.question, "version": a.version, "revision_of": a.version - 1 if a.revise else None,
               "claims_given": len(claims), "per_paper": counts,
               "claims_cited": len(cited), "papers_cited": len({c.rsplit("-c", 1)[0] for c in cited}),
               "y2": {k: (len(v) if isinstance(v, list) else v) for k, v in y2.items()},
               "y2_failing_sentences": {k: v for k, v in y2.items() if isinstance(v, list) and v},
               "seconds": meta["seconds"], "usage": meta["usage"]}
    (run / f"summary_memo_{a.question}.v{a.version}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
