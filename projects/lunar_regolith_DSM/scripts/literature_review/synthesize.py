#!/usr/bin/env python3
"""Y1: draft the memo from the published claims, and check its citations (WORKFLOW.md §3).

    synthesize.py --run <run folder> [--question Q1] [--version 1] [--revise] [--source library|run]

In the line (--source library, the default) the claims are the latest published version of every claim on the
contract's papers (curator.py), for Q1 to Q3, minus any the user rejected in U3 or X1 (u3_decisions.json); the
checklists are the papers' applicability records from this run. The memo answers Q1 and has short Q2 and Q3
sections and a Leads list. With --source run (development, b1-dev2) claims come straight from a run's read.json
and verdicts.json: supported, or overstated with a corrected text, and tagged with --question.

With --revise, version k is a revision of version k-1: the synthesizer gets its own draft, the sentences Y2
failed and the same claims, and fixes only those (WORKFLOW.md §3, Y2 failure -> Y1). The memo goes to
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


QUESTIONS = ("Q1", "Q2", "Q3")


def user_rejected(run):
    """Claim IDs the user rejected (wrong) or left unsure in U3 or X1 (u3_decisions.json, written by u3_review.py).
    u3_review.py: unsure is kept out of the memo too."""
    f = Path(run) / "u3_decisions.json"
    if not f.exists():
        return set()
    return {k for k, d in json.loads(f.read_text()).items() if d.get("verdict") in ("wrong", "unsure")}


def library_claims(run, library=None):
    """The latest published version of every claim on the contract's papers (curator.py), for Q1 to Q3, minus
    those the user rejected; and each paper's applicability checklist from this run's analysis records."""
    import curator
    import run_b1
    library = Path(library or curator.LIBRARY)
    out, lists, counts = [], {}, {}
    rejected = user_rejected(run)
    for paper in run_b1.papers(check_hash=False):
        pdir = library / "papers" / paper["id"]
        latest = {}
        for f in (pdir / "claims").glob("*.v*.json"):
            cid, v = f.name.rsplit(".v", 1)[0], int(f.name.rsplit(".v", 1)[1].split(".")[0])
            if cid not in latest or v > latest[cid][0]:
                latest[cid] = (v, f)
        kept = 0
        for cid, (v, f) in sorted(latest.items()):
            rec = json.loads(f.read_text())
            c, vd = rec["claim"], rec["verdict"]
            if cid in rejected or not set(c.get("questions", [])) & set(QUESTIONS):
                continue
            kept += 1
            out.append({"id": cid, "version": v, "paper": paper["id"], "verdict": vd["text"], "text": c["text"],
                        "questions": c.get("questions", []), "page_pdf": c["page_pdf"], "kind": c["kind"],
                        "basis": c["basis"], "cited_work": c.get("cited_work"), "equation": c.get("equation"),
                        "parameter": c.get("parameter"), "conditions": c.get("conditions", []),
                        "meaning": vd.get("meaning"),
                        "user_verdict": (rec.get("user_decision") or {}).get("verdict")})
        counts[paper["id"]] = {"published": len(latest), "given": kept}
        chk = pdir / "analysis" / Path(run).name / "applicability.json"
        if chk.exists():
            lists[paper["id"]] = json.loads(chk.read_text())["items"]
    return out, lists, counts


def checked_claims(run, question):
    """Development path (b1-dev2): claims straight from a run's read.json and verdicts.json."""
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
    ap.add_argument("--skeptic", action="store_true",
                    help="with --revise: fix the statements Y3 found overstated or unsupported, not Y2's failures")
    ap.add_argument("--source", choices=["library", "run"], default="library",
                    help="library: the published claims (the line); run: read.json and verdicts.json (development)")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    if a.source == "library":
        claims, lists, counts = library_claims(run)
        q = rc.questions()
    else:
        claims, lists, counts = checked_claims(run, a.question)
        q = next(line for line in rc.questions().splitlines() if line.startswith(a.question))
    outdir = run / "untrusted" / "memos"
    outdir.mkdir(parents=True, exist_ok=True)
    if a.revise:
        prev = f"{a.question}.v{a.version - 1}"
        rejected = user_rejected(run)                                       # the draft's claims, minus any
        claims = [c for c in json.loads((outdir / f"{prev}.claims.json").read_text())   # the user rejected since
                  if c["id"] not in rejected]
        draft = (outdir / f"{prev}.md").read_text()
        if a.skeptic:                                   # Y3 found statements the cited claims do not support
            y3 = json.loads((outdir / f"{prev}.y3.json").read_text())
            failures = [f"- {x['verdict']} (section {x['section']}; cites {', '.join(x.get('claim_ids', []))}): "
                        f"{x['summary']} Reason: {x['reason']}"
                        for x in y3["support"]["statements"] if x["verdict"] != "supported"]
            if not failures:
                raise SystemExit(f"Y3 found every statement of v{a.version - 1} supported; nothing to revise")
        else:
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
                            step=("Y1.revise_y3" if a.skeptic else "Y1.revise") if a.revise else "Y1", question=a.question, version=a.version)
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
