#!/usr/bin/env python3
"""Y1: draft the memo from the published claims, and check its citations (WORKFLOW.md §3).

    synthesize.py --run <run folder> [--question Q1] [--version 1] [--revise] [--source library|run]
    synthesize.py --run <run folder> [--question Q1] --after-u3        stage Y1b

In the line (--source library, the default) the claims are the latest published version of every claim on the
contract's papers (curator.py), for Q1 to Q3, minus any the user rejected in U3 or X1 (u3_decisions.json); the
checklists are the papers' applicability records from this run. The memo answers Q1 and has short Q2 and Q3
sections and a Leads list. With --source run (development, b1-dev2) claims come straight from a run's read.json
and verdicts.json: supported, or overstated with a corrected text, and tagged with --question.

With --revise, version k is a revision of version k-1: the synthesizer gets its own draft, the sentences Y2
failed and the same claims, and fixes only those (WORKFLOW.md §3, Y2 failure -> Y1). The memo goes to
<run>/untrusted/memos/<question>.v<k>.md; Y2's result (counts only) is printed and saved.

With --after-u3 (stage Y1b, after L1b has published the user's claim review): the latest memo is compared with the
library. A claim the memo cites that the user rejected (wrong or unsure) is withdrawn; a disputed claim the user
marked right, which the draft never had, is confirmed. If there is neither, nothing is called. Otherwise the
synthesizer revises the latest version into the next one, fixing only those, from the claims the draft cited plus
the confirmed ones (so the revision brings in no claim that X1 or the user has not seen). Then: a prose claim the
revision now cites for a model option that the user has not reviewed opens the Y1b gate (u3_review.py lists it); an
equation or value claim cited that X1 did not read and the user did not decide fails the stage. Counts go to
<run>/summary_redraft_<question>.json.
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
            if cid in rejected or rec.get("status") == "rejected_by_user" or \
                    not set(c.get("questions", [])) & set(QUESTIONS):
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
    ap.add_argument("--after-u3", action="store_true",
                    help="stage Y1b: revise the latest memo for the user's claim review (withdrawn, confirmed claims)")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    if a.after_u3:
        return after_u3(run, a.question)
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
    step = ("Y1.revise_y3" if a.skeptic else "Y1.revise") if a.revise else "Y1"
    print(json.dumps(write_memo(run, a.question, a.version, text, claims, counts, a.revise, step), indent=1))


def cited_ids(memo_md):
    return sorted(set(re.findall(r"\[([a-z0-9_]+-c\d{3})", memo_md)))


def write_memo(run, question, version, text, claims, counts, revise, step):
    """Call the synthesizer, save memo version `version` and its claims, check it (Y2), save and return the summary."""
    outdir = run / "untrusted" / "memos"
    (outdir / f"{question}.v{version}.claims.json").write_text(json.dumps(claims, indent=1))
    with tempfile.TemporaryDirectory() as tmp:
        res, meta = ac.call(tmp, text, "memo", "beth-synthesizer", ledger=run / "calls.jsonl",
                            step=step, question=question, version=version)
    if res is None:
        raise SystemExit(f"Y1 failed: {meta}")
    (outdir / f"{question}.v{version}.md").write_text(res["memo"])
    y2 = cx.check(res["memo"], claims)
    cited = cited_ids(res["memo"])
    summary = {"question": question, "version": version, "revision_of": version - 1 if revise else None,
               "claims_given": len(claims), "per_paper": counts,
               "claims_cited": len(cited), "papers_cited": len({c.rsplit("-c", 1)[0] for c in cited}),
               "y2": {k: (len(v) if isinstance(v, list) else v) for k, v in y2.items()},
               "y2_failing_sentences": {k: v for k, v in y2.items() if isinstance(v, list) and v},
               "seconds": meta["seconds"], "usage": meta["usage"]}
    (run / f"summary_memo_{question}.v{version}.json").write_text(json.dumps(summary, indent=1))
    return summary


def after_u3(run, question):
    """Stage Y1b: revise the latest memo for the user's claim review, if it changed what the memo may cite."""
    import open_gate
    import u3_review
    outdir = run / "untrusted" / "memos"
    k = cx.latest_version(outdir, question)
    memo = (outdir / f"{question}.v{k}.md").read_text()
    given = json.loads((outdir / f"{question}.v{k}.claims.json").read_text())
    current, lists, counts = library_claims(run)
    cur = {c["id"]: c for c in current}
    df = run / "u3_decisions.json"
    decisions = json.loads(df.read_text()) if df.exists() else {}
    cited, given_ids = cited_ids(memo), {c["id"] for c in given}
    withdrawn = [c for c in cited if c in given_ids and c not in cur]
    confirmed = sorted(c for c in cur if c not in given_ids and (decisions.get(c) or {}).get("verdict") == "right")
    out = {"question": question, "from_version": k, "withdrawn": len(withdrawn), "confirmed": len(confirmed),
           "redrafted": False}
    if withdrawn or confirmed:
        pool = [cur[c] for c in cited if c in cur] + [cur[c] for c in confirmed]
        failures = [f"- withdrawn [{c}]: the user rejected this claim in review; it is no longer a checked claim."
                    for c in withdrawn]
        failures += [f"- confirmed [{c}]: the user confirmed this claim in review; you did not have it when you "
                     "drafted." for c in confirmed]
        text = ac.prompt("synthesizer_revise", question=rc.questions(), memo=memo, failures="\n".join(failures),
                         claims=json.dumps(pool), checklists=json.dumps(lists))
        s = write_memo(run, question, k + 1, text, pool, counts, True, "Y1b")
        new = cited_ids((outdir / f"{question}.v{k + 1}.md").read_text())
        comp = run / "untrusted" / "x1" / "compare.json"
        read = {r["claim_id"] for r in json.loads(comp.read_text())} if comp.exists() else set()
        unread = [c for c in new if c in cur and (cur[c].get("equation") or cur[c].get("parameter"))
                  and c not in read and c not in decisions]
        waiting = [c for _, c in u3_review.items(run, question) if c not in decisions]
        out.update({"redrafted": True, "version": k + 1, "claims_cited": s["claims_cited"], "y2_pass": s["y2"]["pass"],
                    "x1_unread": len(unread), "u3_waiting": len(waiting)})
        if unread:
            out["failed"] = f"the revision cites {len(unread)} equation or value claim(s) that X1 did not read"
        elif waiting:
            open_gate.open_gate(run, "Y1b", f"The revised memo cites {len(waiting)} claim(s) you have not reviewed: "
                                "run u3_review.py, then approve Y1b", count=len(waiting))
    (run / f"summary_redraft_{question}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    return 1 if out.get("failed") else 0


if __name__ == "__main__":
    sys.exit(main())
