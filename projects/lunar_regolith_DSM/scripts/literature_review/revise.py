#!/usr/bin/env python3
"""E5: revise the claims the checker did not confirm, at most two rounds (WORKFLOW.md §3, E5).

    revise.py --run <run folder> --paper-id <id> --pdf <file>

Reads read.json, e3.json, checklist.json and verdicts.json from <run>/untrusted/<paper>/ (written by read_check.py).
A claim is confirmed when the checker found its text supported, its quote on the page, its equation and value
unchanged and no dropped condition, and E3 placed its quote on the cited page. Each round, the other claims go to a
fresh morty-reader (E5 task: the claim, the checker's verdict, the cited and condition pages), which revises or
withdraws each one; a revision is quote-checked (E3) and checked by a fresh beth-checker (E4 step 2, with the
paper's checklist). After two rounds a claim still not confirmed is *disputed*: it goes to the user's review queue
(U3) and is not published unless the user decides it. The blind transcription of a revised equation or value is
left to X1, which reads again every equation or value a memo cites.

Writes <run>/untrusted/<paper>/final.json (each claim's latest version, verdict, quote check and status:
confirmed, revised, disputed or withdrawn) and revisions.json (every round), and the counts in
<run>/summary_revise_<paper>.json. Prints counts only.
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call as ac  # noqa: E402
import read_check as rc  # noqa: E402

ROUNDS = 2


def confirmed(verdict, e3):
    if verdict is None or e3 is None:
        return False
    return (verdict["text"] == "supported" and verdict["quote_on_page"] and verdict["equation"] in ("matches", "none")
            and verdict["parameter"] in ("match", "none") and not verdict["dropped_conditions"] and e3["ok"])


def pages_for(claims):
    cited = sorted({c["page_pdf"] for c in claims})
    cond = sorted({k["page_pdf"] for c in claims for k in c.get("conditions", [])} - set(cited))
    return cited + cond


def groups(claims, size):
    """Claims grouped by their cited pages, at most `size` cited pages per group."""
    cited = sorted({c["page_pdf"] for c in claims})
    for i in range(0, len(cited), size):
        g = set(cited[i:i + size])
        yield [c for c in claims if c["page_pdf"] in g]


def add_review(run, items):
    f = Path(run) / "review_queue.json"
    q = json.loads(f.read_text()) if f.exists() else []
    known = {(x["stage"], x.get("claim_id")) for x in q}
    q += [x for x in items if (x["stage"], x.get("claim_id")) not in known]
    f.write_text(json.dumps(q, indent=1))


def run_rounds(state, by_num, textdir, checklist, call, rounds=ROUNDS):
    """state: {claim_id: {claim, verdict, e3, status, history}}. call(kind, files, text, schema, agent, **tag) runs
    one agent call in a folder holding `files` and returns its output or raises. Mutates state; returns rounds run."""
    done = 0
    for r in range(1, rounds + 1):
        failing = [s for s in state.values() if s["status"] == "pending"]
        if not failing:
            break
        done = r
        for grp in groups([s["claim"] for s in failing], rc.GROUP):
            shown = [{"claim": s["claim"], "verdict": state[s["claim"]["id"]]["verdict"],
                      "quote_found_on_cited_page": bool((state[s["claim"]["id"]]["e3"] or {}).get("ok"))}
                     for s in (state[c["id"]] for c in grp)]
            show = [by_num[p] for p in pages_for(grp)[:rc.MAX_IMAGES] if p in by_num]
            res = call(show, ac.prompt("reader_revise", round=str(r), files="{files}", claims=json.dumps(shown)),
                       "revise", "morty-reader", step="E5.read", round=r)
            got = {x["claim_id"]: x for x in res["revisions"]}
            revised = []
            for c in grp:
                s = state[c["id"]]
                x = got.get(c["id"])
                if x is None:
                    s["history"].append({"round": r, "action": "no_answer"})
                    continue
                if x["action"] == "withdrawn" or not x.get("claim"):
                    s["status"] = "withdrawn"
                    s["history"].append({"round": r, "action": "withdrawn", "reason": x["reason"]})
                    continue
                new = dict(x["claim"], id=c["id"], version=s["claim"].get("version", 1) + 1)
                s["history"].append({"round": r, "action": "revised", "reason": x["reason"], "before": s["claim"]})
                s["claim"] = new
                revised.append(new)
            if not revised:
                continue
            e3 = {q["claim_id"]: q for q in rc.check_quotes(revised, textdir)}
            for c in revised:
                state[c["id"]]["e3"] = e3[c["id"]]
            for g2 in groups(revised, rc.GROUP):
                cs = [{k: c[k] for k in c if k != "symbols"} for c in g2]
                for c in cs:
                    c["quote_check"] = "ok" if e3[c["id"]]["ok"] else "failed"
                show = [by_num[p] for p in pages_for(g2)[:rc.MAX_IMAGES] if p in by_num]
                res = call(show, ac.prompt("checker_verdicts", files="{files}",
                                           checklist=json.dumps(checklist["items"]), claims=json.dumps(cs)),
                           "verdicts", "beth-checker", step="E5.check", round=r)
                for v in res["verdicts"]:
                    if v["claim_id"] in state:
                        state[v["claim_id"]]["verdict"] = v
            for c in revised:
                s = state[c["id"]]
                if confirmed(s["verdict"], s["e3"]):
                    s["status"] = "revised"
    for s in state.values():
        if s["status"] == "pending":
            s["status"] = "disputed"
    return done


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--paper-id", required=True)
    ap.add_argument("--pdf", required=True)
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    out = run / "untrusted" / a.paper_id
    ledger = out / "calls.jsonl"

    def load(name):
        f = out / f"{name}.json"
        if not f.exists():
            raise SystemExit(f"{a.paper_id}: no {name}.json; run E1 to E4 first")
        return json.loads(f.read_text())

    claims, e3, checklist = load("read")["claims"], load("e3"), load("checklist")
    verdicts = {v["claim_id"]: v for v in load("verdicts")}
    e3 = {q["claim_id"]: q for q in e3}
    state = {}
    for c in claims:
        v, q = verdicts.get(c["id"]), e3.get(c["id"])
        state[c["id"]] = {"claim": dict(c, version=1), "verdict": v, "e3": q, "history": [],
                          "status": "confirmed" if confirmed(v, q) else "pending"}
    pages = sorted((out / "pages").glob("p*.png"))
    by_num = {int(f.stem[1:]): f for f in pages}

    def call(show, text, schema, agent, **tag):
        with tempfile.TemporaryDirectory() as tmp:
            files = rc.folder_with(show, tmp)
            res, meta = ac.call(tmp, text.replace("{files}", files), schema, agent, ledger=ledger, **tag)
        if res is None:
            raise SystemExit(f"{tag.get('step')} failed: {meta.get('error')}")
        return res

    try:
        rounds = run_rounds(state, by_num, out / "text", checklist, call)
        failed = None
    except SystemExit as e:
        rounds, failed = None, str(e)
    final = [{"claim_id": k, "status": s["status"], "claim": s["claim"], "verdict": s["verdict"], "e3": s["e3"]}
             for k, s in state.items()]
    if failed is None:
        (out / "final.json").write_text(json.dumps(final, indent=1))
        (out / "revisions.json").write_text(json.dumps({k: s["history"] for k, s in state.items() if s["history"]},
                                                       indent=1))
        add_review(run, [{"stage": "E5", "item": a.paper_id, "claim_id": f["claim_id"], "kind": "disputed"}
                         for f in final if f["status"] == "disputed"])
    counts = {st: sum(1 for f in final if f["status"] == st)
              for st in ("confirmed", "revised", "disputed", "withdrawn", "pending")}
    summary = {"paper_id": a.paper_id, "claims": len(final), **counts, "rounds": rounds,
               "calls": len([m for m in ac.ledger_calls(ledger) if str(m.get("step", "")).startswith("E5")])}
    if failed:
        summary["failed"] = failed
    (run / f"summary_revise_{a.paper_id}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
