#!/usr/bin/env python3
"""U3 (and the user's side of X1): the user reviews claims in their own terminal (WORKFLOW.md §3 U3, §5).

    venv_workflows/bin/python scripts/literature_review/u3_review.py --run <run folder> [--question Q1]
    ... --list        counts only: what is waiting, nothing from the papers

Run it in a terminal of your own, not inside Claude Code: it prints quotes from the papers, and output there
would enter the session (it refuses when CLAUDECODE is set). It can be stopped and resumed; decisions are saved
after each answer.

Items, in this order:
1. X1 disagreements: the reads of an equation or value side by side; is the claim's (the reader's) read right?
2. Prose claims a model option rests on (claims cited in the memo's Model options section that are not equation
   or value claims): the claim, its quote and page, and the checker's four meaning fields (what was established,
   by whom, under which assumptions, with what evidence). Are all four right?
3. Disputed claims (E5): the claim and the checker's last reason. Is the claim right?
Answers: r (right), w (wrong), u (unsure), s (skip for now), q (quit); then an optional note. A claim marked
wrong is left out of the memo when Y1 redrafts; unsure is kept out too, and noted. At the end you may add
statements the claims missed (paper, PDF page, statement); they are recorded in u3_missed.json (turning them into
claims through E3 and E4 is not built yet).

--redo CLAIM_ID ... asks again about claims already decided (a mislabeled answer); the earlier answer is kept in the
decision's "previous" list, and stage L1b is marked not done so the line publishes again with the new answer.

Decisions go to <run>/u3_decisions.json and <run>/u3_missed.json. When every item is decided, the tool prints the
run_b1.py approve commands for the U3a and X1 gates; you run them.
"""
import argparse
import datetime
import json
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_citations as cx  # noqa: E402
import synthesize as sy  # noqa: E402

ANSWERS = {"r": "right", "w": "wrong", "u": "unsure"}


def model_option_claims(memo_md):
    """Claim IDs cited in the memo's Model options section."""
    ids, inside = [], False
    for line in memo_md.splitlines():
        if line.startswith("#"):
            inside = bool(re.match(r"^#+\s*(\d+\.\s*)?\**\s*model options", line, re.I))
            continue
        if inside:
            ids += re.findall(r"\[([a-z0-9_]+-c\d{3})", line)
    return list(dict.fromkeys(ids))


def items(run, question):
    """[(kind, claim_id)] waiting for the user, in review order."""
    run = Path(run)
    mdir = run / "untrusted" / "memos"
    k = cx.latest_version(mdir, question)
    memo = (mdir / f"{question}.v{k}.md").read_text()
    claims = {c["id"]: c for c in sy.library_claims(run)[0]}
    out = []
    comp = run / "untrusted" / "x1" / "compare.json"
    if comp.exists():
        out += [("x1", r["claim_id"]) for r in json.loads(comp.read_text()) if not r["agree"]]
    out += [("model_option", cid) for cid in model_option_claims(memo)
            if cid in claims and not (claims[cid].get("equation") or claims[cid].get("parameter"))]
    q = run / "review_queue.json"
    if q.exists():
        out += [("disputed", x["claim_id"]) for x in json.loads(q.read_text())
                if x["stage"] == "E5" and x["kind"] == "disputed"]
    seen, uniq = set(), []
    for kind, cid in out:
        if (kind, cid) not in seen:
            seen.add((kind, cid))
            uniq.append((kind, cid))
    return uniq


def dispute_reasons(f):
    """Each test of revise.confirmed() that the claim's final verdict or quote check fails, in words."""
    v, e3, out = f.get("verdict"), f.get("e3"), []
    if v is None:
        return ["no verdict from the checker"]
    if v.get("text") != "supported":
        out.append(f"statement judged {v.get('text')}" + (f"; checker's wording: {v['corrected_text']}"
                                                          if v.get("corrected_text") else ""))
    if not v.get("quote_on_page"):
        out.append("the checker did not find the quote on the cited page")
    if v.get("equation") not in ("matches", "none"):
        out.append(f"equation {v.get('equation')}: {v.get('corrected_latex', '')}")
    if v.get("parameter") not in ("match", "none"):
        out.append(f"parameter {v.get('parameter')}: {v.get('corrected_parameter', '')}")
    for d in v.get("dropped_conditions") or []:
        out.append(f"dropped condition: {d if isinstance(d, str) else json.dumps(d, ensure_ascii=False)}")
    if e3 is None:
        out.append("quote check (E3) not run on the final version")
    elif not e3.get("ok"):
        out.append(f"quote check (E3) failed: {e3.get('words')} words (8 needed), found on PDF page(s) "
                   f"{e3.get('found_on')} (needed: the cited page only)")
    return out or ["none of the checks failed on the final version (a disputed status here is unexpected)"]


def load(path, default):
    return json.loads(path.read_text()) if path.exists() else default


def show(kind, cid, run, claims):
    c = claims.get(cid, {})
    print(f"\n=== {kind}: {cid} (paper {c.get('paper', '?')}, PDF page {c.get('page_pdf', '?')}) ===")
    if kind == "x1":
        rows = {r["claim_id"]: r for r in json.loads((Path(run) / "untrusted" / "x1" / "compare.json").read_text())}
        for src, val in rows[cid]["reads"].items():
            print(f"  {src:13} {val}")
        print("Is the claim's read (reader) what the page prints?")
        return
    if kind == "disputed":
        final = json.loads((Path(run) / "untrusted" / c.get("paper", cid.rsplit("-c", 1)[0]) / "final.json")
                           .read_text())
        f = next((x for x in final if x["claim_id"] == cid), {})
        cl, v = f.get("claim", {}), f.get("verdict") or {}
        print(f"  (paper {cid.rsplit('-c', 1)[0]}, PDF page {cl.get('page_pdf', '?')})")
        print(f"  claim:  {cl.get('text')}\n  quote:  {cl.get('quote')}\n  checker: statement {v.get('text')}: "
              f"{v.get('reason')}")
        print("  Why it is disputed (what failed after the revision rounds):")
        for r in dispute_reasons(f):
            print(f"   - {r}")
        print("Is the claim right, as stated, on that page? (Judge the claim, including the points above.)")
        return
    print(f"  claim:  {c.get('text')}")
    m = c.get("meaning") or {}
    print(f"  established: {m.get('established')}\n  by whom:     {m.get('by_whom')} {m.get('cited_work') or ''}"
          f"\n  assumptions: {m.get('assumptions')}\n  evidence:    {m.get('evidence')}")
    print(f"  conditions:  {[k.get('text') for k in c.get('conditions', [])]}")
    print("Are all four meaning fields right?")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--question", default="Q1")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--redo", nargs="+", default=[], metavar="CLAIM_ID",
                    help="ask again about these decided claims; the earlier answer is kept under 'previous'")
    a = ap.parse_args(argv)
    run = Path(a.run).resolve()
    dfile, mfile = run / "u3_decisions.json", run / "u3_missed.json"
    decisions = load(dfile, {})
    unknown = [c for c in a.redo if c not in decisions]
    if unknown:
        raise SystemExit(f"--redo: not decided yet, nothing to redo: {unknown}")
    todo = [(k, c) for k, c in items(run, a.question) if c not in decisions or c in a.redo]
    if a.list:
        print(json.dumps({"waiting": len(todo), "decided": len(decisions),
                          "by_kind": {k: sum(1 for x, _ in todo if x == k) for k in ("x1", "model_option", "disputed")}}))
        return 0
    if os.environ.get("CLAUDECODE"):
        raise SystemExit("Run this in your own terminal, not inside Claude Code: it prints paper text.")
    claims = {c["id"]: c for c in sy.library_claims(run)[0]}
    for kind, cid in todo:
        show(kind, cid, run, claims)
        ans = ""
        while ans not in ("r", "w", "u", "s", "q"):
            ans = input("[r]ight / [w]rong / [u]nsure / [s]kip / [q]uit: ").strip().lower()[:1]
        if ans == "q":
            break
        if ans == "s":
            continue
        note = input("note (optional): ").strip()
        rec = {"verdict": ANSWERS[ans], "kind": kind, "note": note, "by": "user",
               "decided": datetime.datetime.now().astimezone().isoformat(timespec="seconds")}
        if cid in decisions:                                        # --redo: keep the earlier answer
            prev = decisions[cid]
            rec["previous"] = prev.get("previous", []) + [{k: v for k, v in prev.items() if k != "previous"}]
            for f in run.glob("summary_publish_u3_*.json"):           # L1b publishes again with the new answer
                f.unlink()
        decisions[cid] = rec
        dfile.write_text(json.dumps(decisions, indent=1, ensure_ascii=False))
    missed = load(mfile, [])
    while input("\nAdd a statement the claims missed? [y/N]: ").strip().lower().startswith("y"):
        missed.append({"paper_id": input("  paper id: ").strip(), "page_pdf": int(input("  PDF page: ").strip()),
                       "statement": input("  statement: ").strip(),
                       "added": datetime.date.today().isoformat()})
        mfile.write_text(json.dumps(missed, indent=1, ensure_ascii=False))
    left = [(k, c) for k, c in items(run, a.question) if c not in decisions]
    print(f"\n{len(decisions)} decided, {len(left)} waiting, {len(missed)} missed statement(s) recorded.")
    if not left:
        py = "venv_workflows/bin/python scripts/literature_review/run_b1.py"
        print("All items decided. Record the gates:")
        print(f'  {py} approve --run {run} U3a --decision "U3 done: see u3_decisions.json"')
        if (run / "gates" / "X1.json").exists():
            print(f'  {py} approve --run {run} X1 --decision "X1 disagreements decided in u3_decisions.json"')
    return 0


if __name__ == "__main__":
    sys.exit(main())
