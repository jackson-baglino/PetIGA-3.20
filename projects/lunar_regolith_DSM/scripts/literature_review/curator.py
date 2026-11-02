#!/usr/bin/env python3
"""L1: publish one paper's checked claims into the library; the only writer of library/ (WORKFLOW.md §3 L1, §6).

    curator.py publish --run <run folder> --paper-id <id>

From the run it takes the copy (summary_copy_<paper>.json and the A2 gate), P1's content for that copy, the
checker's applicability checklist and the reader's symbol list (analysis records, unchecked), and E5's final.json.
It publishes a claim only if its status is confirmed or revised (the checker found it supported and E3 placed its
quote on the cited page); a disputed or withdrawn claim is refused and counted. The user's U3 decisions
(<run>/u3_decisions.json, u3_review.py) take part: a disputed claim the user marked right is published (E3 must
still have placed its quote), a claim marked wrong or unsure is refused, and every decided claim's record carries
the decision (user_decision). A record published before the user marked it wrong stays (retraction is not built);
synthesize.py leaves it out of the memo. U3 comes after L1 in the line, so stage L1b runs this again with
--after-u3 once U3a is approved. Records are written once:
- library/papers/<paper>/card.vN.json        bibliographic data, copies with kind and the user's A2 decision
- library/papers/<paper>/source/<copy>/      the PDF read
- library/papers/<paper>/content/<copy>/     P1's content, copied once, never edited
- library/papers/<paper>/analysis/<run>/     applicability.json, symbols.json
- library/papers/<paper>/claims/<claim>.vN.json
- library/events.jsonl                       one line per publish: who, what, when, counts, hashes (committed)
A record that already exists with the same content is left alone; one that differs gets the next version, and the
event names both. Nothing is deleted. Retraction (user only) is not built yet. Counts go to
<run>/summary_publish_<paper>.json (summary_publish_u3_<paper>.json with --after-u3); nothing printed holds paper
text.
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_b1  # noqa: E402

LIBRARY = Path(os.environ.get("LR_LIBRARY") or run_b1.PROJECT / "library")     # the variable is for tests only
PUBLISHABLE = ("confirmed", "revised")


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def write_versioned(folder, stem, obj):
    """Write obj as <stem>.vN.json unless the latest version already holds it. Returns (path, 'new'|'same'|'next')."""
    folder.mkdir(parents=True, exist_ok=True)
    vs = sorted(folder.glob(f"{stem}.v*.json"), key=lambda p: int(p.name.rsplit(".v", 1)[1].split(".")[0]))
    if vs:
        last = json.loads(vs[-1].read_text())
        if digest({k: v for k, v in last.items() if k != "published"}) == \
                digest({k: v for k, v in obj.items() if k != "published"}):
            return vs[-1], "same"
        n = int(vs[-1].name.rsplit(".v", 1)[1].split(".")[0]) + 1
        kind = "next"
    else:
        n, kind = 1, "new"
    p = folder / f"{stem}.v{n}.json"
    p.write_text(json.dumps(obj, indent=1, ensure_ascii=False))
    return p, kind


def copy_once(src, dst):
    if dst.exists():
        return False
    if src.is_dir():
        shutil.copytree(src, dst)
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
    return True


def user_decisions(run):
    """{claim_id: decision} from the user's U3 review (u3_review.py); empty before U3."""
    f = Path(run) / "u3_decisions.json"
    return json.loads(f.read_text()) if f.exists() else {}


def publish(run, paper_id, library=None):
    library = Path(library or LIBRARY)
    run = Path(run)
    plist = {p["id"]: p for p in run_b1.papers(check_hash=False)}
    if paper_id not in plist:
        raise SystemExit(f"{paper_id} not in the contract")
    paper = plist[paper_id]
    copy = json.loads((run / f"summary_copy_{paper_id}.json").read_text())
    gate = run / "gates" / f"A2_{paper_id}.json"
    decision = json.loads(gate.read_text()) if gate.exists() else None
    if copy["gate_opened"] and (decision or {}).get("status") != "approved":
        raise SystemExit(f"{paper_id}: copy kind {copy['kind']} has no approved A2 decision")
    h = copy["sha256"]
    content = run / "untrusted" / "content" / paper_id / h[:16]
    if not content.is_dir():
        raise SystemExit(f"{paper_id}: no P1 content for copy {h[:16]}")
    u = run / "untrusted" / paper_id
    final = json.loads((u / "final.json").read_text())
    pdir = library / "papers" / paper_id
    card = {"paper_id": paper_id, "doi": copy["doi"], "role": paper["role"], "branches": ["B1"],
            "copies": [{"sha256": h, "kind": copy["kind"], "pages": copy["pages"],
                        "a2_decision": (decision or {}).get("decision"),
                        "file": f"source/{h[:16]}/{Path(paper['pdf']).name}"}],
            "current_copy": h}
    card_path, card_kind = write_versioned(pdir, "card", card)
    copy_once(Path(paper["pdf"]), pdir / "source" / h[:16] / Path(paper["pdf"]).name)
    content_new = copy_once(content, pdir / "content" / h[:16])
    adir = pdir / "analysis" / run.name
    if not adir.exists():
        adir.mkdir(parents=True)
        shutil.copyfile(u / "checklist.json", adir / "applicability.json")
        (adir / "symbols.json").write_text(json.dumps(json.loads((u / "read.json").read_text())["symbols"],
                                                      indent=1, ensure_ascii=False))
    counts = {"published": 0, "same": 0, "next_version": 0, "refused_disputed": 0, "refused_withdrawn": 0,
              "refused_unchecked": 0, "refused_user_wrong": 0, "refused_user_unsure": 0, "by_user_decision": 0}
    decisions = user_decisions(run)
    published = []
    for f in final:
        st, dec = f["status"], decisions.get(f["claim_id"])
        if dec and dec.get("verdict") in ("wrong", "unsure"):
            counts["refused_user_" + dec["verdict"]] += 1
            continue
        by_user = st == "disputed" and (dec or {}).get("verdict") == "right"
        if st not in PUBLISHABLE and not by_user:
            counts["refused_" + (st if st in ("disputed", "withdrawn") else "unchecked")] += 1
            continue
        if not (f["e3"] or {}).get("ok") or (not by_user and (f["verdict"] or {}).get("text") != "supported"):
            counts["refused_unchecked"] += 1
            continue
        rec = {"claim_id": f["claim_id"], "paper_id": paper_id, "copy": h, "run": run.name, "status": st,
               "claim": f["claim"], "verdict": f["verdict"], "quote_check": f["e3"], "published": now()}
        if dec:
            rec["user_decision"] = dec
        counts["by_user_decision"] += by_user
        p, kind = write_versioned(pdir / "claims", f["claim_id"], rec)
        counts["published" if kind == "new" else "same" if kind == "same" else "next_version"] += 1
        published.append(p.name)
    event = {"time": now(), "event": "publish", "by": "curator.py", "paper_id": paper_id, "run": run.name,
             "copy": h, "card": card_path.name, "card_change": card_kind, "content_copied": content_new,
             "claims": sorted(published), **counts}
    library.mkdir(parents=True, exist_ok=True)
    with open(library / "events.jsonl", "a") as ev:
        ev.write(json.dumps(event) + "\n")
    return counts


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("publish")
    p.add_argument("--run", required=True)
    p.add_argument("--paper-id", required=True)
    p.add_argument("--after-u3", action="store_true",
                   help="stage L1b: publish again with the user's U3 decisions; needs the U3a gate approved")
    a = ap.parse_args(argv)
    run = Path(a.run)
    try:
        if a.after_u3:
            g = run / "gates" / "U3a.json"
            if not g.exists() or json.loads(g.read_text()).get("status") != "approved":
                raise SystemExit("U3a is not approved: no user decisions to publish")
        counts, failed = publish(run, a.paper_id), None
    except (SystemExit, FileNotFoundError) as e:
        counts, failed = {}, str(e)
    summary = {"paper_id": a.paper_id, **counts}
    if failed:
        summary["failed"] = failed
    (run / f"summary_publish{'_u3' if a.after_u3 else ''}_{a.paper_id}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
