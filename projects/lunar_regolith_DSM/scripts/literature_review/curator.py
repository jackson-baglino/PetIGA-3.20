#!/usr/bin/env python3
"""L1: publish one paper's checked claims into the library; the only writer of library/ (WORKFLOW.md §3 L1, §6).

    curator.py publish --run <run folder> --paper-id <id>

From the run it takes the copy (summary_copy_<paper>.json and the A2 gate), P1's content for that copy, the
checker's applicability checklist and the reader's symbol list (analysis records, unchecked), and E5's final.json.
It publishes a claim only if its status is confirmed or revised (the checker found it supported and E3 placed its
quote on the cited page); a disputed or withdrawn claim is refused and counted. Records are written once:
- library/papers/<paper>/card.vN.json        bibliographic data, copies with kind and the user's A2 decision
- library/papers/<paper>/source/<copy>/      the PDF read
- library/papers/<paper>/content/<copy>/     P1's content, copied once, never edited
- library/papers/<paper>/analysis/<run>/     applicability.json, symbols.json
- library/papers/<paper>/claims/<claim>.vN.json
- library/events.jsonl                       one line per publish: who, what, when, counts, hashes (committed)
A record that already exists with the same content is left alone; one that differs gets the next version, and the
event names both. Nothing is deleted. Retraction (user only) is not built yet. Counts go to
<run>/summary_publish_<paper>.json; nothing printed holds paper text.
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
              "refused_unchecked": 0}
    published = []
    for f in final:
        st = f["status"]
        if st not in PUBLISHABLE:
            counts["refused_" + (st if st in ("disputed", "withdrawn") else "unchecked")] += 1
            continue
        if not (f["e3"] or {}).get("ok") or (f["verdict"] or {}).get("text") != "supported":
            counts["refused_unchecked"] += 1
            continue
        rec = {"claim_id": f["claim_id"], "paper_id": paper_id, "copy": h, "run": run.name, "status": st,
               "claim": f["claim"], "verdict": f["verdict"], "quote_check": f["e3"], "published": now()}
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
    a = ap.parse_args(argv)
    run = Path(a.run)
    try:
        counts, failed = publish(run, a.paper_id), None
    except (SystemExit, FileNotFoundError) as e:
        counts, failed = {}, str(e)
    summary = {"paper_id": a.paper_id, **counts}
    if failed:
        summary["failed"] = failed
    (run / f"summary_publish_{a.paper_id}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
