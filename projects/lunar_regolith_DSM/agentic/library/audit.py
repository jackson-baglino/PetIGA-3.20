"""Library audit for the human: integrity check + access-log summary, explained in plain language by the librarian-auditor agent.
Writes the report OUTSIDE the library (research/Review/library_audit_<date>.md); the auditor cannot change the library.
  python -m agentic.library.audit"""
from __future__ import annotations

import collections
import json
import time
from pathlib import Path

from agentic.library.store import Library, ROOT
from agentic.library.verify import verify


def summary(lib: Library) -> dict:
    con = lib.ro()
    log = con.execute("SELECT requester, kind, purpose, request FROM access_log").fetchall()
    by = collections.Counter((r["requester"], r["kind"].replace(" DENIED", "")) for r in log)
    denied = collections.Counter((r["requester"], r["kind"]) for r in log if "DENIED" in r["kind"])
    repeats = collections.Counter((r["requester"], r["kind"], r["request"]) for r in log)
    papers = con.execute("SELECT id,status,n_flagged,checker_rate FROM papers").fetchall()
    memos = [(m["name"], (lib.root / m["rel"]).read_text() if m["rel"].endswith(".md") else "") for m in con.execute("SELECT name, rel FROM memos")]
    flagged_papers = {p["id"]: p["n_flagged"] for p in papers if p["n_flagged"]}
    memo_cites = {n: sorted(s for s in flagged_papers if s in t)[:8] for n, t in memos if any(s in t for s in flagged_papers)}
    out = {"requests_total": len(log), "requests_by_requester_and_kind": {f"{a}/{b}": n for (a, b), n in by.most_common(15)},
           "denied": {f"{a}/{b}": n for (a, b), n in denied.items()}, "repeated_identical_requests": [f"{a}/{b}: {n}x" for (a, b, _), n in repeats.items() if n >= 5][:10],
           "flagged_claims_total": con.execute("SELECT COUNT(*) FROM flags").fetchone()[0], "papers_with_flags": len(flagged_papers),
           "papers_on_hold": [p["id"] for p in papers if p["status"] == "hold"], "memos_citing_papers_with_flags": memo_cites}
    con.close(); return out


def main() -> None:
    from agentic.lit import common
    lib = Library(); data = {"integrity": verify(lib), "access": summary(lib)}
    schema = {"type": "object", "required": ["summary", "problems", "unusual_access", "human_should_check"], "properties": {
        "summary": {"type": "string"}, "unusual_access": {"type": "array", "items": {"type": "string"}}, "human_should_check": {"type": "array", "items": {"type": "string"}},
        "problems": {"type": "array", "items": {"type": "object", "required": ["severity", "what", "why_it_matters", "where_to_look"], "properties": {k: {"type": "string"} for k in ("severity", "what", "why_it_matters", "where_to_look")}}}}}
    out, _ = common.call_claude("librarian-auditor", json.dumps(data, indent=1), schema=schema, label="library:audit")
    lines = [f"# Library audit, {time.strftime('%Y-%m-%d')}", "", out["summary"], "", f"Integrity: {data['integrity']['files_checked']} files checked, {len(data['integrity']['problems'])} problems.", ""]
    lines += [f"- **[{p['severity']}] {p['what']}**: {p['why_it_matters']} Look at: {p['where_to_look']}" for p in out["problems"]]
    lines += ["", "## Unusual access"] + [f"- {u}" for u in out["unusual_access"]] + ["", "## Worth your own eyes"] + [f"- {h}" for h in out["human_should_check"]]
    dest = ROOT / "research" / "Review" / f"library_audit_{time.strftime('%Y-%m-%d')}.md"; dest.write_text("\n".join(lines) + "\n"); print(dest)


if __name__ == "__main__":
    main()
