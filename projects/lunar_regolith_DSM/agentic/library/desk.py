"""The reference desk: the only way agents get information out of the library.  Every request is checked against the policy and logged.
Deterministic retrieval comes first (SQL, full-text search, hashed files); the librarian-reference agent only composes an answer from what was retrieved.

  python -m agentic.library.desk ask --as human --why "checking R8 claims" claims q="contact angle" limit=5
  python -m agentic.library.desk stats
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from agentic.library.store import Library, Curator, sha256, now

# who may ask for what (an agent not listed gets nothing)
POLICY = {
    "human": {"stats", "list", "paper", "claims", "record", "page", "memo", "errata", "candidates", "log"},
    "topic-synthesizer": {"list", "paper", "claims", "record", "errata"},
    "claim-checker": {"page", "paper", "claims"},
    "image-checker": {"page"},
    "formulation-reviewer": {"paper", "claims", "memo", "errata"},
    "math-spec-writer": {"paper", "claims", "memo"},
    "math-deriver": {"claims"},
    "math-adversary": {"claims", "errata"},
    "report-outliner": {"stats"}, "report-drafter": {"stats"}, "report-critic": set(), "report-rubric-auditor": set(),
    "query-planner": {"list", "candidates", "stats"},
    "librarian-reference": {"claims", "paper", "record", "memo", "errata", "list"},
    "librarian-auditor": {"stats", "log", "errata", "list"},
}


class Denied(PermissionError):
    pass


class Desk(Library):
    def __init__(self, root: Path | None = None):
        super().__init__(root)

    # ---- helpers
    def _log(self, requester: str, purpose: str, kind: str, request: dict, returned: list) -> None:
        cur = Curator(self.root, actor="desk")           # the desk may append to the access log only (an event row), nothing else
        with cur.rw() as con:
            con.execute("INSERT INTO access_log VALUES (?,?,?,?,?,?)", (now(), requester, purpose, kind, json.dumps(request), json.dumps(returned)[:2000]))

    def _read_verified(self, rel: str) -> str:
        con = self.ro()
        row = con.execute("SELECT sha256, version FROM files WHERE rel=?", (rel,)).fetchone(); con.close()
        if not row:
            raise FileNotFoundError(rel)
        p = self.root / rel
        if sha256(p) != row["sha256"]:
            raise RuntimeError(f"library file {rel} does not match its recorded hash: run verify and tell the human")
        return p.read_text()

    # ---- requests
    def ask(self, requester: str, purpose: str, kind: str, **params) -> dict:
        if requester not in POLICY or kind not in POLICY[requester]:
            self._log(requester, purpose, kind + " DENIED", params, [])
            raise Denied(f"{requester!r} may not request {kind!r}")
        if not purpose.strip():
            raise Denied("a purpose is required for every request")
        con = self.ro()
        out: dict = {"kind": kind, "provenance": []}
        try:
            if kind == "stats":
                out["papers"] = con.execute("SELECT COUNT(*) FROM papers").fetchone()[0]; out["claims"] = con.execute("SELECT COUNT(*) FROM claims").fetchone()[0]
                out["flagged_claims"] = con.execute("SELECT COUNT(*) FROM flags").fetchone()[0]; out["candidates"] = con.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
                out["memos"] = con.execute("SELECT COUNT(*) FROM memos").fetchone()[0]
            elif kind == "list":
                q = "SELECT id,year,title,wave,topics,status,checker_rate,n_flagged FROM papers WHERE 1=1"; a = []
                for k in ("status", "wave"):
                    if params.get(k): q += f" AND {k}=?"; a.append(params[k])
                if params.get("topic"): q += " AND topics LIKE ?"; a.append(f"%{params['topic']}%")
                out["papers"] = [dict(r) for r in con.execute(q + " ORDER BY year DESC LIMIT ?", a + [int(params.get("limit", 100))])]
            elif kind == "paper":
                r = con.execute("SELECT * FROM papers WHERE id=?", (params["paper_id"],)).fetchone()
                out["paper"] = dict(r) if r else None
                out["flags"] = [dict(x) for x in con.execute("SELECT claim_id,source,verdict,note FROM flags WHERE paper_id=?", (params["paper_id"],))]
            elif kind == "claims":
                if params.get("paper_id") and params.get("page"):          # every claim on one page (used when checking a page)
                    out["claims"] = [dict(r) for r in con.execute("SELECT paper_id,claim_id,page,claim,quote,flag FROM claims WHERE paper_id=? AND page=?", (params["paper_id"], int(params["page"])))]
                    self._log(requester, purpose, kind, params, []); con.close(); return out
                q = params.get("q", ""); terms = " ".join(re.findall(r"[A-Za-z0-9]+", q))
                if not terms:
                    out["claims"] = []; out["note"] = "empty query"; self._log(requester, purpose, kind, params, []); con.close(); return out
                sql = "SELECT c.paper_id,c.claim_id,c.page,c.claim,c.quote,c.flag,p.title,p.year FROM claims_fts f JOIN claims c ON c.paper_id=f.paper_id AND c.claim_id=f.claim_id JOIN papers p ON p.id=c.paper_id WHERE claims_fts MATCH ?"
                a = [" OR ".join(terms.split())]
                if params.get("paper_id"): sql += " AND c.paper_id=?"; a.append(params["paper_id"])
                sql += " ORDER BY rank LIMIT ?"; a.append(int(params.get("limit", 10)))
                out["claims"] = [dict(r) for r in con.execute(sql, a)]
                out["note"] = "claims with a non-empty 'flag' were judged only partly supported by the independent checker: use the quote, not the claim wording"
            elif kind == "record":
                rel = f"papers/{params['paper_id']}/record.md"; out["text"] = self._read_verified(rel)
                out["flags"] = [dict(x) for x in con.execute("SELECT claim_id,verdict,note FROM flags WHERE paper_id=?", (params["paper_id"],))]
                out["provenance"].append(rel)
            elif kind == "page":
                from agentic.lit.verify_record import split_pages
                text = self._read_verified(f"papers/{params['paper_id']}/text.txt"); pg = split_pages(text)
                out["text"] = pg.get(int(params["page"]), ""); out["provenance"].append(f"papers/{params['paper_id']}/text.txt p.{params['page']}")
            elif kind == "memo":
                r = con.execute("SELECT rel FROM memos WHERE name=?", (params["name"],)).fetchone()
                out["text"] = self._read_verified(r["rel"]) if r else None; out["provenance"].append(r["rel"] if r else "not found")
            elif kind == "errata":
                out["flags"] = [dict(x) for x in con.execute("SELECT * FROM flags WHERE paper_id=?", (params.get("paper_id", ""),))] if params.get("paper_id") else \
                               [dict(x) for x in con.execute("SELECT * FROM flags LIMIT ?", (int(params.get("limit", 50)),))]
            elif kind == "candidates":
                out["candidates"] = [dict(r) for r in con.execute("SELECT id,doi,title,year,best_topic,in_library FROM candidates WHERE title LIKE ? LIMIT ?", (f"%{params.get('q', '')}%", int(params.get("limit", 25))))]
            elif kind == "log":
                out["log"] = [dict(r) for r in con.execute("SELECT ts,requester,purpose,kind FROM access_log ORDER BY ts DESC LIMIT ?", (int(params.get("limit", 50)),))]
        finally:
            con.close()
        self._log(requester, purpose, kind, params, [out.get("provenance")] if out.get("provenance") else [])
        return out

    def ask_natural(self, requester: str, purpose: str, question: str, limit: int = 12) -> dict:
        """Deterministic retrieval, then the librarian-reference agent selects and cites; it may only use what was retrieved."""
        from agentic.lit import common
        hits = self.ask(requester if requester in POLICY else "human", purpose, "claims", q=question, limit=limit)["claims"]
        schema = {"type": "object", "required": ["items", "not_found", "caveats", "human_should_check"], "properties": {
            "items": {"type": "array", "items": {"type": "object", "required": ["paper_id", "claim_id", "page", "quote", "why"], "properties": {k: {"type": ["string", "integer"]} for k in ("paper_id", "claim_id", "page", "quote", "why")}}},
            "not_found": {"type": "boolean"}, "caveats": {"type": "array", "items": {"type": "string"}}, "human_should_check": {"type": "array", "items": {"type": "string"}}}}
        out, _ = common.call_claude("librarian-reference", json.dumps({"requester": requester, "purpose": purpose, "question": question, "retrieved": hits}, indent=1), schema=schema, label="library:reference")
        self._log("librarian-reference", purpose, "answer", {"question": question}, [[i["paper_id"], i["claim_id"]] for i in out["items"]])
        return out


def main() -> None:
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("ask"); a.add_argument("--as", dest="who", required=True); a.add_argument("--why", required=True); a.add_argument("kind"); a.add_argument("params", nargs="*")
    sub.add_parser("stats")
    n = ap.parse_args()
    d = Desk()
    if n.cmd == "stats":
        print(json.dumps(d.ask("human", "stats", "stats"), indent=1)); return
    kv = dict(p.split("=", 1) for p in n.params)
    print(json.dumps(d.ask(n.who, n.why, n.kind, **kv), indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
