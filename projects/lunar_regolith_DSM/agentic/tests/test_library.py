"""Library tests on a throw-away library: only the curator writes, the desk enforces the policy and logs, tampering is detected."""
import json, os, stat, sys, tempfile, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from agentic.library.store import Curator, Library
from agentic.library.desk import Desk, Denied
from agentic.library.verify import verify

TMP = pathlib.Path(tempfile.mkdtemp()) / "Library"
REC = {"paper": {"title": "A test paper", "year": 2020, "doi": "10.1/x"}, "claims": [
    {"id": "c1", "page": 1, "claim": "Contact angle is prescribed on the wall.", "quote": "the contact angle is prescribed on the wall", "type": "method"},
    {"id": "c2", "page": 2, "claim": "Error scales with width.", "quote": "error scales with the interface width", "type": "result"}]}
VER = {"metrics": {"quote_match_rate": 1.0}, "checker": {"supported_rate": 0.5, "flagged": {"c2": {"verdict": "partially_supported", "problem": "only in 1D"}}}}

def build():
    cur = Curator(TMP, actor="test")
    txt = TMP.parent / "t.txt"; txt.write_text("=== PAGE 1 ===\nthe contact angle is prescribed on the wall\n=== PAGE 2 ===\nerror scales with the interface width\n")
    cur.publish_paper("test_2020_x", REC, VER, "# test_2020_x\n", txt, None, wave="w", topics="R4")
    return Desk(TMP)

def test_publish_and_read_only():
    d = build()
    assert not os.access(TMP / "papers/test_2020_x/record.json", os.W_OK)
    assert not os.access(TMP / "catalog.sqlite", os.W_OK)
    import sqlite3
    con = Library(TMP).ro()
    try: con.execute("INSERT INTO events VALUES ('a','b','c','d','e')"); raise AssertionError("catalog should be read-only")
    except sqlite3.OperationalError: pass

def test_desk_policy_and_log():
    d = build()
    assert d.ask("human", "testing", "stats")["papers"] == 1
    r = d.ask("topic-synthesizer", "write memo", "claims", q="contact angle")
    assert r["claims"][0]["claim_id"] == "c1"
    for who, kind in (("image-checker", "claims"), ("stranger", "stats"), ("report-critic", "stats")):
        try: d.ask(who, "x", kind); raise AssertionError((who, kind))
        except Denied: pass
    try: d.ask("human", "", "stats"); raise AssertionError("purpose required")
    except Denied: pass
    log = d.ask("human", "see log", "log")["log"]
    assert any("DENIED" in l["kind"] for l in log)

def test_flags_visible():
    d = build()
    c = d.ask("human", "t", "claims", q="error width")["claims"]
    assert c[0]["flag"] == "partially_supported"

def test_page_read_and_tamper_detection():
    d = build()
    assert "contact angle" in d.ask("claim-checker", "check c1", "page", paper_id="test_2020_x", page=1)["text"]
    assert verify(Library(TMP))["ok"]
    p = TMP / "papers/test_2020_x/text.txt"; p.parent.chmod(0o755); p.chmod(0o644); p.write_text("tampered")
    assert not verify(Library(TMP))["ok"]
    try: d.ask("claim-checker", "check", "page", paper_id="test_2020_x", page=1); raise AssertionError("tampered file must be refused")
    except RuntimeError: pass

def test_published_paper_is_frozen_and_versions_kept():
    from agentic.library.store import Frozen
    d = build()
    cur = Curator(TMP, actor="test"); REC2 = dict(REC, claims=REC["claims"][:1])
    try: cur.publish_paper("test_2020_x", REC2, VER, None, None, None); raise AssertionError("must be frozen")
    except Frozen: pass
    try: cur.publish_paper("test_2020_x", REC2, VER, None, None, None, supersede_reason="typo", human_approved=False); raise AssertionError("needs human approval")
    except Frozen: pass
    cur.publish_paper("test_2020_x", REC2, VER, None, None, None, supersede_reason="fix", human_approved=True)
    assert (TMP / "history/papers/test_2020_x/record.json.v1").exists()

if __name__ == "__main__":
    f = 0
    for n, fn in list(globals().items()):
        if n.startswith("test_"):
            try:
                import shutil
                if TMP.exists():
                    os.system(f"chmod -R u+w '{TMP}'"); shutil.rmtree(TMP)
                fn(); print("PASS", n)
            except Exception as e: f += 1; print("FAIL", n, type(e).__name__, e)
    sys.exit(1 if f else 0)
