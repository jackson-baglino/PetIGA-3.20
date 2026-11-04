"""Y1 reads the library (latest version, minus user rejections); Y1b redrafts after the claim review; Y2's run mode
checks the latest memo."""
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import check_citations as cx  # noqa: E402
import curator  # noqa: E402
import run_b1  # noqa: E402
import synthesize as sy  # noqa: E402


def rec(cid, text, q=("Q1",)):
    return {"claim_id": cid, "claim": {"id": cid, "text": text, "questions": list(q), "page_pdf": 2, "kind": "result",
                                       "basis": "derived", "conditions": []},
            "verdict": {"text": "supported"}}


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        (root / "r" / "Literature").mkdir(parents=True)
        (root / "r" / "Literature" / "pa.pdf").write_bytes(b"a")
        h = hashlib.sha256(b"a").hexdigest()[:12]
        contract = root / "c.md"
        contract.write_text("## Papers\n\n| ID | DOI | File | Pages | sha | Status | Role |\n|---|---|---|---|---|---|---|\n"
                            f"| pa | 10.1/a | `Literature/pa.pdf` | 1 | {h} | ok | core |\n\n## Next\n")
        self.lib = root / "library"
        cl = self.lib / "papers" / "pa" / "claims"
        cl.mkdir(parents=True)
        (cl / "pa-c001.v1.json").write_text(json.dumps(rec("pa-c001", "old")))
        (cl / "pa-c001.v2.json").write_text(json.dumps(rec("pa-c001", "new")))
        (cl / "pa-c002.v1.json").write_text(json.dumps(rec("pa-c002", "rejected one", ("Q2",))))
        (cl / "pa-c003.v1.json").write_text(json.dumps(rec("pa-c003", "off topic", ("Q9",))))
        self.run = root / "run"
        self.run.mkdir()
        (self.run / "u3_decisions.json").write_text(json.dumps({"pa-c002": {"verdict": "wrong"}}))
        self.p = [mock.patch.object(run_b1, "CONTRACT", contract), mock.patch.object(run_b1, "RESEARCH", root / "r"),
                  mock.patch.object(curator, "LIBRARY", self.lib)]
        for x in self.p:
            x.start()

    def tearDown(self):
        for x in self.p:
            x.stop()
        self.tmp.cleanup()


class Memo(Base):
    def test_library_claims(self):
        claims, _, counts = sy.library_claims(self.run)
        self.assertEqual([(c["id"], c["text"], c["version"]) for c in claims], [("pa-c001", "new", 2)])
        self.assertEqual(counts["pa"], {"published": 3, "given": 1})

    def test_y2_run_mode_checks_latest_memo(self):
        m = self.run / "untrusted" / "memos"
        m.mkdir(parents=True)
        claims = [{"id": "pa-c001", "text": "new"}, {"id": "pa-c002", "text": "rejected one"}]
        for k, body in ((1, "## 1. Summary\nThe method works well here [pa-c002].\n"),
                        (2, "## 1. Summary\nThe method works well here [pa-c001].\n")):
            (m / f"Q1.v{k}.md").write_text(body)
            (m / f"Q1.v{k}.claims.json").write_text(json.dumps(claims))
        s = cx.run_mode(self.run, "Q1", revise=False)
        self.assertEqual((s["version"], s["pass"], s["rejected_claims_excluded"]), (2, True, 1))
        (m / "Q1.v3.md").write_text("## 1. Summary\nThe method works well here [pa-c002].\n")
        (m / "Q1.v3.claims.json").write_text(json.dumps(claims))
        s = cx.run_mode(self.run, "Q1", revise=False)               # cites a claim the user rejected
        self.assertEqual((s["version"], s["pass"]), (3, False))
        self.assertTrue((self.run / "summary_y2_Q1.json").exists())


class Redraft(Base):
    """Y1b (synthesize.py --after-u3): withdrawn and confirmed claims, the gate, the X1 guard."""

    def setUp(self):
        super().setUp()
        cl = self.lib / "papers" / "pa" / "claims"
        (cl / "pa-c004.v1.json").write_text(json.dumps(rec("pa-c004", "confirmed by the user")))
        eq = rec("pa-c005", "an equation")
        eq["claim"]["equation"] = {"latex": "x", "number": "1"}
        (cl / "pa-c005.v1.json").write_text(json.dumps(eq))
        self.m = self.run / "untrusted" / "memos"
        self.m.mkdir(parents=True)

    def draft(self, body, given):
        (self.m / "Q1.v1.md").write_text(body)
        (self.m / "Q1.v1.claims.json").write_text(json.dumps([{"id": c} for c in given]))

    def redraft(self, new_body):
        calls = []
        with mock.patch.object(sy.ac, "prompt", side_effect=lambda name, **kw: calls.append(kw) or "order"), \
                mock.patch.object(sy.ac, "call", return_value=({"memo": new_body}, {"seconds": 1, "usage": {}})), \
                mock.patch.object(sy.rc, "questions", return_value="Q1 ..."), mock.patch("builtins.print"):
            code = sy.main(["--run", str(self.run), "--after-u3"])
        return code, json.loads((self.run / "summary_redraft_Q1.json").read_text()), calls

    def test_nothing_changed_makes_no_call(self):
        self.draft("## 1. Summary\nA [pa-c001].\n", ["pa-c001"])
        code, s, calls = self.redraft("unused")
        self.assertEqual((code, s["redrafted"], s["withdrawn"], s["confirmed"], calls), (0, False, 0, 0, []))
        self.assertFalse((self.m / "Q1.v2.md").exists())

    def test_withdrawn_and_confirmed_give_the_next_version(self):
        (self.run / "u3_decisions.json").write_text(json.dumps({"pa-c002": {"verdict": "wrong"},
                                                                "pa-c004": {"verdict": "right"}}))
        self.draft("## 1. Summary\nA [pa-c001]. B [pa-c002].\n", ["pa-c001", "pa-c002", "pa-c005"])
        code, s, calls = self.redraft("## 1. Summary\nA [pa-c001]. C [pa-c004].\n")
        self.assertEqual((code, s["redrafted"], s["version"], s["withdrawn"], s["confirmed"]), (0, True, 2, 1, 1))
        self.assertIn("withdrawn [pa-c002]", calls[0]["failures"])
        self.assertIn("confirmed [pa-c004]", calls[0]["failures"])
        pool = [c["id"] for c in json.loads((self.m / "Q1.v2.claims.json").read_text())]
        self.assertEqual(pool, ["pa-c001", "pa-c004"])          # cited, minus withdrawn, plus confirmed; not c005
        self.assertFalse((self.run / "gates" / "Y1b.json").exists())

    def test_unreviewed_model_option_claim_opens_the_gate(self):
        (self.run / "u3_decisions.json").write_text(json.dumps({"pa-c004": {"verdict": "right"}}))
        self.draft("## 1. Summary\nA [pa-c001].\n", ["pa-c001"])
        code, s, _ = self.redraft("## 1. Summary\nC [pa-c004].\n## 3. Model options\n**Option A.** Use A [pa-c001].\n")
        self.assertEqual((code, s["u3_waiting"]), (0, 1))
        self.assertEqual(json.loads((self.run / "gates" / "Y1b.json").read_text())["status"], "waiting")

    def test_an_equation_claim_x1_did_not_read_fails(self):
        (self.run / "u3_decisions.json").write_text(json.dumps({"pa-c004": {"verdict": "right"}}))
        self.draft("## 1. Summary\nA [pa-c001]. E [pa-c005].\n", ["pa-c001", "pa-c005"])
        code, s, _ = self.redraft("## 1. Summary\nE [pa-c005]. C [pa-c004].\n")
        self.assertEqual((code, s["x1_unread"]), (1, 1))
        self.assertIn("failed", s)


if __name__ == "__main__":
    unittest.main()
