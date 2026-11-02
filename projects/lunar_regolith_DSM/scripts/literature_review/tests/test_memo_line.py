"""Y1 reads the library (latest version, minus user rejections); Y2's run mode checks the latest memo."""
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


class Memo(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
