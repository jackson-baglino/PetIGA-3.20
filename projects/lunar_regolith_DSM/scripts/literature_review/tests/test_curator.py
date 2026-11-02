"""Tests for curator.py (L1): publish only checked claims, write once, version on change, refuse without A2."""
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import curator  # noqa: E402
import run_b1  # noqa: E402


def final_rec(cid, status, ok=True, text="supported"):
    return {"claim_id": cid, "status": status, "claim": {"id": cid, "text": "t", "version": 1},
            "verdict": {"claim_id": cid, "text": text}, "e3": {"ok": ok}}


class Publish(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        lit = root / "research" / "Literature"
        lit.mkdir(parents=True)
        (lit / "pa.pdf").write_bytes(b"pdf-a")
        h = hashlib.sha256(b"pdf-a").hexdigest()
        self.contract = root / "c.md"
        self.contract.write_text("## Papers\n\n| ID | DOI | File | Pages | sha | Status | Role |\n|---|---|---|---|---|---|---|\n"
                                 f"| pa | 10.1/a | `Literature/pa.pdf` | 1 | {h[:12]} | ok | core |\n\n## Next\n")
        self.run = root / "run"
        u = self.run / "untrusted" / "pa"
        u.mkdir(parents=True)
        (self.run / "untrusted" / "content" / "pa" / h[:16] / "images").mkdir(parents=True)
        (u / "checklist.json").write_text(json.dumps({"items": []}))
        (u / "read.json").write_text(json.dumps({"symbols": [{"symbol": "x"}]}))
        self.final = [final_rec("pa-c001", "confirmed"), final_rec("pa-c002", "revised"),
                      final_rec("pa-c003", "disputed", text="overstated"), final_rec("pa-c004", "withdrawn"),
                      final_rec("pa-c005", "confirmed", ok=False)]
        (u / "final.json").write_text(json.dumps(self.final))
        self.copy = {"doi": "10.1/a", "kind": "arxiv", "pages": 1, "sha256": h, "gate_opened": True}
        (self.run / "summary_copy_pa.json").write_text(json.dumps(self.copy))
        self.lib = root / "library"
        self.p = [mock.patch.object(run_b1, "CONTRACT", self.contract),
                  mock.patch.object(run_b1, "RESEARCH", root / "research")]
        for x in self.p:
            x.start()

    def tearDown(self):
        for x in self.p:
            x.stop()
        self.tmp.cleanup()

    def approve_copy(self):
        (self.run / "gates").mkdir(exist_ok=True)
        (self.run / "gates" / "A2_pa.json").write_text(json.dumps({"status": "approved", "decision": "read it"}))

    def test_refuses_without_a2_decision(self):
        with self.assertRaises(SystemExit):
            curator.publish(self.run, "pa", self.lib)

    def test_publishes_only_checked_claims_once(self):
        self.approve_copy()
        c = curator.publish(self.run, "pa", self.lib)
        self.assertEqual((c["published"], c["refused_disputed"], c["refused_withdrawn"], c["refused_unchecked"]),
                         (2, 1, 1, 1))
        claims = sorted(p.name for p in (self.lib / "papers" / "pa" / "claims").glob("*.json"))
        self.assertEqual(claims, ["pa-c001.v1.json", "pa-c002.v1.json"])
        card = json.loads((self.lib / "papers" / "pa" / "card.v1.json").read_text())
        self.assertEqual(card["copies"][0]["a2_decision"], "read it")
        self.assertTrue((self.lib / "papers" / "pa" / "analysis" / "run" / "applicability.json").exists())
        again = curator.publish(self.run, "pa", self.lib)                 # same content: nothing new
        self.assertEqual((again["published"], again["same"]), (0, 2))
        self.final[0]["claim"]["text"] = "corrected"
        (self.run / "untrusted" / "pa" / "final.json").write_text(json.dumps(self.final))
        third = curator.publish(self.run, "pa", self.lib)
        self.assertEqual(third["next_version"], 1)
        self.assertTrue((self.lib / "papers" / "pa" / "claims" / "pa-c001.v2.json").exists())
        events = (self.lib / "events.jsonl").read_text().splitlines()
        self.assertEqual(len(events), 3)


if __name__ == "__main__":
    unittest.main()
