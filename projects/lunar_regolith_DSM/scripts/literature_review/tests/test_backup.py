"""backup.py and restore_test.py (BK, T-restore) on a small fake run and library."""
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import backup as bk  # noqa: E402
import restore_test as rt  # noqa: E402
import run_b1  # noqa: E402


class Backup(unittest.TestCase):
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
        (self.lib / "papers" / "pa" / "claims").mkdir(parents=True)
        (self.lib / "papers" / "pa" / "claims" / "pa-c001.v1.json").write_text("{}")
        (self.lib / "papers" / "pa" / "card.v1.json").write_text("{}")
        (self.lib / "events.jsonl").write_text("{}\n")
        self.run = root / "run1"
        m = self.run / "untrusted" / "memos"
        m.mkdir(parents=True)
        (m / "Q1.v1.md").write_text("## 1. Summary\nThe method converges in this setting [pa-c001].\n")
        (m / "Q1.v1.claims.json").write_text(json.dumps([{"id": "pa-c001", "text": "t"}]))
        self.dest = root / "archive"
        self.p = [mock.patch.object(run_b1, "CONTRACT", contract), mock.patch.object(run_b1, "RESEARCH", root / "r")]
        for x in self.p:
            x.start()

    def tearDown(self):
        for x in self.p:
            x.stop()
        self.tmp.cleanup()

    def test_backup_then_restore(self):
        s = bk.backup(self.run, self.lib, self.dest)
        self.assertGreaterEqual(s["files"], 5)
        with self.assertRaises(SystemExit):
            bk.backup(self.run, self.lib, self.dest)                     # written once
        r = rt.restore_test(self.run)
        self.assertTrue(r["restore_ok"], r)
        dest = Path(s["dest"])
        (dest / "library" / "papers" / "pa" / "claims" / "pa-c001.v1.json").write_text('{"x": 1}')
        r = rt.restore_test(self.run)
        self.assertFalse(r["restore_ok"])
        self.assertTrue(any("hash differs" in m for m in r["restore_mismatches"]))


if __name__ == "__main__":
    unittest.main()
