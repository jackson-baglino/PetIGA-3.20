"""Tests of PROVENANCE.md: the chain verifies, and any change after writing is found and located."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import provenance  # noqa: E402


class Chain(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "PROVENANCE.md"
        p = provenance.Provenance(self.path, "PROVENANCE: test", "Written only by the launcher.")
        p.add("run started", ["agent: beth-mathematician", "task: T01"])
        p.add("grade", ["pass: True"])
        self.final = p.add("verdict", ["**PASS**"])

    def tearDown(self):
        self.tmp.cleanup()

    def test_intact(self):
        ok, msg, last = provenance.verify(self.path)
        self.assertTrue(ok, msg)
        self.assertEqual(last, self.final)

    def test_edited_entry_is_located(self):
        self.path.write_text(self.path.read_text().replace("pass: True", "pass: False"))
        ok, msg, _ = provenance.verify(self.path)
        self.assertFalse(ok)
        self.assertIn("entry 2", msg)

    def test_forged_entry_without_chain(self):
        with self.path.open("a") as f:
            f.write("\n## 4. note\n\n- all fine\n")
        self.assertFalse(provenance.verify(self.path)[0])

    def test_rewritten_tail_is_caught_by_the_manifest_hash(self):
        """Someone who rewrites the last entry and recomputes its hash still changes the final hash."""
        text = self.path.read_text()
        head = text[:text.index("## 3. verdict")]
        prev = provenance.CHAIN.search(head.strip().split("\n")[-1]).group(1)
        entry = "## 3. verdict\n\n- **FAIL**\n"
        forged = head + entry + f"\n`chain: {provenance._hash(prev, entry)}`\n"
        self.path.write_text(forged)
        ok, _, last = provenance.verify(self.path)
        self.assertTrue(ok)                       # internally consistent ...
        self.assertNotEqual(last, self.final)     # ... but no longer matches manifest.json


if __name__ == "__main__":
    unittest.main()
