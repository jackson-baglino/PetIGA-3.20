"""Y3: one Y1 revision when statements are not supported; the revision must pass Y2; the rest to the queue."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import check_citations  # noqa: E402
import skeptic  # noqa: E402


def out(*verdicts):
    return {"support": {"statements": [{"verdict": v, "section": "2", "summary": "s", "reason": "r"}
                                       for v in verdicts]},
            "omissions": {"omissions": [{"kind": "gap", "importance": "high", "what": "w", "claim_ids": []}]}}


class Y3Revision(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run = Path(self.tmp.name)
        m = self.run / "untrusted" / "memos"
        m.mkdir(parents=True)
        (m / "Q1.v2.md").write_text("memo")

    def tearDown(self):
        self.tmp.cleanup()

    def fake_y3(self, results):
        def f(run, question, version, source="library"):
            o = results.pop(0)
            return o, {"version": version}
        return f

    def test_revises_once_then_queues_the_rest(self):
        def fake_synth(argv):
            v = argv[argv.index("--version") + 1]
            (self.run / "untrusted" / "memos" / f"Q1.v{v}.md").write_text("revised")
            self.assertIn("--skeptic", argv)
        with mock.patch.object(skeptic, "y3", self.fake_y3([out("supported", "unsupported"),
                                                            out("supported", "overstated")])), \
                mock.patch.object(skeptic.sy, "main", side_effect=fake_synth) as synth, \
                mock.patch.object(check_citations, "run_mode", return_value={"pass": True}):
            skeptic.main(["--run", str(self.run)])
        self.assertEqual(synth.call_count, 1)                                  # once, even if issues remain
        s = json.loads((self.run / "summary_y3_Q1.json").read_text())
        self.assertEqual((s["version"], s["revised_after_y3"], s["first_check"]["not_supported"]), (3, True, 1))
        q = json.loads((self.run / "review_queue.json").read_text())
        self.assertEqual(sorted(x["kind"] for x in q), ["gap", "overstated"])
        self.assertTrue(all(x["item"] == "Q1.v3" for x in q))

    def test_no_revision_when_all_supported(self):
        with mock.patch.object(skeptic, "y3", self.fake_y3([out("supported")])), \
                mock.patch.object(skeptic.sy, "main") as synth:
            skeptic.main(["--run", str(self.run)])
        synth.assert_not_called()

    def test_revision_failing_y2_stops(self):
        with mock.patch.object(skeptic, "y3", self.fake_y3([out("unsupported")])), \
                mock.patch.object(skeptic.sy, "main"), \
                mock.patch.object(check_citations, "run_mode", return_value={"pass": False, "failing_sentences": {}}):
            with self.assertRaises(SystemExit):
                skeptic.main(["--run", str(self.run)])


if __name__ == "__main__":
    unittest.main()
