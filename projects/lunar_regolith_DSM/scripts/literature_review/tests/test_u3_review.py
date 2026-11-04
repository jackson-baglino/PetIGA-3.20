"""u3_review.py: which claims the user reviews, and it never prints paper text inside Claude Code."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import u3_review as u3  # noqa: E402

MEMO = """## 1. Summary
A method works [a-c001].
## 3. Model options
**Option A.** Use the sum constraint [a-c002] with the wall term [a-c003] and [a-c002].
## 4. Library silent
Nothing on X [a-c004].
"""


class DisputeReasons(unittest.TestCase):
    OK = {"text": "supported", "quote_on_page": True, "equation": "none", "parameter": "none", "dropped_conditions": []}

    def test_each_failed_test_is_named(self):
        f = {"verdict": dict(self.OK, equation="corrected", corrected_latex="x"),
             "e3": {"ok": False, "words": 10, "found_on": [5, 7]}}
        r = u3.dispute_reasons(f)
        self.assertEqual(len(r), 2)
        self.assertIn("equation corrected", r[0])
        self.assertIn("[5, 7]", r[1])

    def test_supported_statement_with_a_dropped_condition(self):
        r = u3.dispute_reasons({"verdict": dict(self.OK, dropped_conditions=["isothermal"]), "e3": {"ok": True}})
        self.assertEqual(r, ["dropped condition: isothermal"])

    def test_nothing_failed_is_said_plainly(self):
        self.assertIn("unexpected", u3.dispute_reasons({"verdict": self.OK, "e3": {"ok": True}})[0])


class Review(unittest.TestCase):
    def test_model_option_claims(self):
        self.assertEqual(u3.model_option_claims(MEMO), ["a-c002", "a-c003"])

    def test_items_and_refusal(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            m = run / "untrusted" / "memos"
            m.mkdir(parents=True)
            (m / "Q1.v1.md").write_text(MEMO)
            (run / "untrusted" / "x1").mkdir()
            (run / "untrusted" / "x1" / "compare.json").write_text(json.dumps(
                [{"claim_id": "a-c009", "agree": False, "reads": {}}, {"claim_id": "a-c008", "agree": True}]))
            (run / "review_queue.json").write_text(json.dumps([{"stage": "E5", "kind": "disputed",
                                                                "claim_id": "a-c007"}]))
            claims = [{"id": "a-c002"}, {"id": "a-c003", "equation": {"latex": "x"}}]
            with mock.patch.object(u3.sy, "library_claims", return_value=(claims, {}, {})):
                self.assertEqual(u3.items(run, "Q1"),
                                 [("x1", "a-c009"), ("model_option", "a-c002"), ("disputed", "a-c007")])
                with mock.patch.dict(os.environ, {"CLAUDECODE": "1"}):
                    with self.assertRaises(SystemExit):
                        u3.main(["--run", str(run)])
                    self.assertEqual(u3.main(["--run", str(run), "--list"]), 0)     # counts only: allowed

    def test_redo_keeps_the_earlier_answer_and_reopens_l1b(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            m = run / "untrusted" / "memos"
            m.mkdir(parents=True)
            (m / "Q1.v1.md").write_text(MEMO)
            (run / "untrusted" / "a").mkdir()
            (run / "untrusted" / "a" / "final.json").write_text("[]")
            (run / "review_queue.json").write_text(json.dumps([{"stage": "E5", "kind": "disputed",
                                                                "claim_id": "a-c007"}]))
            (run / "u3_decisions.json").write_text(json.dumps({"a-c007": {"verdict": "right", "kind": "disputed"}}))
            (run / "summary_publish_u3_a.json").write_text("{}")
            env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}
            with mock.patch.object(u3.sy, "library_claims", return_value=([], {}, {})), \
                    mock.patch.dict(os.environ, env, clear=True), \
                    mock.patch("builtins.input", side_effect=["w", "mislabeled", "n"]), \
                    mock.patch("builtins.print"):
                with self.assertRaises(SystemExit):
                    u3.main(["--run", str(run), "--redo", "a-c999"])            # not decided: refused
                u3.main(["--run", str(run), "--redo", "a-c007"])
            d = json.loads((run / "u3_decisions.json").read_text())["a-c007"]
            self.assertEqual((d["verdict"], d["previous"][0]["verdict"]), ("wrong", "right"))
            self.assertFalse((run / "summary_publish_u3_a.json").exists())


if __name__ == "__main__":
    unittest.main()
