"""blind_read.py (X1): a reused read is counted, and calls go to the run's ledger."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import blind_read as br  # noqa: E402


class Reuse(unittest.TestCase):
    def test_reused_read_is_counted_and_new_call_uses_the_run_ledger(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            x1 = run / "untrusted" / "x1"
            x1.mkdir(parents=True)
            (x1 / "pa.json").write_text(json.dumps({"memo_version": 2, "reads": []}))   # made earlier, same memo
            (x1 / "pb.json").write_text(json.dumps({"memo_version": 1, "reads": []}))   # an older memo: redo
            for p in ("pa", "pb"):
                (run / "untrusted" / p / "pages").mkdir(parents=True)
            items = [("pa", {"id": "pa-c001", "page_pdf": 1, "parameter": {"quantity": "x"}}),
                     ("pb", {"id": "pb-c001", "page_pdf": 1, "parameter": {"quantity": "y"}})]
            calls = []

            def fake_call(folder, text, schema, agent, ledger=None, **tag):
                calls.append(Path(ledger))
                return {"reads": []}, {}
            with mock.patch.object(br, "cited_items", return_value=(2, items)), \
                    mock.patch.object(br.ac, "call", side_effect=fake_call), \
                    mock.patch.object(br.ac, "prompt", return_value="task"):
                br.main(["--run", str(run)])
            s = json.loads((run / "summary_x1_reads_Q1.json").read_text())
            self.assertEqual((s["calls"], s["reused"], s["reused_papers"]), (1, 1, ["pa"]))
            self.assertEqual(calls, [run.resolve() / "calls.jsonl"])


if __name__ == "__main__":
    unittest.main()
