"""Tests for agent_call.py: token usage is kept for every call, failed ones too, in the ledger."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import agent_call as ac  # noqa: E402

USAGE = {"input_tokens": 10, "cache_creation_input_tokens": 200, "cache_read_input_tokens": 3000, "output_tokens": 40}


def fake(code, stdout, stderr=""):
    return mock.patch.object(ac.subprocess, "run",
                             return_value=subprocess.CompletedProcess([], code, stdout=stdout, stderr=stderr))


class Ledger(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ledger = Path(self.tmp.name) / "calls.jsonl"

    def tearDown(self):
        self.tmp.cleanup()

    def test_success_failure_and_no_json(self):
        ok = json.dumps({"structured_output": {"a": 1}, "usage": USAGE, "num_turns": 3, "is_error": False})
        bad = json.dumps({"result": "usage limit reached", "usage": USAGE, "is_error": True})
        with fake(0, ok):
            r, m = ac.call(self.tmp.name, "x", "read", "morty-reader", ledger=self.ledger, step="E1")
        self.assertEqual(r, {"a": 1})
        self.assertEqual(m["usage"], USAGE)
        with fake(1, bad):
            r, m = ac.call(self.tmp.name, "x", "read", "beth-checker", ledger=self.ledger, step="E4.0")
        self.assertIsNone(r)
        self.assertEqual(m["usage"], USAGE)
        self.assertIn("usage limit", m["error"])
        with fake(1, "", "crash"):
            r, m = ac.call(self.tmp.name, "x", "read", "beth-checker", ledger=self.ledger, step="E4.0")
        self.assertIsNone(r)
        self.assertIsNone(m["usage"])
        calls = ac.ledger_calls(self.ledger)
        self.assertEqual([c["step"] for c in calls], ["E1", "E4.0", "E4.0"])
        tot = ac.usage_total(calls)
        self.assertEqual(tot["cache_read_input_tokens"], 6000)
        self.assertEqual(tot["output_tokens"], 80)
        self.assertEqual((tot["calls"], tot["calls_without_usage"]), (3, 1))

    def test_no_ledger_writes_nothing(self):
        with fake(0, json.dumps({"structured_output": {}, "usage": USAGE})):
            ac.call(self.tmp.name, "x", "read", "morty-reader")
        self.assertFalse(self.ledger.exists())
        self.assertEqual(ac.ledger_calls(self.ledger), [])


if __name__ == "__main__":
    unittest.main()
