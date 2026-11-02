"""Tests for agent_call.py on scripts/agentlib.py: definitions load and match MODELS; token usage is kept for every
call, failed ones too, in the ledger; a failed attempt is retried once; runs.db and transcripts are written."""
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import agent_call as ac  # noqa: E402
import agentlib  # noqa: E402

USAGE = {"input_tokens": 10, "cache_creation_input_tokens": 200, "cache_read_input_tokens": 3000, "output_tokens": 40}
B1_AGENTS = ("morty-reader", "beth-checker", "morty-blind-reader", "beth-synthesizer", "beth-skeptic",
             "beth-recall-judge")


def stream(model, result):
    lines = [{"type": "system", "subtype": "init", "model": model}]
    if result is not None:
        lines.append(result)
    return "\n".join(json.dumps(x) for x in lines)


def ok(model, out):
    return stream(model, {"type": "result", "subtype": "success", "is_error": False, "structured_output": out,
                          "usage": USAGE, "num_turns": 3})


def bad(model):
    return stream(model, {"type": "result", "subtype": "error_during_execution", "is_error": True,
                          "result": "usage limit reached", "usage": USAGE, "num_turns": 1})


def fake(*outs):
    return mock.patch.object(agentlib.subprocess, "run", side_effect=[
        subprocess.CompletedProcess([], 0 if '"success"' in o else 1, stdout=o, stderr="") for o in outs])


class Definitions(unittest.TestCase):
    def test_every_b1_agent_loads_with_its_pinned_model(self):
        for name in B1_AGENTS:
            with self.subTest(name=name):
                a = ac.load(name)
                self.assertEqual(a["model"], ac.MODELS[name.split("-")[0]])
                self.assertIn("## Your place in the line", a["system_prompt"])
                self.assertTrue((ac.AGENTS / f"{name}.settings.json").exists(), "no permission file")

    def test_foreman_is_a_session_agent_with_its_guard(self):
        a = ac.load("beth-foreman")
        self.assertTrue(a["session"])
        self.assertIn("foreman_bash_guard.py", json.dumps(a["hooks"]))
        self.assertTrue({"Edit", "Write", "Agent"} <= set(a["disallowed"]))
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import foreman
        d = foreman.definition()["beth-foreman"]
        self.assertEqual(d["model"], ac.MODELS["beth"])
        self.assertIn("hooks", d)
        self.assertEqual(d["skills"], ["lr-status"])

    def test_work_orders_have_no_job_text(self):
        """prompts/ holds only work orders: each starts 'Work order:' and is short."""
        for p in ac.PROMPTS.glob("*.md"):
            with self.subTest(prompt=p.name):
                t = p.read_text()
                self.assertTrue(t.startswith("Work order:"))
                self.assertLess(len(t), 600)


class Ledger(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run = Path(self.tmp.name) / "run"
        self.run.mkdir()
        self.ledger = self.run / "calls.jsonl"
        self.folder = Path(self.tmp.name) / "clean"
        self.folder.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def test_success_retry_and_failure(self):
        m = ac.MODELS["morty"]
        with fake(ok(m, {"items": []})):
            r, meta = ac.call(self.folder, "x", "blind", "morty-blind-reader", ledger=self.ledger, step="X1")
        self.assertEqual(r, {"items": []})
        self.assertEqual((meta["usage"], meta["attempts"]), (USAGE, 1))
        with fake(bad(m), ok(m, {"items": []})):
            r, meta = ac.call(self.folder, "x", "blind", "morty-blind-reader", ledger=self.ledger, step="X1")
        self.assertEqual((r, meta["attempts"]), ({"items": []}, 2))
        with fake(bad(m), "not json"):
            r, meta = ac.call(self.folder, "x", "blind", "morty-blind-reader", ledger=self.ledger, step="X1")
        self.assertIsNone(r)
        self.assertEqual(meta["exit"], 1)
        self.assertTrue(meta["error"])
        calls = ac.ledger_calls(self.ledger)
        self.assertEqual(len(calls), 3)
        tot = ac.usage_total(calls)
        self.assertEqual((tot["calls"], tot["calls_without_usage"]), (3, 1))
        con = sqlite3.connect(self.run / "runs.db")
        rows = con.execute("SELECT status FROM calls ORDER BY id").fetchall()
        self.assertEqual([r[0] for r in rows], ["ok", "failed", "ok", "failed", "failed"])
        self.assertEqual(len(list((self.run / "transcripts" / "X1").glob("*.jsonl"))), 5)

    def test_wrong_model_or_schema_fails(self):
        with fake(ok("claude-opus-5-5", {"items": []}), ok("claude-opus-5-5", {"items": []})):
            r, meta = ac.call(self.folder, "x", "blind", "morty-blind-reader", ledger=self.ledger, step="X1")
        self.assertIsNone(r)
        self.assertIn("session model", meta["error"])
        m = ac.MODELS["morty"]
        with fake(ok(m, {"wrong": 1}), ok(m, {"wrong": 1})):
            r, meta = ac.call(self.folder, "x", "blind", "morty-blind-reader", ledger=self.ledger, step="X1")
        self.assertIsNone(r)

    def test_command_has_system_prompt_tools_and_settings(self):
        cmd = agentlib.command(ac.load("beth-skeptic"), "task", {}, settings="s.json")
        self.assertEqual(cmd[cmd.index("--tools") + 1], "")
        self.assertIn("--system-prompt", cmd)
        self.assertEqual(cmd[cmd.index("--settings") + 1], "s.json")
        self.assertEqual(cmd[cmd.index("--setting-sources") + 1], "")

    def test_no_ledger_writes_nothing(self):
        with fake(ok(ac.MODELS["morty"], {"items": []})):
            ac.call(self.folder, "x", "blind", "morty-blind-reader")
        self.assertFalse(self.ledger.exists())
        self.assertEqual(ac.ledger_calls(self.ledger), [])


if __name__ == "__main__":
    unittest.main()
