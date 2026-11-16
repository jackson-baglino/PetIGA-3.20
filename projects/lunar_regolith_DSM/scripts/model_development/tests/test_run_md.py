"""Tests for run_md.py's control flow without agents: usage-limit stop, gates, approve, reopen."""
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_md  # noqa: E402


class RunMdTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.runs = Path(self.tmp.name)
        self.p = mock.patch.object(run_md, "RUNS", self.runs)
        self.p.start()
        self.n = mock.patch.object(run_md, "notify", lambda *a: None)
        self.n.start()
        rdir = self.runs / "t1"
        for sub in ("derivation", "checks", "gates", "triage", "tests", "transcripts"):
            (rdir / sub).mkdir(parents=True)
        st = {"run_id": "t1", "problem": "toy", "dev": True,
              "stages": {s: {"status": "pending"} for s in run_md.STAGES}, "feedback": {}, "model_name": None,
              "worktree": None, "branch": None, "layers_done": [], "triage_n": 0}
        st["stages"]["M8"] = {"status": "skipped"}
        run_md.save_state(rdir, st)
        (rdir / "decisions.json").write_text("[]")
        self.rdir = rdir

    def tearDown(self):
        self.p.stop()
        self.n.stop()
        self.tmp.cleanup()

    def run_cli(self, *argv):
        with redirect_stdout(io.StringIO()) as buf:
            code = run_md.main(list(argv))
        return code, buf.getvalue()

    def test_usage_limit_stops_cleanly(self):
        def boom(*a, **k):
            raise run_md.UsageLimit("You've hit your 5-hour limit")
        with mock.patch.dict(run_md.STATIONS, {"M1": boom}):
            code, out = self.run_cli("run", "--run", "t1")
        self.assertEqual(code, 3)
        self.assertIn("Resume later", out)
        st = run_md.load_state(self.rdir)
        self.assertEqual(st["stages"]["M1"]["status"], "pending")

    def test_gate_stops_and_approve_moves_on(self):
        def fake_m1(rdir, st):
            st["stages"]["M1"] = {"status": "gate_waiting"}
            run_md.save_state(rdir, st)
            return False
        calls = []
        with mock.patch.dict(run_md.STATIONS, {"M1": fake_m1, "M2": lambda r, s: calls.append("M2") or False}):
            self.run_cli("run", "--run", "t1")
            self.assertEqual(calls, [])                       # stopped at G:Scope
            code, _ = self.run_cli("approve", "G:Scope", "--run", "t1", "--decision", "ok, insulating grains")
            self.run_cli("run", "--run", "t1")
            self.assertEqual(calls, ["M2"])
        dec = json.loads((self.rdir / "decisions.json").read_text())
        self.assertEqual(dec[0]["gate"], "G:Scope")

    def test_selection_needs_model_name(self):
        st = run_md.load_state(self.rdir)
        st["stages"]["M2"] = {"status": "gate_waiting"}
        run_md.save_state(self.rdir, st)
        with self.assertRaises(SystemExit):
            self.run_cli("approve", "G:Selection", "--run", "t1", "--decision", "DDM")
        self.run_cli("approve", "G:Selection", "--run", "t1", "--decision", "DDM", "--model-name", "ddm")
        self.assertEqual(run_md.load_state(self.rdir)["model_name"], "ddm")

    def test_reopen_resets_later_stages(self):
        st = run_md.load_state(self.rdir)
        for s in ("M1", "M2", "M3"):
            st["stages"][s] = {"status": "approved"}
        run_md.save_state(self.rdir, st)
        self.run_cli("reopen", "M2", "--run", "t1", "--note", "add a level-set candidate")
        st = run_md.load_state(self.rdir)
        self.assertEqual(st["stages"]["M1"]["status"], "approved")
        self.assertEqual(st["stages"]["M2"]["status"], "pending")
        self.assertEqual(st["stages"]["M3"]["status"], "pending")
        self.assertEqual(st["stages"]["M8"]["status"], "skipped")
        self.assertIn("level-set", run_md.feedback_text(st, "M2"))

    def test_question_stop_and_answer(self):
        def asks(rdir, st):
            run_md.ask(rdir, st, "M1", {"questions": ["Insulating or conducting grains?"], "summary": "s"})
            return False
        with mock.patch.dict(run_md.STATIONS, {"M1": asks}):
            self.run_cli("run", "--run", "t1")
        self.assertIn("question at M1", run_md.status_line("t1"))
        self.run_cli("approve", "M1", "--run", "t1", "--decision", "insulating")
        self.assertEqual(run_md.load_state(self.rdir)["stages"]["M1"]["status"], "pending")


if __name__ == "__main__":
    unittest.main()
