"""Tests for run_b1.py, the B1 line, with stub stations (no agent is called)."""
import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import run_b1  # noqa: E402
import wf_launch  # noqa: E402

STUBS = {
    "s_ok.py": "import sys, json, pathlib\na = sys.argv\nrun = pathlib.Path(a[a.index('--run') + 1])\n"
               "name = a[a.index('--out') + 1]\n(run / name).write_text(json.dumps({'n': 1}))\n",
    "s_fail.py": "import sys\nsys.exit(2)\n",
    "s_quiet.py": "pass\n",
    "s_gate.py": "import sys, json, pathlib\na = sys.argv\nrun = pathlib.Path(a[a.index('--run') + 1])\n"
                 "(run / 'c.json').write_text('{}')\n(run / 'gates').mkdir(exist_ok=True)\n"
                 "(run / 'gates' / 'C.json').write_text(json.dumps({'status': 'waiting'}))\n",
}


def spec(nodes):
    return {"title": "t", "summary": "", "status": "", "lanes": [], "edges": [],
            "nodes": [{"id": i, "title": i, "lane": "script", "kind": "script", "run": r} for i, r in nodes]}


def ok(out, per_paper=False):
    return [["s_ok.py", "--run", "{run}", "--out", out]]


class Line(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.here = root / "stations"
        self.here.mkdir()
        for n, code in STUBS.items():
            (self.here / n).write_text(code)
        research = root / "research"
        (research / "Literature").mkdir(parents=True)
        rows = []
        for pid in ("pa", "pb"):
            f = research / "Literature" / f"{pid}.pdf"
            f.write_bytes(pid.encode())
            h = hashlib.sha256(pid.encode()).hexdigest()[:12]
            rows.append(f"| {pid} | 10.1/x | `Literature/{pid}.pdf` | 1 | {h} | on disk | core |")
        self.contract = root / "B1_CONTRACT.md"
        self.contract.write_text("# c\n\n## Papers\n\n| ID | DOI | File | Pages | sha256 | Status | Role |\n"
                                 "|---|---|---|---|---|---|---|\n" + "\n".join(rows) + "\n\n## Outputs\n")
        self.perm = root / "perm"
        (self.perm / "profiles").mkdir(parents=True)
        (self.perm / "launcher.allow").write_text("s_ok.py offline\ns_fail.py offline\ns_quiet.py offline\n"
                                                  "s_gate.py offline\n")
        real = Path(wf_launch.PROJECT) / ".claude/permissions/literature-review/profiles/offline.sb"
        (self.perm / "profiles" / "offline.sb").write_text(real.read_text())
        self.spec = root / "flowchart.json"
        self.run = root / "run"
        self.patches = [mock.patch.object(run_b1, "HERE", self.here), mock.patch.object(run_b1, "SPEC", self.spec),
                        mock.patch.object(run_b1, "CONTRACT", self.contract),
                        mock.patch.object(run_b1, "RESEARCH", research),
                        mock.patch.object(run_b1, "PYTHON", Path(sys.executable)),
                        mock.patch.object(run_b1, "FREEZE", root / "freeze"),
                        mock.patch.object(wf_launch, "PERM", self.perm)]
        for p in self.patches:
            p.start()
        stages_ = [
            ("U", {"per": "once", "requires": [], "station": None, "gate": "blocking", "ask": "go?"}),
            ("A", {"per": "paper", "requires": ["U"], "station": "s_ok.py", "gate": "none",
                   "steps": ok("a_{paper}.json"), "done": [{"path": "a_{paper}.json", "has": ["n"]}]}),
            ("B", {"per": "once", "requires": ["A"], "station": "s_ok.py", "gate": "none",
                   "steps": ok("b.json"), "done": [{"path": "b.json"}]}),
            ("C", {"per": "once", "requires": ["B"], "station": "s_gate.py", "gate": "blocking",
                   "steps": [["s_gate.py", "--run", "{run}"]], "done": [{"path": "c.json"}]}),
            ("D", {"per": "once", "requires": ["C"], "station": "missing.py", "gate": "none",
                   "steps": [["missing.py"]], "done": [{"path": "d.json"}]}),
        ]
        self.spec.write_text(json.dumps(spec(stages_)))

    def tearDown(self):
        for p in self.patches:
            p.stop()
        self.tmp.cleanup()

    def cli(self, *args):
        return run_b1.main(list(args))

    def test_walks_in_order_stops_at_gates_and_unbuilt(self):
        run = str(self.run)
        self.cli("init", "--run", run, "--dev")
        self.assertEqual(self.cli("run", "--run", run), run_b1.WAITING)          # user gate U
        self.assertFalse((self.run / "a_pa.json").exists())
        self.cli("approve", "--run", run, "U", "--decision", "yes")
        self.assertEqual(self.cli("run", "--run", run), run_b1.WAITING)          # C opened its gate
        for f in ("a_pa.json", "a_pb.json", "b.json", "c.json"):
            self.assertTrue((self.run / f).exists(), f)
        self.cli("approve", "--run", run, "C", "--decision", "accepted")
        self.assertEqual(self.cli("run", "--run", run), run_b1.NOT_BUILT)        # D's station is missing
        rows = {r["stage"]: r for r in run_b1.status_rows(run)}
        self.assertEqual((rows["A"]["done"], rows["A"]["total"]), (2, 2))
        self.assertFalse(rows["D"]["built"])

    def test_requires_are_enforced(self):
        run = str(self.run)
        self.cli("init", "--run", run, "--dev")
        self.cli("approve", "--run", run, "U", "--decision", "yes")
        plist = run_b1.papers()
        s, i = run_b1.next_item(self.run, plist, None)
        self.assertEqual((s["id"], i), ("A", "pa"))
        self.assertEqual(self.cli("run", "--run", run, "--max", "1"), 0)
        s, i = run_b1.next_item(self.run, plist, None)
        self.assertEqual((s["id"], i), ("A", "pb"))                              # B waits for every paper

    def test_failure_and_missing_output_stop_the_line(self):
        d = json.loads(self.spec.read_text())
        d["nodes"][1]["run"]["steps"] = [["s_fail.py"]]
        self.spec.write_text(json.dumps(d))
        run = str(self.run)
        self.cli("init", "--run", run, "--dev")
        self.cli("approve", "--run", run, "U", "--decision", "yes")
        self.assertEqual(self.cli("run", "--run", run), run_b1.FAILED)
        d["nodes"][1]["run"]["steps"] = [["s_quiet.py"]]                         # exit 0, no output
        self.spec.write_text(json.dumps(d))
        self.assertEqual(self.cli("run", "--run", run), run_b1.FAILED)

    def test_station_not_on_the_allow_list_is_refused(self):
        with self.assertRaises(wf_launch.LaunchError):
            wf_launch.command("missing.py", [], self.run)
        cmd = wf_launch.command("s_ok.py", ["--run", "x"], self.run)
        self.assertEqual(cmd[:2], ["sandbox-exec", "-f"])
        self.assertTrue(cmd[2].endswith("offline.sb"))

    def test_sandbox_blocks_writes_outside_the_run(self):
        d = json.loads(self.spec.read_text())
        out_dir = Path(tempfile.mkdtemp(dir="/private/tmp"))              # outside the run and the temp folder
        self.addCleanup(lambda: (out_dir / "outside.json").unlink(missing_ok=True) or out_dir.rmdir())
        outside = out_dir / "outside.json"
        d["nodes"][1]["run"]["steps"] = [["s_ok.py", "--run", str(outside.parent), "--out", "outside.json"]]
        self.spec.write_text(json.dumps(d))
        run = str(self.run)
        self.cli("init", "--run", run, "--dev")
        self.cli("approve", "--run", run, "U", "--decision", "yes")
        self.assertEqual(self.cli("run", "--run", run), run_b1.FAILED)
        self.assertFalse(outside.exists())

    def test_dev_run_on_one_paper_with_a_skipped_stage(self):
        run = str(self.run)
        self.cli("init", "--run", run, "--dev", "--papers", "pb", "--skip", "C")
        self.cli("approve", "--run", run, "U", "--decision", "yes")
        self.assertEqual(self.cli("run", "--run", run), run_b1.NOT_BUILT)        # past C (skipped) to D
        self.assertTrue((self.run / "a_pb.json").exists())
        self.assertFalse((self.run / "a_pa.json").exists())
        self.assertFalse((self.run / "c.json").exists())
        rows = {r["stage"]: r for r in run_b1.status_rows(run)}
        self.assertEqual((rows["A"]["total"], rows["C"]["skipped"]), (1, True))

    def test_papers_and_skip_only_in_dev(self):
        with self.assertRaises(SystemExit):
            self.cli("init", "--run", str(self.run), "--papers", "pa")

    def test_execute_run_needs_a_freeze(self):
        with self.assertRaises(SystemExit):
            self.cli("init", "--run", str(self.run))

    def test_contract_hash_mismatch_fails(self):
        (Path(self.tmp.name) / "research" / "Literature" / "pa.pdf").write_bytes(b"changed")
        with self.assertRaises(SystemExit):
            run_b1.papers()

    def test_approve_only_blocking_gates(self):
        self.cli("init", "--run", str(self.run), "--dev")
        with self.assertRaises(SystemExit):
            self.cli("approve", "--run", str(self.run), "B", "--decision", "x")


class Freeze(unittest.TestCase):
    def test_phase_and_diff(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(run_b1, "FREEZE", Path(tmp)), \
                mock.patch.object(run_b1, "snapshot", return_value={"a": "1", "b": "2"}):
            self.assertEqual(run_b1.phase(), ("develop", None))
            (Path(tmp) / "b1-v1.json").write_text(json.dumps({"version": "b1-v1", "files": {"a": "1", "b": "2"}}))
            self.assertEqual(run_b1.phase(), ("execute", "b1-v1"))
            self.assertEqual(run_b1.freeze_diff(), [])
            (Path(tmp) / "b1-v2.json").write_text(json.dumps({"version": "b1-v2", "files": {"a": "1", "b": "3"},
                                                              "thawed": {"reason": "x"}}))
            self.assertEqual(run_b1.phase(), ("develop", "b1-v2"))
            self.assertEqual(run_b1.freeze_diff(), ["b"])


class ParserInFreeze(unittest.TestCase):
    def test_parser_state_in_snapshot_and_missing_env_is_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            ref = Path(tmp) / "main"
            ref.write_text("abc123\n")
            with mock.patch.object(run_b1, "DOCPARSE_PY", Path(tmp) / "no_python"), \
                    mock.patch.object(run_b1, "PARSER_MODEL", ref):
                st = run_b1.parser_state()
                self.assertEqual(st, {"venv_docparse#installed": "missing", "parser_model#snapshot": "abc123"})
                self.assertEqual(run_b1.snapshot()["parser_model#snapshot"], "abc123")
                ref.write_text("def456\n")                  # a new model snapshot changes the freeze record
                self.assertEqual(run_b1.snapshot()["parser_model#snapshot"], "def456")
        self.assertIn("requirements-docparse.lock", run_b1.FROZEN)


class RealGraph(unittest.TestCase):
    def test_b1_graph_order_matches_the_contract(self):
        """In the real flowchart: calibration before reading; V1 (X1, U3a) before Y2; U4 before the backup."""
        by = {s["id"]: s for s in run_b1.stages()}
        self.assertIn("C0", by["E1"]["requires"])
        self.assertEqual(set(by["Y2"]["requires"]), {"X1", "U3a"})
        self.assertIn("X1", by["U3a"]["requires"])      # the user's U3 pass includes X1's disagreements
        c0 = [st[0] + (" --baseline" if "--baseline" in st else "") for st in by["C0"]["steps"]]
        self.assertLess(c0.index("baseline.py"), c0.index("calibrate.py"))   # the baseline is scored even
        self.assertLess(c0.index("calibrate.py --baseline"), c0.index("calibrate.py"))  # if the checker fails
        self.assertEqual(by["Y1"]["requires"], ["L1"])
        self.assertEqual(by["BK"]["requires"], ["U4"])
        for gate in ("U0", "C0", "A2", "X1", "U3a", "U4"):
            self.assertEqual(by[gate]["gate"], "blocking", gate)


if __name__ == "__main__":
    unittest.main()
