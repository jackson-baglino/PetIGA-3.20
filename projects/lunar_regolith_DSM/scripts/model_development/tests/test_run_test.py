"""Tests for run_test.py's caps, lock and locked list (python runner and a fake project; no solver runs)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_test as rt  # noqa: E402


class RunTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)
        (self.p / "inputs" / "geometry" / "fam").mkdir(parents=True)
        (self.p / "inputs" / "experiment").mkdir(parents=True)
        (self.p / "inputs" / "solver.opts").write_text("-dof 3\n")
        (self.p / "inputs" / "geometry" / "fam" / "g128.opts").write_text("-Nx 128\n-Ny 128\n")
        (self.p / "inputs" / "geometry" / "g1d.opts").write_text("-Nx 2000\n")
        (self.p / "inputs" / "geometry" / "gfile.opts").write_text("-geom_file x.dat\n# DOF_GRID: 300 300\n")
        (self.p / "ok.py").write_text("print('ran')\n")
        (self.p / "sleep.py").write_text("import time\ntime.sleep(30)\n")
        (self.p / "check.py").write_text("import json\nprint(json.dumps({'pass': True, 'value': 1.0}))\n")
        self.lock = mock.patch.object(rt, "LOCK", self.p / ".run_lock")
        self.lock.start()

    def tearDown(self):
        self.lock.stop()
        self.tmp.cleanup()

    def write_tests(self, tests, lock=False):
        path = self.p / "tests.json"
        path.write_text(json.dumps({"tests": tests}))
        if lock:
            (self.p / "tests.json.lock").write_text(json.dumps({"files": rt.lock_hashes(path, tests, self.p)}))
        return path

    def caps(self, **t):
        return rt.check_caps(dict({"runner": "solver", "geometry": "g128"}, **t), self.p)

    def test_dofs_per_rank(self):
        # 3 * 128 * 128 = 49152 DOFs
        self.assertEqual(self.caps(ranks=1)[0], [])
        self.assertEqual(self.caps(ranks=2)[0], [])
        self.assertTrue(self.caps(ranks=4)[0])            # 12288 per rank < 20k
        self.assertEqual(self.caps(geometry="g1d", ranks=1)[0], [])  # 6000 on one rank is allowed
        self.assertTrue(self.caps(flags=["-Nx", "512", "-Ny", "512"], ranks=1)[0])  # 786k on 1 rank
        self.assertEqual(self.caps(flags=["-Nx", "512", "-Ny", "512"], ranks=8)[0], [])  # 98k per rank

    def test_rank_and_wall_limits(self):
        self.assertTrue(self.caps(ranks=9)[0])
        self.assertTrue(self.caps(ranks=1, wall_minutes=31)[0])

    def test_geom_file_dof_grid(self):
        _, info = self.caps(geometry="gfile", ranks=4)
        self.assertEqual(info["total_dofs"], 3 * 300 * 300)

    def test_python_runner_pass_and_record(self):
        tests = self.write_tests([{"id": "P1", "runner": "python", "script": "ok.py", "check": "check.py",
                                   "wall_minutes": 1}])
        code = rt.main(["--tests", str(tests), "--id", "P1", "--project", str(self.p), "--out",
                        str(self.p / "out"), "--draft"])
        self.assertEqual(code, 0)
        rec = json.loads(next((self.p / "out" / "P1").glob("*.json")).read_text())
        self.assertTrue(rec["pass"])
        self.assertFalse((self.p / ".run_lock").exists())

    def test_unlocked_list_refused_without_draft(self):
        tests = self.write_tests([{"id": "P1", "runner": "python", "script": "ok.py"}])
        with self.assertRaises(SystemExit):
            rt.main(["--tests", str(tests), "--id", "P1", "--project", str(self.p), "--out", str(self.p / "o")])

    def test_locked_list_change_refused(self):
        tests = self.write_tests([{"id": "P1", "runner": "python", "script": "ok.py", "check": "check.py"}],
                                 lock=True)
        self.assertEqual(rt.main(["--tests", str(tests), "--id", "P1", "--project", str(self.p), "--out",
                                  str(self.p / "o")]), 0)
        (self.p / "check.py").write_text("import json\nprint(json.dumps({'pass': True}))\n# loosened\n")
        with self.assertRaises(SystemExit):
            rt.main(["--tests", str(tests), "--id", "P1", "--project", str(self.p), "--out", str(self.p / "o")])

    def test_live_lock_refuses_second_run(self):
        (self.p / ".run_lock").mkdir()
        (self.p / ".run_lock" / "pid").write_text(str(__import__("os").getpid()))
        tests = self.write_tests([{"id": "P1", "runner": "python", "script": "ok.py"}])
        with self.assertRaises(SystemExit):
            rt.main(["--tests", str(tests), "--id", "P1", "--project", str(self.p), "--out", str(self.p / "o"),
                     "--draft"])

    def test_stale_lock_is_cleared(self):
        (self.p / ".run_lock").mkdir()
        (self.p / ".run_lock" / "pid").write_text("999999")
        tests = self.write_tests([{"id": "P1", "runner": "python", "script": "ok.py", "check": "check.py"}])
        self.assertEqual(rt.main(["--tests", str(tests), "--id", "P1", "--project", str(self.p), "--out",
                                  str(self.p / "o"), "--draft"]), 0)

    def test_timeout_kills_and_fails(self):
        tests = self.write_tests([{"id": "S1", "runner": "python", "script": "sleep.py", "check": "check.py",
                                   "wall_minutes": 0.02}])
        code = rt.main(["--tests", str(tests), "--id", "S1", "--project", str(self.p), "--out", str(self.p / "o"),
                        "--draft"])
        self.assertEqual(code, 1)
        rec = json.loads(next((self.p / "o" / "S1").glob("*.json")).read_text())
        self.assertTrue(rec["timed_out"])


if __name__ == "__main__":
    unittest.main()
