"""Tests for the model-line hooks (.claude/hooks/md_guard.py, md_stop.py), run as the harness runs them."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
GUARD = ROOT / ".claude" / "hooks" / "md_guard.py"
STOP = ROOT / ".claude" / "hooks" / "md_stop.py"


def run(script, payload, env):
    full = dict(os.environ, **env)
    p = subprocess.run([sys.executable, str(script)], input=json.dumps(payload), capture_output=True, text=True,
                       env=full)
    return p.returncode, p.stderr


class GuardTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        t = Path(os.path.realpath(self.tmp.name))
        self.read, self.write, self.deny = t / "src", t / "run" / "derivation", t / "docs" / "planning"
        for d in (self.read, self.write, self.deny):
            d.mkdir(parents=True)
        self.env = {"MD_AGENT": "beth-modeler", "MD_CWD": str(t),
                    "MD_READ_ROOTS": json.dumps([str(self.read), str(self.write)]),
                    "MD_WRITE_ROOTS": json.dumps([str(self.write)]),
                    "MD_DENY_ROOTS": json.dumps([str(self.deny)]),
                    "MD_BASH_ALLOW": json.dumps([r"make( -C \S+)?", r"\S*venv_workflows/bin/python \S+\.py( .*)?"])}

    def tearDown(self):
        self.tmp.cleanup()

    def tool(self, name, **inp):
        return run(GUARD, {"tool_name": name, "tool_input": inp}, self.env)[0]

    def test_no_agent_is_noop(self):
        self.assertEqual(run(GUARD, {"tool_name": "Write", "tool_input": {"file_path": "/etc/x"}}, {"MD_AGENT": ""})[0], 0)

    def test_write_inside_and_outside(self):
        self.assertEqual(self.tool("Write", file_path=str(self.write / "ledger.json")), 0)
        self.assertEqual(self.tool("Write", file_path=str(self.read / "assembly.c")), 2)
        self.assertEqual(self.tool("Edit", file_path=str(self.write / ".." / "checks" / "c.py")), 2)

    def test_read_roots_and_deny(self):
        self.assertEqual(self.tool("Read", file_path=str(self.read / "assembly.c")), 0)
        self.assertEqual(self.tool("Read", file_path=str(self.deny / "18_f2.md")), 2)
        self.assertEqual(self.tool("Read", file_path="/etc/passwd"), 2)
        self.assertEqual(self.tool("Grep", path=str(self.deny)), 2)

    def test_relative_path_resolves_against_cwd(self):
        self.assertEqual(self.tool("Read", file_path="docs/planning/18_f2.md"), 2)
        self.assertEqual(self.tool("Read", file_path="src/assembly.c"), 0)

    def test_bash_allowlist_and_chaining(self):
        self.assertEqual(self.tool("Bash", command="make -C build"), 0)
        self.assertEqual(self.tool("Bash", command="venv_workflows/bin/python checks/m3_a.py --x 1"), 0)
        self.assertEqual(self.tool("Bash", command="make && rm -rf /"), 2)
        self.assertEqual(self.tool("Bash", command="curl http://example.com"), 2)
        self.assertEqual(self.tool("Bash", command="make $(echo x)"), 2)


class StopTest(unittest.TestCase):
    def test_blocks_twice_then_allows(self):
        with tempfile.TemporaryDirectory() as state:
            env = {"MD_AGENT": "beth-modeler", "MD_STOP_CHECK": "echo missing ledger step S7; exit 1",
                   "MD_HOOK_STATE": state}
            codes = [run(STOP, {}, env) for _ in range(3)]
            self.assertEqual([c for c, _ in codes], [2, 2, 0])
            self.assertIn("missing ledger step S7", codes[0][1])

    def test_passing_check_allows(self):
        with tempfile.TemporaryDirectory() as state:
            env = {"MD_AGENT": "beth-modeler", "MD_STOP_CHECK": "exit 0", "MD_HOOK_STATE": state}
            self.assertEqual(run(STOP, {}, env)[0], 0)


if __name__ == "__main__":
    unittest.main()
