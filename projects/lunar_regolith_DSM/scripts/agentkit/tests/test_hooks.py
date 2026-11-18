"""Tests of the two hooks, run as Claude Code runs them: JSON on stdin, the decision in the exit code.

guard.py: allowed calls exit 0; forbidden ones exit 2 with a reason; any error exits 2 (fails safe).
stop_check.py: blocks a stop once while a required file is missing.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

KIT = Path(__file__).resolve().parents[1]
ALLOW = ["python scratch/[A-Za-z0-9_-]+\\.py"]


def run_hook(script, payload, env_extra, raw=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith("AGENTKIT_")}
    env.update(env_extra)
    p = subprocess.run([sys.executable, str(KIT / script)], input=raw if raw is not None else json.dumps(payload),
                       capture_output=True, text=True, env=env)
    return p.returncode, p.stderr


class Guard(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ws = Path(self.tmp.name) / "workspace"
        (self.ws / "scratch").mkdir(parents=True)
        self.log = Path(self.tmp.name) / "guard.jsonl"
        self.env = {"AGENTKIT_WORKSPACE": str(self.ws), "AGENTKIT_BASH_ALLOW": json.dumps(ALLOW),
                    "AGENTKIT_GUARD_LOG": str(self.log)}

    def tearDown(self):
        self.tmp.cleanup()

    def call(self, tool, inp, env=None, raw=None):
        return run_hook("guard.py", {"tool_name": tool, "tool_input": inp, "cwd": str(self.ws)},
                        self.env if env is None else env, raw)

    def allowed(self, tool, inp):
        code, err = self.call(tool, inp)
        self.assertEqual(code, 0, f"{tool} {inp} was blocked: {err}")

    def blocked(self, tool, inp, why=None):
        code, err = self.call(tool, inp)
        self.assertEqual(code, 2, f"{tool} {inp} was allowed")
        if why:
            self.assertIn(why, err)

    def test_files_inside_the_workspace(self):
        self.allowed("Read", {"file_path": "task.json"})
        self.allowed("Write", {"file_path": str(self.ws / "derivation.md"), "content": "x"})
        self.allowed("Edit", {"file_path": "scratch/check.py"})
        self.allowed("Glob", {"pattern": "scratch/*.py"})
        self.allowed("Grep", {"pattern": "rho", "path": str(self.ws)})

    def test_files_outside_the_workspace(self):
        self.blocked("Read", {"file_path": "/etc/hosts"}, "outside your working folder")
        self.blocked("Read", {"file_path": "../grade.py"}, "outside your working folder")
        self.blocked("Read", {"file_path": str(Path.home() / ".ssh" / "id_ed25519")})
        self.blocked("Glob", {"pattern": "../**/*.py"})
        self.blocked("Glob", {"pattern": "/Users/**"})
        self.blocked("Grep", {"pattern": "x", "path": "/etc"})

    def test_symlink_out_of_the_workspace(self):
        (self.ws / "escape").symlink_to("/etc")
        self.blocked("Read", {"file_path": "escape/hosts"}, "outside your working folder")

    def test_provenance_never_writable(self):
        self.blocked("Write", {"file_path": str(self.ws / "PROVENANCE.md"), "content": "x"}, "launcher")
        self.blocked("Edit", {"file_path": "../PROVENANCE.md"}, "launcher")

    def test_bash_allow_list(self):
        self.allowed("Bash", {"command": "python scratch/check_weak_form.py"})
        self.blocked("Bash", {"command": "curl https://example.com"}, "not allowed")
        self.blocked("Bash", {"command": "python scratch/a.py; rm -rf ~"}, "one plain command")
        self.blocked("Bash", {"command": "python scratch/a.py && ls"}, "one plain command")
        self.blocked("Bash", {"command": "python $(echo scratch/a.py)"}, "one plain command")
        self.blocked("Bash", {"command": "python scratch/a.py > ../x"}, "one plain command")
        self.blocked("Bash", {"command": "python ../tasks/grade.py"}, "not allowed")
        self.blocked("Bash", {"command": "python -c 'import os'"}, "not allowed")

    def test_unknown_tools_blocked(self):
        self.blocked("WebFetch", {"url": "https://example.com"}, "not allowed for this agent")
        self.blocked("Agent", {"prompt": "x"})

    def test_typed_answer_allowed(self):
        self.allowed("StructuredOutput", {"status": "done"})

    def test_fails_safe(self):
        code, err = self.call("Read", {"file_path": "task.json"}, raw="not json")
        self.assertEqual(code, 2)
        self.assertIn("guard error", err)
        code, err = self.call("Read", {"file_path": "task.json"}, env={})     # no policy at all
        self.assertEqual(code, 2)
        bad = dict(self.env, AGENTKIT_BASH_ALLOW="[not json")
        code, _ = self.call("Bash", {"command": "python scratch/a.py"}, env=bad)
        self.assertEqual(code, 2)

    def test_every_decision_logged(self):
        self.call("Read", {"file_path": "task.json"})
        self.call("Read", {"file_path": "/etc/hosts"})
        rows = [json.loads(s) for s in self.log.read_text().splitlines()]
        self.assertEqual([r["allowed"] for r in rows], [True, False])


class StopCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.env = {"AGENTKIT_WORKSPACE": self.tmp.name, "AGENTKIT_REQUIRED": json.dumps(["derivation.md"])}

    def tearDown(self):
        self.tmp.cleanup()

    def test_blocks_once_while_missing(self):
        code, err = run_hook("stop_check.py", {"stop_hook_active": False}, self.env)
        self.assertEqual(code, 2)
        self.assertIn("derivation.md", err)
        code, _ = run_hook("stop_check.py", {"stop_hook_active": True}, self.env)
        self.assertEqual(code, 0)

    def test_empty_file_counts_as_missing(self):
        (Path(self.tmp.name) / "derivation.md").write_text("")
        self.assertEqual(run_hook("stop_check.py", {"stop_hook_active": False}, self.env)[0], 2)

    def test_allows_when_present(self):
        (Path(self.tmp.name) / "derivation.md").write_text("S1\n")
        self.assertEqual(run_hook("stop_check.py", {"stop_hook_active": False}, self.env)[0], 0)


if __name__ == "__main__":
    unittest.main()
