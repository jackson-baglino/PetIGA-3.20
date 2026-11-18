"""Tests of the launcher with the fake `claude` (no model is called), on beth-mathematician and task T01.

Each scenario checks the verdict and the run record: PASS, FAIL, ABSTAIN, ERROR; the hooks really run through
the composed settings; the audit catches a barrier that failed; PROVENANCE.md verifies against manifest.json;
and a real model call is refused while any enforcer is unbuilt or no red-team PASS matches the current barriers.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

KIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(KIT))
import launch  # noqa: E402
import provenance  # noqa: E402

FAKE = KIT / "tests" / "fake_claude.py"
AGENT, TASK = "beth-mathematician", "T01_heat_weak_form"


class Scenarios(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.log = Path(self.tmp.name) / "fake.jsonl"

    def tearDown(self):
        self.tmp.cleanup()

    def run_scenario(self, scenario):
        with mock.patch.dict(os.environ, {"AGENTKIT_FAKE_SCENARIO": scenario, "AGENTKIT_FAKE_LOG": str(self.log)}):
            verdict, run_dir = launch.run_agent(AGENT, TASK, runs_root=self.tmp.name, fake=str(FAKE))
        manifest = json.loads((run_dir / "manifest.json").read_text())
        notes = [json.loads(s) for s in self.log.read_text().splitlines()]
        return verdict, run_dir, manifest, notes

    def check_record(self, run_dir, manifest):
        ok, msg, last = provenance.verify(run_dir / "PROVENANCE.md")
        self.assertTrue(ok, msg)
        self.assertEqual(last, manifest["provenance_chain"])
        for f in ("settings.json", "transcript.jsonl", "workspace/task.json"):
            self.assertTrue((run_dir / f).is_file(), f)

    def test_pass(self):
        verdict, run_dir, manifest, notes = self.run_scenario("pass")
        self.assertEqual(verdict, "PASS", manifest["problems"])
        self.check_record(run_dir, manifest)
        self.assertTrue(json.loads((run_dir / "grade.json").read_text())["pass"])
        argv = notes[0]["argv"]
        self.assertEqual(argv[argv.index("--setting-sources") + 1], "")      # no CLAUDE.md, no other settings
        self.assertEqual(argv[argv.index("--model") + 1], "claude-opus-5-5")
        self.assertIn("You are the mathematician", argv[argv.index("--system-prompt") + 1])
        self.assertEqual(manifest["guard_blocks"], 0)
        schema = json.loads(argv[argv.index("--json-schema") + 1])
        self.assertNotIn("$schema", schema)      # found in the first real call: Claude Code rejects it

    def test_wrong_answer_fails(self):
        verdict, run_dir, manifest, _ = self.run_scenario("wrong")
        self.assertEqual(verdict, "FAIL")
        self.check_record(run_dir, manifest)

    def test_blocked_abstains(self):
        verdict, run_dir, manifest, _ = self.run_scenario("blocked")
        self.assertEqual(verdict, "ABSTAIN")
        self.assertFalse((run_dir / "grade.json").exists())
        self.check_record(run_dir, manifest)

    def test_no_answer_is_an_error(self):
        verdict, _, manifest, _ = self.run_scenario("no_answer")
        self.assertEqual(verdict, "ERROR")
        self.assertIn("no structured answer", manifest["problems"])

    def test_forbidden_attempts_are_blocked_and_recorded(self):
        verdict, run_dir, manifest, notes = self.run_scenario("probe")
        self.assertEqual(verdict, "PASS", manifest["problems"])     # blocked attempts are not a failure ...
        self.assertEqual(manifest["guard_blocks"], 4)                # ... but each one is recorded
        self.assertEqual([n["blocked"] for n in notes if n.get("tool")][:4], [True] * 4)
        self.assertFalse((run_dir / "PROVENANCE.md").read_text().count("forged"))
        self.assertIn("guard blocked Bash", (run_dir / "PROVENANCE.md").read_text())
        self.check_record(run_dir, manifest)

    def test_a_failed_barrier_is_caught_by_the_audit(self):
        verdict, _, manifest, _ = self.run_scenario("bypass")
        self.assertEqual(verdict, "ERROR")
        self.assertTrue(any(p.startswith("BARRIER FAILURE: Read") for p in manifest["problems"]))

    def test_command_line_with_relative_paths(self):
        """Found in the first hand run: a relative --fake path broke once the agent ran in its own folder."""
        with tempfile.TemporaryDirectory() as d:
            env = dict(os.environ, AGENTKIT_FAKE_SCENARIO="pass")
            p = subprocess.run([sys.executable, "scripts/agentkit/launch.py", AGENT, "--task", TASK,
                                "--runs-root", os.path.relpath(d, launch.ROOT),
                                "--fake", "scripts/agentkit/tests/fake_claude.py"],
                               cwd=launch.ROOT, capture_output=True, text=True, env=env)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(p.stdout.startswith("PASS"))

    def test_silent_permission_mode_switch_is_an_error(self):
        with mock.patch.dict(os.environ, {"AGENTKIT_FAKE_MODE": "default"}):
            verdict, _, manifest, _ = self.run_scenario("pass")
        self.assertEqual(verdict, "ERROR")
        self.assertTrue(any("permission mode 'default'" in p for p in manifest["problems"]))

    def test_allowed_tools_passed_explicitly(self):
        _, _, _, notes = self.run_scenario("pass")
        argv = notes[0]["argv"]
        self.assertIn("Bash(python scratch/*)", argv[argv.index("--allowedTools") + 1])

    def test_stop_hook_sends_the_agent_back(self):
        verdict, _, manifest, notes = self.run_scenario("forgets")
        stops = [n for n in notes if n.get("stop")]
        self.assertEqual([s["blocked"] for s in stops], [True, False])
        self.assertEqual(verdict, "PASS", manifest["problems"])


class Gates(unittest.TestCase):
    def test_real_call_refused_while_enforcers_unbuilt(self):
        planned = ["N9 (a rule): planned: step 3"]
        with mock.patch.object(launch, "unbuilt_enforcers", return_value=planned), \
                mock.patch.object(subprocess, "run", side_effect=AssertionError("no process may start")):
            with self.assertRaises(launch.LaunchError) as cm:
                launch.run_agent(AGENT, TASK, runs_root=tempfile.mkdtemp())
        self.assertIn("unbuilt", str(cm.exception))

    def test_unbuilt_enforcers_reads_the_status(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "tests").mkdir()
            (Path(d) / "tests" / "enforcers.json").write_text(json.dumps({"enforcers": [
                {"id": "N1", "rule": "a", "status": "built (step 2)"},
                {"id": "N2", "rule": "b", "status": "prompt only"},
                {"id": "N3", "rule": "c", "status": "tools: built; sandbox: planned step 3"}]}))
            self.assertEqual(launch.unbuilt_enforcers(Path(d)), ["N3 (c): tools: built; sandbox: planned step 3"])
            self.assertEqual(launch.unbuilt_enforcers(Path(d) / "missing"), ["tests/enforcers.json is missing"])

    def write_redteam(self, runs, run_id, barrier, verdict="PASS", fake=False):
        d = Path(runs) / f"{run_id}_redteam"
        d.mkdir(parents=True)
        (d / "manifest.json").write_text(json.dumps({"verdict": verdict, "fake": fake, "barrier": barrier}))

    def test_redteam_gap_needs_a_real_pass_on_the_same_barriers(self):
        fp = {"scripts/agentkit/guard.py": "aaa", "(Claude Code) version": "2.1.326"}
        with tempfile.TemporaryDirectory() as runs:
            self.assertIn("no red-team PASS", launch.redteam_gap(runs, fp)[0])
            self.write_redteam(runs, "20261117-000001", fp, verdict="ERROR")
            self.write_redteam(runs, "20261117-000002", fp, fake=True)
            self.assertIn("no red-team PASS", launch.redteam_gap(runs, fp)[0])   # a failed or fake run never counts
            self.write_redteam(runs, "20261117-000003", fp)
            self.assertEqual(launch.redteam_gap(runs, fp), [])
            changed = dict(fp, **{"scripts/agentkit/guard.py": "bbb"})
            gap = launch.redteam_gap(runs, changed)
            self.assertIn("20261117-000003", gap[0])
            self.assertIn("scripts/agentkit/guard.py", gap[0])
            self.assertNotIn("version", gap[0])

    def test_fingerprint_covers_barriers_not_prompt_text(self):
        definition, folder, _ = launch.find_agent(AGENT)
        agent = launch.agentlib.load_agent(definition)
        fp = launch.barrier_fingerprint(folder, agent, "v1")
        for key in ("scripts/agentkit/guard.py", "scripts/agentkit/launch.py", "scripts/agentkit/base_settings.json",
                    f"{launch.agentlib.rel(folder)}/settings.json", f"{launch.agentlib.rel(folder)}/guard.json",
                    f"{launch.agentlib.rel(folder)}/env.allow", "(Claude Code) version"):
            self.assertIn(key, fp)
        self.assertNotIn(launch.agentlib.rel(definition), fp)
        self.assertEqual(fp, launch.barrier_fingerprint(folder, dict(agent, system_prompt="other text"), "v1"))
        self.assertNotEqual(fp, launch.barrier_fingerprint(folder, dict(agent, tools=agent["tools"] + ["WebFetch"]), "v1"))
        self.assertNotEqual(fp, launch.barrier_fingerprint(folder, agent, "v2"))

    def test_real_call_refused_without_a_current_redteam_pass(self):
        """A stand-in not named fake_claude.py is treated as the real thing: the gates apply."""
        with tempfile.TemporaryDirectory() as d:
            stub = Path(d) / "claude_stub.py"
            stub.write_text(FAKE.read_text())
            runs = Path(d) / "runs"
            with mock.patch.object(launch, "unbuilt_enforcers", return_value=[]):
                with self.assertRaises(launch.LaunchError) as cm:
                    launch.run_agent(AGENT, TASK, runs_root=runs, fake=str(stub))
                self.assertIn("red-team", str(cm.exception))
                self.assertFalse(runs.exists())                      # refused before any run folder or call
                definition, folder, _ = launch.find_agent(AGENT)
                fp = launch.barrier_fingerprint(folder, launch.agentlib.load_agent(definition),
                                                "0.0.0 (fake_claude, tests only)")
                self.write_redteam(runs, "20261117-000001", fp)
                verdict, _ = launch.run_agent(AGENT, TASK, runs_root=runs, fake=str(stub))
            self.assertEqual(verdict, "PASS")

    def test_name_must_be_unique(self):
        with tempfile.TemporaryDirectory() as d:
            for wf in ("a", "b"):
                (Path(d) / "workflows" / wf / "agents").mkdir(parents=True)
                (Path(d) / "workflows" / wf / "agents" / "beth-x.md").write_text("---\nname: beth-x\n---\n")
            with self.assertRaises(launch.LaunchError):
                launch.find_agent("beth-x", root=Path(d))

    def test_unknown_task_refused(self):
        with self.assertRaises(launch.LaunchError):
            launch.run_agent(AGENT, "T99_none", runs_root=tempfile.mkdtemp(), fake=str(FAKE))


class Environment(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.rdir = Path(self.tmp.name)
        self.ws = self.rdir / "workspace"
        self.ws.mkdir()
        _, self.folder, _ = launch.find_agent(AGENT)

    def tearDown(self):
        self.tmp.cleanup()

    def test_only_listed_variables_pass(self):
        secrets = {"SSH_AUTH_SOCK": "/tmp/agent.sock", "GIT_ASKPASS": "/x", "ANTHROPIC_API_KEY": "sk-test",
                   "GITHUB_TOKEN": "t", "HOME": "/Users/test"}
        with mock.patch.dict(os.environ, secrets):
            env = launch.build_env(self.ws, self.rdir, {"bash_allow": []}, self.folder)
        for k in ("SSH_AUTH_SOCK", "GIT_ASKPASS", "ANTHROPIC_API_KEY", "GITHUB_TOKEN"):
            self.assertNotIn(k, env)
        self.assertEqual(env["HOME"], "/Users/test")
        self.assertTrue(env["PATH"].startswith(str(launch.PY.parent)))
        self.assertTrue(env["TMPDIR"].startswith(str(self.ws)))

    def test_env_allow_may_not_set_launcher_variables(self):
        f = Path(self.tmp.name) / "beth-y"
        f.mkdir()
        (f / "env.allow").write_text("HOME\nPATH\n")
        with self.assertRaises(launch.LaunchError):
            launch.env_names(f)

    def test_env_allow_is_required(self):
        with self.assertRaises(launch.LaunchError):
            launch.env_names(Path(self.tmp.name))


class Settings(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name) / "beth-x"
        self.folder.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def compose(self, own):
        (self.folder / "settings.json").write_text(json.dumps(own))
        return launch.compose_settings(self.folder, Path(self.tmp.name) / "out.json")

    def test_base_denies_always_present(self):
        s = self.compose({"permissions": {"allow": ["Read(./**)"]}})
        deny = s["permissions"]["deny"]
        self.assertIn("Read(~/.ssh/**)", deny)
        self.assertIn(f"Read(/{launch.ROOT}/docs/planning/**)", deny)
        self.assertIn("Edit(**/PROVENANCE.md)", deny)
        self.assertEqual(s["permissions"]["allow"], ["Read(./**)"])
        self.assertIn("guard.py", s["hooks"]["PreToolUse"][0]["hooks"][0]["command"])
        self.assertIn("stop_check.py", s["hooks"]["Stop"][0]["hooks"][0]["command"])

    def test_agent_cannot_set_hooks_or_modes(self):
        with self.assertRaises(launch.LaunchError):
            self.compose({"permissions": {}, "hooks": {}})
        with self.assertRaises(launch.LaunchError):
            self.compose({"permissions": {"defaultMode": "bypassPermissions"}})

    def test_sandbox_comes_from_the_base_and_is_filled_in(self):
        ws = Path(self.tmp.name) / "run" / "workspace"
        s = self.compose_with_workspace({"permissions": {}}, ws)
        sb = s["sandbox"]
        self.assertTrue(sb["enabled"] and sb["failIfUnavailable"])
        self.assertFalse(sb["allowUnsandboxedCommands"])
        self.assertEqual(sb["network"]["allowedDomains"], [])
        self.assertIn(str(ws), sb["filesystem"]["allowRead"])
        self.assertIn(str(launch.git_root()), sb["filesystem"]["denyRead"])
        self.assertNotRegex(json.dumps(sb), r"\{(ROOT|GITROOT|WORKSPACE)\}|//ROOT/")   # all filled

    def test_agent_cannot_set_its_own_sandbox(self):
        with self.assertRaises(launch.LaunchError):
            self.compose({"permissions": {}, "sandbox": {"enabled": False}})

    def compose_with_workspace(self, own, ws):
        (self.folder / "settings.json").write_text(json.dumps(own))
        return launch.compose_settings(self.folder, Path(self.tmp.name) / "out.json", workspace=ws)

    def test_real_agent_settings_compose(self):
        _, folder, _ = launch.find_agent(AGENT)
        s = launch.compose_settings(folder, Path(self.tmp.name) / "out.json")
        self.assertIn("Bash(python scratch/*)", s["permissions"]["allow"])
        self.assertIn("WebFetch", s["permissions"]["deny"])


if __name__ == "__main__":
    unittest.main()
