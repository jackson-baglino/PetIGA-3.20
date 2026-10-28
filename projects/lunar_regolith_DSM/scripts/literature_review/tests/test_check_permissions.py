"""Tests for check_permissions.py: the real files pass, and planted faults fail.

    venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import check_permissions as cp  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
PERM = os.path.join(".claude", "permissions", "literature-review")


class CheckPermissions(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        shutil.copytree(os.path.join(ROOT, PERM), os.path.join(self.tmp, PERM))
        for f in ("settings.json", "settings.local.json"):
            src = os.path.join(ROOT, ".claude", f)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(self.tmp, ".claude", f))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def edit(self, rel, fn):
        p = os.path.join(self.tmp, rel)
        with open(p) as f:
            d = json.load(f)
        fn(d)
        with open(p, "w") as f:
            json.dump(d, f)

    def faults(self):
        return cp.collect(self.tmp)[1]

    def test_real_files_pass(self):
        self.assertEqual(self.faults(), [])

    def test_reader_with_web_and_private_inputs_fails(self):
        self.edit(os.path.join(PERM, "agents", "morty-reader.settings.json"),
                  lambda d: d["permissions"]["allow"].append("WebSearch"))
        self.edit(os.path.join(PERM, "inputs.json"),
                  lambda d: d["agents"]["morty-reader"].append("user_private"))
        self.assertTrue(any("morty-reader: holds all three" in f for f in self.faults()))

    def test_scout_given_the_brief_fails(self):
        self.edit(os.path.join(PERM, "inputs.json"),
                  lambda d: d["agents"]["morty-scout"].append("user_private"))
        self.assertTrue(any("morty-scout: holds all three" in f for f in self.faults()))

    def test_sandbox_that_fails_open_fails(self):
        self.edit(os.path.join(PERM, "agents", "beth-checker.settings.json"),
                  lambda d: d["sandbox"].pop("failIfUnavailable"))
        self.assertTrue(any("fail open" in f for f in self.faults()))

    def test_agent_allowed_to_write_fails(self):
        self.edit(os.path.join(PERM, "agents", "jerry-screener.settings.json"),
                  lambda d: d["permissions"]["deny"].remove("Write"))
        self.assertTrue(any("jerry-screener: Write must be denied" in f for f in self.faults()))

    def test_review_mode_without_hook_fails(self):
        self.edit(os.path.join(PERM, "sessions", "review-mode.json"), lambda d: d.pop("hooks"))
        self.assertTrue(any("allow-list hook" in f for f in self.faults()))

    def test_review_mode_missing_store_deny_fails(self):
        self.edit(os.path.join(PERM, "sessions", "review-mode.json"),
                  lambda d: d["permissions"]["deny"].remove("Read(/library/papers/**)"))
        self.assertTrue(any("review-mode: missing deny Read(/library/papers/**)" in f
                            for f in self.faults()))

    def test_normal_profile_without_acceptance_fails(self):
        self.edit(os.path.join(PERM, "accepted_risks.json"), lambda d: d.pop("normal"))
        self.assertTrue(any("normal: holds all three" in f for f in self.faults()))

    def test_normal_profile_missing_store_deny_fails(self):
        self.edit(os.path.join(".claude", "settings.json"),
                  lambda d: d["permissions"]["deny"].remove("Read(/index/**)"))
        self.assertTrue(any("normal: missing deny Read(/index/**)" in f for f in self.faults()))

    def test_agent_without_input_classes_fails(self):
        shutil.copy(os.path.join(self.tmp, PERM, "agents", "jerry-screener.settings.json"),
                    os.path.join(self.tmp, PERM, "agents", "jerry-extra.settings.json"))
        self.assertTrue(any("jerry-extra: no input classes" in f for f in self.faults()))

    def test_launcher_profile_must_start_closed(self):
        p = os.path.join(self.tmp, PERM, "profiles", "offline.sb")
        with open(p, "w") as f:
            f.write("(version 1)\n(allow default)\n")
        self.assertTrue(any("profiles/offline.sb" in f for f in self.faults()))


if __name__ == "__main__":
    unittest.main()
