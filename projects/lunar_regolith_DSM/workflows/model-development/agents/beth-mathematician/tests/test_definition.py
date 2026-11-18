"""Tests of beth-mathematician's definition (doc 28, step 1).

The definition loads with the shared loader; its frontmatter grants exactly the intended tools; the body has
the house sections; and every "Never" in the body has an enforcer entry in enforcers.json (CLAUDE.md, "Where
agents live": every rule needs an enforcer outside the prompt, or is marked "prompt only").

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import json
import re
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
AGENT_DIR = HERE.parent
DEFINITION = AGENT_DIR.parent / "beth-mathematician.md"
ROOT = AGENT_DIR.parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import agentlib  # noqa: E402

TOOLS = ["Read", "Write", "Edit", "Bash"]
MUST_DISALLOW = {"WebFetch", "WebSearch", "Agent", "NotebookEdit"}


def never_bullets(body):
    """The bullets under '- **Never:**', as plain text without the trailing ':' or ';'."""
    m = re.search(r"^- \*\*Never:\*\*\n((?:  - .*\n(?:    .*\n)*)+)", body, re.M)
    if not m:
        return []
    items = re.split(r"^  - ", m.group(1), flags=re.M)
    return [" ".join(i.split()).rstrip(";.") for i in items if i.strip()]


class Definition(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.agent = agentlib.load_agent(DEFINITION)     # checks name = file, full model ID, prefix = family
        cls.meta, cls.body = agentlib._frontmatter(DEFINITION.read_text(), DEFINITION)

    def test_tools_exactly(self):
        self.assertEqual(self.agent["tools"], TOOLS)

    def test_dangerous_tools_disallowed(self):
        self.assertTrue(MUST_DISALLOW <= set(self.agent["disallowed"]), self.agent["disallowed"])

    def test_worker_settings(self):
        self.assertEqual(self.agent["model"], "claude-opus-5-5")
        self.assertEqual(self.meta.get("permissionMode"), "dontAsk")
        self.assertLessEqual(self.agent["max_turns"], 60)
        for key in ("memory", "effort"):          # a fresh worker stays fresh; workflow agents keep default effort
            self.assertNotIn(key, self.meta)

    def test_house_sections(self):
        for heading in ("## Your place in the line", "## When invoked", "- **Never:**"):
            self.assertIn(heading, self.body)

    def test_name_is_unique_in_project(self):
        hits = [p for p in list(ROOT.glob("workflows/*/agents/*.md")) + list(ROOT.glob(".claude/agents/*.md"))
                if p.stem == "beth-mathematician"]
        self.assertEqual([p.relative_to(ROOT).as_posix() for p in hits],
                         ["workflows/model-development/agents/beth-mathematician.md"])


class EveryNeverHasAnEnforcer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, body = agentlib._frontmatter(DEFINITION.read_text(), DEFINITION)
        cls.nevers = never_bullets(body)
        cls.entries = json.loads((HERE / "enforcers.json").read_text())["enforcers"]

    def test_found_the_never_list(self):
        self.assertGreaterEqual(len(self.nevers), 5)

    def test_each_never_has_one_entry(self):
        for rule in self.nevers:
            with self.subTest(rule=rule):
                hits = [e for e in self.entries if e["rule"] in rule]
                self.assertEqual(len(hits), 1, f"'Never {rule}' needs exactly one enforcers.json entry")

    def test_each_entry_has_a_never(self):
        for e in self.entries:
            with self.subTest(id=e["id"]):
                self.assertTrue(any(e["rule"] in rule for rule in self.nevers), f"{e['id']} matches no 'Never'")

    def test_entries_are_complete(self):
        layers = {"environment", "harness", "repository", "prompt"}
        for e in self.entries:
            with self.subTest(id=e["id"]):
                self.assertTrue(e["enforcer"].strip())
                self.assertTrue(set(e["layer"]) <= layers and e["layer"])
                if e["layer"] == ["prompt"]:
                    self.assertIn("prompt only", e["enforcer"])   # says so, and names its backstop
                    self.assertIn("backstop", e["enforcer"])

    def test_built_tool_enforcers_hold(self):
        """N3 and N4 are enforced by the frontmatter today; check that they really are there."""
        agent = agentlib.load_agent(DEFINITION)
        self.assertIn("WebFetch", agent["disallowed"])
        self.assertIn("WebSearch", agent["disallowed"])
        self.assertIn("Agent", agent["disallowed"])


if __name__ == "__main__":
    unittest.main()
