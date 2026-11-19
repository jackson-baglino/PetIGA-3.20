"""Tests of render_derivation.py: the layout is fixed, and the lint finds each kind of format problem.

Run: venv_workflows/bin/python -m unittest discover -s scripts/agentkit/tests
"""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import render_derivation as rd  # noqa: E402

TASK = {"task_id": "T99_demo", "statement": "Derive the demo.",
        "symbols": [{"name": "phi(x, t)", "latex": "\\phi", "meaning": "phase field", "units": "1"},
                    {"name": "x", "latex": "x", "meaning": "position", "units": "m"},
                    {"name": "t", "latex": "t", "meaning": "time", "units": "s"}],
        "given": ["Derivative(phi(x, t), t)"]}


def step(**kw):
    s = {"id": "S1", "title": "Relax", "why": "The phase field $\\phi$ relaxes at rate $M$.", "from": ["given"],
         "operation": "substitution", "latex": "\\partial_t\\phi = -M\\phi", "check": {"function": "check_s1",
                                                                                     "printed": "0"}}
    s.update(kw)
    return s


def answer(steps=None, **kw):
    a = {"status": "done", "result_expr": "x", "result_latex": "R = \\partial_t\\phi + M\\phi",
         "symbols": [{"name": "M", "latex": "M", "meaning": "mobility", "units": "s^-1", "typical": 1.0}],
         "steps": steps if steps is not None else [step()], "assumptions_used": ["$M$ is constant"],
         "neglected": [{"effect": "noise", "why": "small"}], "questions": []}
    a.update(kw)
    return a


class Lint(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        (Path(self.tmp.name) / "scratch").mkdir()
        (Path(self.tmp.name) / "scratch" / "check.py").write_text("def check_s1():\n    print(0)\n")

    def tearDown(self):
        self.tmp.cleanup()

    def problems(self, a):
        return rd.lint(a, TASK, self.tmp.name)

    def test_clean(self):
        self.assertEqual(self.problems(answer()), [])

    def test_undeclared_symbol(self):
        self.assertTrue(any("undefined" in p and "K" in p
                            for p in self.problems(answer([step(latex="\\partial_t\\phi = -M K")]))))

    def test_undeclared_greek(self):
        self.assertTrue(any("\\gamma" in p for p in self.problems(answer([step(latex="\\gamma\\phi")]))))

    def test_unicode_in_prose(self):
        self.assertTrue(any("Unicode" in p for p in self.problems(answer([step(why="The field φ relaxes here.")]))))

    def test_ascii_name_in_prose(self):
        self.assertTrue(any("outside $...$" in p and "phi" in p
                            for p in self.problems(answer([step(why="The field phi relaxes at rate $M$.")]))))

    def test_missing_check_function(self):
        bad = step(check={"function": "check_s9", "printed": "0"})
        self.assertTrue(any("check_s9" in p for p in self.problems(answer([bad]))))

    def test_no_check_needs_a_reason(self):
        self.assertTrue(any("no check" in p for p in self.problems(answer([step(check={"function": "",
                                                                                     "printed": ""})]))))

    def test_dollar_in_latex(self):
        self.assertTrue(any("$ delimiters" in p for p in self.problems(answer([step(latex="$\\phi$")]))))


class Layout(unittest.TestCase):
    def test_sections_in_order(self):
        text = rd.render(answer(), TASK, "run1")
        order = ["# T99_demo", "## Work order", "## Symbols", "## Given", "## Assumptions", "## Derivation",
                 "### S1. Relax", "## Result", "## Neglected effects"]
        positions = [text.index(h) for h in order]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("| $M$ | mobility | $\\mathrm{s^{-1}}$ | 1 | introduced |", text)
        self.assertIn("- **A1.** $M$ is constant", text)


if __name__ == "__main__":
    unittest.main()
