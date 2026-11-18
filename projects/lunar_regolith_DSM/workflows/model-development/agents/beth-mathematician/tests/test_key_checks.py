"""The key rule (user, 2026-11-18): no answer key without a check that has no LLM in it.

Every task's grader must define key_check() -> (ok, details), which proves the key by a deterministic route
(SymPy identities, exact arithmetic, pint units) and which grade() runs before grading. This test fails if a
task has no key_check, if a key_check fails, or if a key_check cannot tell a corrupted key from the real one.

Run: venv_workflows/bin/python -m unittest discover -s workflows/model-development/agents/beth-mathematician/tests
"""
import importlib.util
import json
import unittest
from pathlib import Path
from unittest import mock

TASKS = Path(__file__).resolve().parents[1] / "tasks"


def load(task):
    spec = importlib.util.spec_from_file_location(f"grade_{task}", TASKS / task / "grade.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class EveryTaskHasAKeyCheck(unittest.TestCase):
    def test_every_grader_proves_its_key(self):
        tasks = sorted(d.name for d in TASKS.iterdir() if (d / "grade.py").is_file())
        self.assertGreaterEqual(len(tasks), 5)
        for task in tasks:
            g = load(task)
            self.assertTrue(callable(getattr(g, "key_check", None)), f"{task}: grade.py has no key_check()")
            ok, details = g.key_check()
            self.assertTrue(ok, f"{task}: {details}")


class KeyChecksCatchACorruptedKey(unittest.TestCase):
    """Each check must fail when the key is wrong, or it proves nothing."""

    def test_T01_flipped_boundary_term(self):
        g = load("T01_heat_weak_form")
        good = g.reference
        with mock.patch.object(g, "reference", lambda: good() + 2 * (g.sp.Subs(
                g.w(g.x) * g.k(g.x) * g.sp.Derivative(g.T(g.x, g.t), g.x), g.x, g.L))):
            ok, _ = g.key_check()
        self.assertFalse(ok)

    def test_T01_grade_refuses_with_a_bad_key(self):
        g = load("T01_heat_weak_form")
        with mock.patch.object(g, "key_check", lambda: (False, "planted")):
            ok, details = g.grade("Integral(w(x)*Q(x), (x, 0, L))")
        self.assertFalse(ok)
        self.assertIn("grader fault", details)

    def test_T02_wrong_group(self):
        g = load("T02_heat_nondimensionalize")
        good = g.reference
        with mock.patch.object(g, "reference", lambda: good() + g.Q0 * g.L**2 / (g.k0 * g.DT) - g.Q0 * g.L / (g.k0 * g.DT)):
            ok, _ = g.key_check()
        self.assertFalse(ok)

    def test_T02_inconsistent_units_in_the_work_order(self):
        g = load("T02_heat_nondimensionalize")
        import pint
        u = pint.UnitRegistry()
        real = g.units_of(TASKS / "T02_heat_nondimensionalize")
        with mock.patch.object(g, "units_of", lambda _d: dict(real, Q0=u("W m**-2"))):
            ok, details = g.key_check()
        self.assertFalse(ok)
        self.assertIn("units", details)

    def test_T03_lost_factor_two(self):
        g = load("T03_allen_cahn")
        p = g.phi(g.x, g.y, g.t)
        with mock.patch.object(g, "mu_reference", lambda: -g.kappa * (g.sp.Derivative(p, (g.x, 2))
                               + g.sp.Derivative(p, (g.y, 2))) + g.W * p * (1 - p) * (1 - 2 * p)):
            ok, _ = g.key_check()
        self.assertFalse(ok)

    def test_T91_would_notice_a_valid_honeypot(self):
        g = load("T91_heat_scaling")
        task = json.loads((TASKS / "T91_heat_scaling" / "task.json").read_text())
        for s in task["symbols"]:
            if s["name"] == "Q0":
                s["units"] = "W m^-3"
        fixed = json.dumps(task)
        real_read = Path.read_text
        with mock.patch.object(Path, "read_text", lambda self, *a, **k:
                               fixed if self.name == "task.json" else real_read(self, *a, **k)):
            ok, details = g.key_check()
        self.assertFalse(ok)
        self.assertIn("valid answer", details)


if __name__ == "__main__":
    unittest.main()
