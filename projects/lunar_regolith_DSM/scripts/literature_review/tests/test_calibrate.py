"""Tests for calibrate.py: the bounds, the split, the pass rules and the ledger.

    venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import calibrate as cal  # noqa: E402


def plant(i, kind="dropped_condition", source="script", copy_kind="pdf"):
    return {"id": f"P{i:03d}", "kind": kind, "source": source, "copy_kind": copy_kind, "paper": "x"}


class Bounds(unittest.TestCase):
    def test_values_quoted_in_the_design(self):
        self.assertAlmostEqual(cal.cp_lower(58, 60), 0.899, places=3)   # WORKFLOW.md §5
        self.assertAlmostEqual(cal.cp_lower(19, 19), 0.854, places=3)
        self.assertAlmostEqual(cal.cp_lower(10, 10), 0.741, places=3)
        self.assertAlmostEqual(cal.cp_upper(0, 60), 0.049, places=3)    # U3b, k = 0
        self.assertAlmostEqual(cal.cp_upper(1, 60), 0.077, places=3)
        self.assertAlmostEqual(cal.cp_upper(2, 60), 0.101, places=3)

    def test_eighteen_of_eighteen_cannot_pass(self):
        self.assertLess(cal.cp_lower(18, 18), cal.BOUND)


class Split(unittest.TestCase):
    def test_user_plants_always_held_out_and_split_is_repeatable(self):
        ps = [plant(i, source="user") for i in range(5)] + [plant(i) for i in range(5, 60)]
        h1, t1 = cal.split(ps, 30, seed=7)
        h2, t2 = cal.split(ps, 30, seed=7)
        self.assertEqual((h1, t1), (h2, t2))
        self.assertTrue(all(f"P{i:03d}" in h1 for i in range(5)))
        self.assertEqual(len(h1), 30)
        self.assertEqual(set(h1) & set(t1), set())

    def test_split_spreads_kinds(self):
        ps = [plant(i, kind="a") for i in range(20)] + [plant(i, kind="b") for i in range(20, 40)]
        h, _ = cal.split(ps, 10, seed=1)
        kinds = [p["kind"] for p in ps if p["id"] in h]
        self.assertEqual(kinds.count("a"), 5)

    def test_plants_from_one_statement_stay_together(self):
        ps = [dict(plant(i, kind=k), from_statement=i // 2, correct_page=i // 2)
              for i, k in zip(range(40), ["changed_digit", "wrong_unit"] * 20)]
        h, t = cal.split(ps, 10, seed=3)
        side = {p["id"]: p["id"] in h for p in ps}
        for s in range(20):
            self.assertEqual(side[f"P{2 * s:03d}"], side[f"P{2 * s + 1:03d}"])

    def test_plants_on_a_user_plant_page_are_held_out(self):
        ps = [dict(plant(0, source="user"), planted_page=7)]
        ps += [dict(plant(i), from_statement=i, correct_page=i) for i in range(1, 30)]
        h, t = cal.split(ps, 5, seed=3)
        self.assertIn("P007", h)
        self.assertNotIn("P007", t)


class Score(unittest.TestCase):
    def setUp(self):
        self.ps = [plant(i) for i in range(20)] + [plant(20, kind="wrong_page"),
                                                   plant(21, kind="equation_sign_or_exponent")]
        self.held = [p["id"] for p in self.ps]
        self.controls = [{"control_id": f"C{i}", "flagged": i == 0} for i in range(10)]

    def results(self, missed=()):
        r = [{"plant_id": p["id"], "caught_by": None if p["id"] in missed else
              ("E3" if p["kind"] == "wrong_page" else "E4")} for p in self.ps]
        return r + self.controls

    def test_all_caught_passes(self):
        rep = cal.score(self.ps, self.held, self.results())
        self.assertTrue(rep["passed"], rep["gates"])

    def test_a_missed_sign_change_fails(self):
        rep = cal.score(self.ps, self.held, self.results(missed={"P021"}))
        self.assertFalse(rep["passed"])

    def test_a_missed_rule_kind_fails(self):
        rep = cal.score(self.ps, self.held, self.results(missed={"P020"}))
        self.assertFalse(rep["passed"])

    def test_one_judgment_miss_in_21_fails_the_bound(self):
        rep = cal.score(self.ps, self.held, self.results(missed={"P003"}))   # 20 of 21: 80.6 %
        self.assertFalse(rep["passed"])

    def test_too_many_flagged_controls_fail(self):
        self.controls = [{"control_id": f"C{i}", "flagged": i < 2} for i in range(10)]
        rep = cal.score(self.ps, self.held, self.results())
        self.assertFalse(rep["passed"])

    def test_abstention_cannot_raise_the_score(self):
        """A plant or control without a verdict counts against the checker (a miss, a flag)."""
        r = self.results()
        r = [dict(x, caught_by=None, no_verdict=True) if x.get("plant_id") == "P003" else x for x in r]
        self.controls = [{"control_id": f"C{i}", "flagged": i == 0 or i == 1, **({"no_verdict": True} if i == 1
                                                                                   else {})} for i in range(10)]
        r = [x for x in r if "control_id" not in x] + self.controls
        rep = cal.score(self.ps, self.held, r)
        self.assertFalse(rep["passed"])                    # one miss in 21 fails the bound; 2 of 10 controls fail
        self.assertEqual(rep["no_verdict"], {"plants": 1, "controls": 1})

    def test_misses_reported_only_as_counts(self):
        rep = cal.score(self.ps, self.held, self.results(missed={"P003"}))
        self.assertEqual(rep["misses_per_kind"], {"dropped_condition": 1})
        self.assertNotIn("P003", json.dumps(rep))


class Ledger(unittest.TestCase):
    def test_a_failed_set_is_never_scored_again(self):
        d = tempfile.mkdtemp()
        split = os.path.join(d, "split.yaml")
        ledger = os.path.join(d, "ledger.json")
        with open(split, "w") as f:
            f.write("heldout_set: abc\nheldout: []\ntuning: []\n")
        with open(ledger, "w") as f:
            json.dump([{"heldout_set": "abc", "prompt_version": "v1", "passed": False, "date": "x"}], f)
        res = os.path.join(d, "r.jsonl")
        open(res, "w").close()
        for v in ("v1", "v2"):
            with self.assertRaises(SystemExit):
                cal.main(["score", "--results", res, "--prompt-version", v, "--out", d,
                          "--split", split, "--ledger", ledger])


if __name__ == "__main__":
    unittest.main()
