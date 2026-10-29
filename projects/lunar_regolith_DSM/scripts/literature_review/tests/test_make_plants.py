"""Tests for make_plants.py: each plant changes exactly one thing, and identifiers are kept.

    venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests
"""
import os
import random
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import make_plants as mp  # noqa: E402


class MakePlants(unittest.TestCase):
    def test_identifiers_and_references_are_not_values(self):
        for t in ["in the SA81 model", "Fig. 1 shows", "equation (31)", "Table 2 lists",
                  "on the order of $10^{90}$ years"]:
            self.assertEqual(mp.value_numbers(t), [], t)

    def test_values_are_found(self):
        got = [n[2] + n[3] for n in mp.value_numbers("50 microns at $-55^\\circ$C after ~1,000 years")]
        self.assertEqual(got, ["50", "-55", "1,000"])

    def test_changed_digit_changes_one_digit(self):
        rng = random.Random(1)
        for _ in range(50):
            new, _ = mp.changed_digit("a grain size of 50 microns", rng)
            a, b = "a grain size of 50 microns", new
            self.assertEqual(len(a), len(b))
            self.assertEqual(sum(x != y for x, y in zip(a, b)), 1)

    def test_value_digit_keeps_the_exponent(self):
        rng = random.Random(2)
        for _ in range(50):
            new, _ = mp.changed_value_digit("1.8e6", rng)
            self.assertTrue(new.endswith("e6") and new != "1.8e6" and not new.startswith("0"))

    def test_typo_is_carried_to_every_plant(self):
        e = {"paper": "li_2009", "page": 85, "kind": "equation", "statement": "a term $-2 u$ with $10^{3}$",
             "typo": "the sign should be +"}
        plants = mp.plants_for(1, e, {}, random.Random(4))
        self.assertTrue(plants)
        self.assertTrue(all(p["source_typo"] == "the sign should be +" for p in plants))

    def test_exponent_and_sign_flips(self):
        rng = random.Random(3)
        t, v, note = mp.sign_or_exponent("order of $10^{90}$ years", None, rng)
        self.assertIn("10^{-90}", t)
        t, v, note = mp.sign_or_exponent("text", "1.8e6", rng)
        self.assertEqual(v, "1.8e-6")
        t, v, note = mp.sign_or_exponent("around $-55^\\circ$C", None, rng)
        self.assertIn("$55^\\circ$C", t)

    def test_unit_swap_keeps_latex_balanced(self):
        t, _, _ = mp.wrong_unit("at $T=-5^\\circ$C and", None, random.Random(4))
        self.assertEqual(t.count("$") % 2, 0)
        self.assertIn("$T=-5$ K", t)

    def test_unit_field_swap(self):
        _, u, _ = mp.wrong_unit("Table I gives it.", "J m^-3 K^-1", random.Random(5))
        self.assertNotEqual(u, "J m^-3 K^-1")

    def test_wrong_page_moves_the_page(self):
        p, _ = mp.wrong_page(250, None, random.Random(6))
        self.assertIn(abs(p - 250), (1, 2))


if __name__ == "__main__":
    unittest.main()
