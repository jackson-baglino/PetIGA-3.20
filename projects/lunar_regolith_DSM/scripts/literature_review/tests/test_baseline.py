"""baseline.py: the script-only C0 baseline (WORKFLOW.md §5)."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import baseline as b  # noqa: E402
import calibrate as cal  # noqa: E402

PAGE = "Surface diffusion dominates above −6 ◦C and for necks smaller than 2 µm; d0 = 1.0152 × 10−9 m (Eq. 2.23)."


class Numbers(unittest.TestCase):
    def test_notations_normalize_to_one_value(self):
        for text in ("1.0152e-9 m", "1.0152 × 10−9 m", "1.0152\\times10^{-9} m", "1.0152 x 10^-9 m"):
            v, u = b.numbers(text)[0]
            self.assertAlmostEqual(v, 1.0152e-9, delta=1e-20, msg=text)
            self.assertEqual(u, "m", text)

    def test_units_and_unicode(self):
        self.assertEqual(b.numbers("−6 ◦C and 2 microns"), [(-6.0, "degc"), (2.0, "um")])
        self.assertEqual(b.numbers(b.plain("$-6^\\circ$C and $2\\,\\mu$m")), [(-6.0, "degc"), (2.0, "um")])


class Judge(unittest.TestCase):
    def j(self, text, equations=None, parameter=None):
        return b.judge({"text": text, "parameter": parameter}, PAGE, equations)

    def test_correct_claim_passes(self):
        r = self.j("Surface diffusion dominates above $-6^\\circ$C for necks under 2 microns.")
        self.assertFalse(r["flag"], r)

    def test_changed_digit_and_unit_are_caught(self):
        self.assertTrue(self.j("Surface diffusion dominates above $-7^\\circ$C.")["flag"])
        r = self.j("necks smaller than 2 mm")
        self.assertTrue(r["flag"])
        self.assertEqual(r["checks"], {"numbers": True, "units": False})

    def test_parameter_field_is_checked(self):
        r = self.j("The capillary length is given.", parameter={"quantity": "d_0", "value": "1.0152e-8",
                                                               "unit_printed": "m"})
        self.assertTrue(r["flag"])

    def test_prose_without_numbers_is_not_judged(self):
        r = self.j("The method holds only in one dimension.")
        self.assertEqual(r, {"flag": False, "checks": {}})

    def test_equations_need_p1_content(self):
        claim = "The equation is $u_t = \\nabla^2 u$."
        self.assertEqual(self.j(claim)["checks"]["equations"], "not checkable")
        self.assertFalse(self.j(claim)["flag"])
        from tres import canon
        self.assertFalse(self.j(claim, equations=[canon("u_t = \\nabla^2 u")])["flag"])
        self.assertTrue(self.j(claim, equations=[canon("u_t = -\\nabla^2 u")])["flag"])


class ScoreBaseline(unittest.TestCase):
    def test_baseline_score_never_touches_the_ledger(self):
        d = tempfile.mkdtemp()
        split, ledger, res = (os.path.join(d, n) for n in ("split.yaml", "ledger.json", "r.jsonl"))
        with open(split, "w") as f:
            f.write("heldout_set: abc\nheldout: []\ntuning: []\n")
        with open(res, "w") as f:
            f.write(json.dumps({"control_id": "G001", "flagged": False}) + "\n")
        rc = cal.main(["score", "--baseline", "--results", res, "--prompt-version", "baseline", "--out", d,
                       "--split", split, "--ledger", ledger])
        self.assertEqual(rc, 0)
        self.assertFalse(os.path.exists(ledger))
        self.assertTrue(os.path.exists(os.path.join(d, "baseline.json")))


if __name__ == "__main__":
    unittest.main()
