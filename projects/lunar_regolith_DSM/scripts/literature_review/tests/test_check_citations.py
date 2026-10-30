"""Tests for check_citations.py (Y2)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import check_citations as cx  # noqa: E402

CLAIMS = [{"id": "li_2009-c001", "text": "The interface width is 0.5 in the tests."},
          {"id": "li_2009-c002", "text": "The method converges at second order."}]


class Citations(unittest.TestCase):
    def test_good_memo_passes(self):
        md = ("## 1. Summary\n\nThe diffuse domain method converges at second order [li_2009-c002]. "
              "Its tests used a width of 0.5 [li_2009-c001 interface width].\n\n"
              "## 3. Model options\n\nWe could try the method with a width of 3.\n")
        self.assertTrue(cx.check(md, CLAIMS)["pass"])

    def test_failures_are_found(self):
        md = ("## 1. Summary\n\nThe method is exact in every case.\n\n"
              "The width was 0.7 [li_2009-c001 interface width] in tests [li_2009-c009].\n")
        r = cx.check(md, CLAIMS)
        self.assertEqual(r["uncited"], [1])
        self.assertEqual(r["unknown_ids"], [2])
        self.assertEqual(r["tag_mismatch"], [2])
        self.assertFalse(r["pass"])

    def test_equation_numbers_are_not_values(self):
        md = "## 2. Methods\n\nThe constraint in Eq. (2.21) sums the fields [li_2009-c002].\n"
        self.assertTrue(cx.check(md, CLAIMS)["pass"])


if __name__ == "__main__":
    unittest.main()
