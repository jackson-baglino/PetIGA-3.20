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

    def test_number_found_at_its_place_not_inside_a_citation_id(self):
        # "2" also occurs in "2009"; the tag must be looked for after this "2", not after the one in the id.
        claims = [{"id": "li_2009-c003", "text": "The width was 16.74 grid spacings and the deviation less than 2 deg."}]
        md = ("## 2. Methods\n\nWith a width of 16.74 [li_2009-c003 width] spacings the angle deviated "
              "by less than 2 [li_2009-c003 deviation] degrees.\n")
        r = cx.check(md, claims)
        self.assertEqual((r["untagged_numbers"], r["tag_mismatch"]), ([], []))

    def test_bold_labels_and_intervals_are_not_checked(self):
        md = ("## 2. Methods\n\n**Folch and Plapp: what it does.**\n\n"
              "- The bounds [0,1] are not imposed [li_2009-c002].\n")
        r = cx.check(md, CLAIMS)
        self.assertTrue(r["pass"], r)
        self.assertEqual(r["sentences"], 1)

    def test_bold_text_inside_a_sentence_is_still_checked(self):
        md = "## 2. Methods\n\n**Note.** The method is exact in every case.\n"
        self.assertEqual(cx.check(md, CLAIMS)["uncited"], [1])


if __name__ == "__main__":
    unittest.main()
