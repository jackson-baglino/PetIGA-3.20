"""check_plants.py: a gold statement listing several pages is checked on all of them."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import check_plants as cp  # noqa: E402


class Pages(unittest.TestCase):
    def test_every_listed_page_is_shown_and_named(self):
        ctl = {"id": "G004", "text": "t", "page_pdf": 5, "pages_pdf": [5, 3, 6]}
        plant = {"id": "P1", "text": "p", "page_pdf": 8}
        self.assertEqual(cp.pages_shown([ctl, plant]), [3, 5, 6, 8])
        c = cp.shown_claim(ctl)
        self.assertEqual((c["page_pdf"], c["location"]), (5, "the statement also draws on PDF page(s) 3, 6"))
        self.assertEqual(cp.shown_claim(plant)["location"], "")

    def test_a_single_page_gold_statement(self):
        self.assertEqual(cp.pages_shown([{"id": "G1", "page_pdf": 4, "pages_pdf": [4]}]), [4])


if __name__ == "__main__":
    unittest.main()
