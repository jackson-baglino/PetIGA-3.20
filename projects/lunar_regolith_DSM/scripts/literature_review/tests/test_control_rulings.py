"""control_rulings.py: ruling categories and note detection, without reading any text out."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import control_rulings as cr  # noqa: E402

MD = """# Controls the checker flagged (tuning set, calib-dev6)

## G001 (paper_a, printed p. 4)

**Reason:** quoted model text that mentions Ruling: and Notes: inline.

Ruling: gold wrong. Needs the supplement.

Notes: the fact is only in the supplementary material.

## G004 (paper_b, printed p. 7)

Ruling: checker wrong

Notes:
first line of a note on the next line

## G007 (paper_c, printed p. 2)

Ruling: **Gold** wrong

Notes:

## G014 (paper_d, printed p. 9)

Ruling:
"""


class TestControlRulings(unittest.TestCase):
    def test_categories(self):
        self.assertEqual(cr.rulings(MD), {"G001": "gold wrong", "G004": "checker wrong", "G007": "gold wrong",
                                          "G014": "blank"})

    def test_notes_same_line_or_below(self):
        self.assertEqual(cr.noted(MD), {"G001", "G004"})

    def test_header_variants(self):
        for head in ("**Notes:** x", "**Notes**: x", "Note: x", "### Notes\nx", "  notes: x", "NOTES: x"):
            md = "## G001 (p, printed p. 1)\n\nRuling: gold wrong\n\n" + head + "\n"
            self.assertEqual(cr.noted(md), {"G001"}, head)

    def test_note_that_in_reason_is_not_a_header(self):
        md = "## G001 (p, printed p. 1)\n\nNote that the page omits it.\n\nRuling: gold wrong\n"
        self.assertEqual(cr.noted(md), set())

    def test_empty_notes_and_old_files(self):
        old = MD.replace("Notes:", "")
        self.assertEqual(cr.noted(old), set())


if __name__ == "__main__":
    unittest.main()
