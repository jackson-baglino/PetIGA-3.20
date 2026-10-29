"""Tests for prepare.py (P1): tiles on synthetic page images, and Marker output parsing.

    venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests
"""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import prepare as pr  # noqa: E402

W, H = 2550, 3300                   # a letter page at 300 dpi


def page():
    return np.full((H, W), 255, dtype=np.uint8)


def lines(img, x0, x1, y0, y1, pitch=50, height=30):
    """Text-like lines: dark bars of the given height every pitch rows."""
    for y in range(y0, y1, pitch):
        img[y:y + height, x0:x1] = 0
    return img


def boxes_ok(test, img, boxes):
    ink = img < pr.INK
    test.assertEqual(pr.coverage(img, [{"box_px": b} for b in boxes]), 1.0)
    for x0, y0, x1, y1 in boxes:
        test.assertGreater(y1, y0)
        test.assertLessEqual(y1 - y0, pr.MAX_EDGE)
        # A cut never falls inside a line: the row just inside each edge has no ink that
        # continues across the edge, unless the tile was cut with an overlap.
        if y0 > 0 and ink[y0, x0:x1].any() and ink[y0 - 1, x0:x1].any():
            test.assertTrue(any(b[3] - pr.OVERLAP <= y0 < b[3] for b in boxes if b[0] == x0),
                            f"cut through a line at row {y0}")


class Tiles(unittest.TestCase):
    def test_single_column_is_cut_in_halves_at_gaps(self):
        img = lines(page(), 300, 2250, 300, 3000)
        boxes = pr.plan_tiles(img)
        self.assertEqual(len(boxes), 2)
        boxes_ok(self, img, boxes)

    def test_two_columns_are_tiled_per_column(self):
        img = lines(page(), 300, 1200, 300, 3000)
        lines(img, 1350, 2250, 325, 3000)
        boxes = pr.plan_tiles(img)
        self.assertEqual(len(boxes), 4)
        self.assertTrue(all(b[2] <= 1350 or b[0] >= 1200 for b in boxes))
        boxes_ok(self, img, boxes)

    def test_title_across_the_gutter_is_a_wide_band(self):
        img = lines(page(), 300, 2250, 200, 300)
        lines(img, 300, 1200, 400, 3000)
        lines(img, 1350, 2250, 425, 3000)
        boxes = pr.plan_tiles(img)
        wide = [b for b in boxes if b[0] < 1275 < b[2]]
        self.assertEqual(len(wide), 1)
        self.assertLess(wide[0][3], 400)
        boxes_ok(self, img, boxes)

    def test_centered_equations_do_not_make_a_gutter(self):
        img = lines(page(), 300, 2250, 300, 1500)
        lines(img, 900, 1650, 1600, 1800, pitch=100, height=60)       # displayed equations
        lines(img, 300, 2250, 1900, 3000)
        ink = img < pr.INK
        self.assertIsNone(pr.find_gutter(ink, 296, 2254))

    def test_a_block_without_gaps_is_cut_with_overlap(self):
        img = page()
        img[100:3200, 300:2250:7] = 0                                 # a figure: ink on every row
        boxes = pr.plan_tiles(img)
        self.assertGreater(len(boxes), 1)
        for a, b in zip(boxes, boxes[1:]):
            self.assertEqual(a[3] - b[1], pr.OVERLAP)
        boxes_ok(self, img, boxes)

    def test_blank_page_has_no_tiles(self):
        self.assertEqual(pr.plan_tiles(page()), [])


class MarkerOutput(unittest.TestCase):
    def test_pages_equations_and_headings(self):
        md = ("\n\n{0}" + "-" * 48 + "\n\n# Title\n\ntext\n\n"
              "{1}" + "-" * 48 + "\n\n## 2. Model\n\n$$a = b \\tag{3}$$\n\n$$c = d$$ (4)\n\n$$e$$\n\n$$f$$ (2.14)\n\n$$g$$ (A.3)\n\n$$h$$\n in $\\Omega$ , (2.21)\n\n$$k$$\n\nsee the results in (2.5) and the text that follows it\n")
        pages = pr.split_pages(md)
        self.assertEqual(sorted(pages), [1, 2])
        eqs = pr.equations(pages)
        self.assertEqual([(e["page_pdf"], e["number"]) for e in eqs], [(2, "3"), (2, "4"), (2, None), (2, "2.14"), (2, "A.3"), (2, "2.21"), (2, None)])
        self.assertEqual([(s["page_pdf"], s["level"]) for s in pr.sections(pages)], [(1, 1), (2, 2)])


if __name__ == "__main__":
    unittest.main()
