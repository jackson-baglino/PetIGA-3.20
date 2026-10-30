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
import docparse_run as dp  # noqa: E402

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


def block(kind, x0, y0, x1, y1):
    """A Marker block from a box in pixels (the layout is in points, 72 per inch)."""
    k = 72 / pr.DPI
    return {"type": kind, "bbox": [x0 * k, y0 * k, x1 * k, y1 * k]}


def layout(*blocks):
    return {"bbox": [0, 0, W * 72 / pr.DPI, H * 72 / pr.DPI], "blocks": list(blocks)}


class LayoutTiles(unittest.TestCase):
    def test_figure_keeps_its_axis_labels_and_caption(self):
        img = page()
        img[300:1500, 300:2250:9] = 0                                 # the plot
        img[1530:1560, 600:1900] = 0                                  # axis labels below its box
        lines(img, 300, 2250, 1650, 1800)                             # caption, three lines
        lines(img, 300, 2250, 2000, 3000)                             # text after it
        lay = layout(block("Figure", 300, 300, 2250, 1500), block("Caption", 300, 1650, 2250, 1780),
                     block("Text", 300, 2000, 2250, 2980))
        boxes = pr.plan_layout_tiles(img, lay)
        fig = [b for b in boxes if b[1] <= 300]
        self.assertEqual(len(fig), 1)
        self.assertGreaterEqual(fig[0][3], 1780)                      # labels and caption in it
        self.assertEqual(pr.coverage(img, [{"box_px": b} for b in boxes]), 1.0)

    def test_tall_table_is_one_tile(self):
        img = lines(page(), 300, 2250, 200, 3100, pitch=40, height=20)
        boxes = pr.plan_layout_tiles(img, layout(block("Table", 300, 200, 2250, 3090)))
        self.assertEqual(len(boxes), 1)

    def test_text_is_packed_per_column_between_blocks(self):
        img = page()
        cols = []
        for x0, x1 in ((300, 1200), (1350, 2250)):
            for y0 in range(300, 3000, 450):
                lines(img, x0, x1, y0, y0 + 380)
                cols.append(block("Text", x0, y0, x1, y0 + 380))
        boxes = pr.plan_layout_tiles(img, layout(*cols))
        self.assertTrue(all(b[2] <= 1360 or b[0] >= 1190 for b in boxes))      # no tile spans both
        self.assertTrue(all(b[3] - b[1] <= pr.TARGET for b in boxes))
        self.assertLess(len(boxes), len(cols))
        for b in boxes:                                                         # cuts between blocks
            self.assertFalse(any(c["bbox"][1] * pr.DPI / 72 < b[1] < c["bbox"][3] * pr.DPI / 72
                                 and b[0] < c["bbox"][2] * pr.DPI / 72 for c in cols))


    def test_margin_stamp_does_not_swallow_the_page(self):
        img = lines(page(), 300, 2250, 300, 2400)
        img[100:3200, 2450:2470] = 0                                  # a download stamp down the margin
        lay = layout(block("Text", 300, 300, 2250, 1200), block("Text", 300, 1250, 2250, 2400))
        boxes = pr.plan_layout_tiles(img, lay)
        body = [b for b in boxes if b[0] < 2250]
        self.assertTrue(all(b[2] <= 2260 for b in body))
        self.assertEqual(pr.coverage(img, [{"box_px": b} for b in boxes]), 1.0)


class MinerUOutput(unittest.TestCase):
    def test_numbers_join_their_equation_and_reach_prepare(self):
        blocks = [{"type": "paragraph_title", "bbox": [0.1, 0.05, 0.5, 0.07], "content": "2. Model"},
                  {"type": "equation", "bbox": [0.1, 0.10, 0.40, 0.14], "content": "\\[\na = b\n\\]"},
                  {"type": "formula_number", "bbox": [0.44, 0.11, 0.48, 0.13], "content": "(2.14)"},
                  {"type": "equation", "bbox": [0.6, 0.10, 0.88, 0.14], "content": "\\[\nc = d \\tag {7}\n\\]"},
                  {"type": "formula_number", "bbox": [0.91, 0.11, 0.95, 0.13], "content": "(7)"},
                  {"type": "equation", "bbox": [0.1, 0.20, 0.40, 0.24], "content": "\\[ e \\]"}]
        out = dp.attach_numbers([dict(b) for b in blocks])
        self.assertEqual([b["type"] for b in out], ["paragraph_title", "equation", "equation", "equation"])
        self.assertAlmostEqual(out[1]["bbox"][2], 0.48)                 # the number joins the left equation
        eqs = pr.equations({1: dp.markdown(out)})
        self.assertEqual([e["number"] for e in eqs], ["2.14", "7", None])
        self.assertEqual(eqs[0]["latex"], "a = b \\tag{2.14}")
        self.assertEqual([s["level"] for s in pr.sections({1: dp.markdown(out)})], [2])

    def test_numbered_lines_read_as_text_become_equations(self):
        blocks = [{"type": "text", "bbox": [0.13, 0.53, 0.72, 0.55], "content": "\\(u = g\\) on \\(\\partial\\Omega\\), (2.21)"},
                  {"type": "text", "bbox": [0.11, 0.60, 0.72, 0.62], "content": "and we find from \\(u\\) the result we need in (2.5)"},
                  {"type": "equation", "bbox": [0.2, 0.70, 0.7, 0.74], "content": "\\[ v = 0. (2. 4 0) \\]"}]
        eqs = pr.equations({1: dp.markdown(dp.attach_numbers([dict(b) for b in blocks]))})
        self.assertEqual([e["number"] for e in eqs], ["2.21", "2.40"])

    def test_every_mineru_type_maps_to_a_layout_type(self):
        types = {"text", "title", "doc_title", "paragraph_title", "table", "equation", "formula_number", "code",
                 "algorithm", "aside_text", "ref_text", "index", "phonetic", "list_item", "table_caption",
                 "image_caption", "code_caption", "caption", "table_footnote", "image_footnote", "footnote",
                 "header", "footer", "page_number", "page_footnote", "image", "chart", "list", "image_block",
                 "equation_block", "unknown"}
        self.assertEqual(types - set(dp.TYPES), set())


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
