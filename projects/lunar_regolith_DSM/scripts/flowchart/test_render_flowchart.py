"""Tests for render_flowchart.py: the overlap lint must catch overlaps, and the committed charts must be clean.

    venv_lunar/bin/python -m unittest scripts/flowchart/test_render_flowchart.py
"""
import json
import re
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_flowchart as rf  # noqa: E402

SPECS = sorted((rf.ROOT / "workflows").glob("*/flowchart.json"))


class SegmentRectangle(unittest.TestCase):

    def test_crossing_inside_and_outside(self):
        r = (0, 0, 10, 10)
        self.assertTrue(rf._seg_hits((-5, 5), (15, 5), r))      # passes straight through
        self.assertTrue(rf._seg_hits((2, 2), (3, 3), r))        # wholly inside
        self.assertFalse(rf._seg_hits((-5, 12), (15, 12), r))   # passes above
        self.assertTrue(rf._seg_hits((-1000, 5), (1000, 5), r))  # long segment: sampling could miss this


class Lint(unittest.TestCase):

    def _render(self, spec_path):
        spec = json.loads(spec_path.read_text())
        pos, labels, total = rf.layout(spec)
        import subprocess
        r = subprocess.run(["neato", "-n2", "-Tsvg"], input=rf.to_dot(spec, pos, labels, total),
                           capture_output=True, text=True, check=True)
        svg = re.sub(r"^.*?(<svg)", r"\1", r.stdout, flags=re.S)
        svg, titles = rf.add_lanes(svg, spec, total)
        return spec, rf.place_labels(spec, svg, titles), titles

    def test_committed_charts_are_clean(self):
        for p in SPECS:
            spec, svg, titles = self._render(p)
            self.assertEqual(rf.lint(svg, spec, titles)["issues"], [], p)

    def test_label_moved_onto_its_edge_is_flagged(self):
        spec, svg, titles = self._render(SPECS[0])
        nodes, edges, _ = rf.parse_svg(svg)
        m = re.search(r'<g id="(edge\d+)" class="edge">.*?<text class="elabel" x="[-\d.]+" y="[-\d.]+"', svg, re.S)
        gid = m.group(1)
        pl = next(lines for (g, lines) in edges.values() if g == gid)[0]
        (x1, y1), (x2, y2) = pl[0], pl[-1]
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        bad = re.sub(rf'(<g id="{gid}" class="edge">.*?<text class="elabel") x="[-\d.]+" y="[-\d.]+"',
                     lambda k: f'{k.group(1)} x="{mx - 5:.1f}" y="{my + 3:.1f}"', svg, count=1, flags=re.S)
        types = {i["type"] for i in rf.lint(bad, spec, titles)["issues"]}
        self.assertTrue(types & {"label_over_edge", "label_over_node"}, types)


if __name__ == "__main__":
    unittest.main()
