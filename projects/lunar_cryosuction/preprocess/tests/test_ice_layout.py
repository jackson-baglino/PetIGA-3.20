"""Tests for preprocess/build_ice_layout.py (plan type 7 ice layouts).
Run from the project root:  python -m unittest discover -s preprocess/tests -t . -v"""
import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "preprocess"))
import build_ice_layout as bl  # noqa: E402

FLAT = {"seed_bottom": 1, "seed_top": 2, "n_bumps": 10, "R_range": [2.1e-5, 2.1e-5],
        "height_ratio_range": [1e-9, 2e-9], "overlap_range": [0.2, 0.2]}


class Pieces(unittest.TestCase):
    def setUp(self):
        self.spec = {"Lx": 3e-4, "mean_gap": 8e-5, "walls": FLAT}
        self.pore = bl.Pore(self.spec)
        self.rs = bl.Raster(self.pore, 2.5e-7)

    def test_plug_in_a_flat_channel_has_concave_faces_at_theta(self):
        th = 60.0
        ice, r = bl.plug(None, self.rs, 1.0e-4, 2.0e-4, th)
        a = 0.5 * 8e-5
        self.assertAlmostEqual(r * math.cos(math.radians(th)) / a, 1.0, delta=1e-3)
        seg = 2 * a * r - a * math.sqrt(r * r - a * a) - r * r * math.asin(a / r)
        want = 2 * a * 1.0e-4 + 2 * seg           # rectangle + two concave segments
        got = ice.sum() * self.rs.dx ** 2
        self.assertAlmostEqual(got / want, 1.0, delta=5e-3)

    def test_cap_on_a_flat_wall(self):
        R, th = 3e-5, 40.0
        # cell counting misses up to one row (dx) along the 39 um base: use a fine raster
        rs = bl.Raster(self.pore, 5e-8)
        ice = bl.cap(rs, 1.5e-4, "bottom", R, th)
        t = math.radians(th)
        want = R * R * (t - math.sin(t) * math.cos(t))   # circular segment of half-angle theta
        self.assertAlmostEqual(ice.sum() * rs.dx ** 2 / want, 1.0, delta=1e-2)

    def test_phi_file_roundtrip_and_layout(self):
        with tempfile.TemporaryDirectory() as d:
            surf = bl.build_surface(self.pore, 40, 10)
            ice = bl.grain(self.rs, 1.5e-4, "axis", 2e-5, 2e-5)[0]
            rec = bl.node_field(surf, self.rs, ice, 1.5e-6)
            p = Path(d) / "f.bin"
            bl.write_phi(p, rec)
            back = bl.read_phi(p)
            np.testing.assert_array_equal(back, rec)
            self.assertEqual(back.shape, (12, 42, 3))            # (my, mx, 3), x fastest
            np.testing.assert_allclose(back[0, :, 0], surf.control[:, 0, 0])
            self.assertGreater(back[..., 2].max(), 0.99)
            self.assertLess(back[..., 2].min(), 0.01)

    def test_main_writes_files(self):
        spec = dict(self.spec, name="t", theta_ice_init=40,
                    plugs=[{"x0": 2.2e-4, "x1": 3e-4}], grains=[{"x": 1e-4, "y": "axis", "R": 1.5e-5}])
        with tempfile.TemporaryDirectory() as d:
            j = Path(d) / "l.json"
            j.write_text(json.dumps(spec))
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                rc = bl.main([str(j), "--eps", "2e-6", "--outdir", d, "--no-pdf"])
            self.assertEqual(rc, 0, buf.getvalue())
            text = (Path(d) / "t.opts").read_text()
            self.assertIn("-ic_type phi_file", text)
            self.assertIn("# DOF_GRID:", text)
            self.assertTrue((Path(d) / "t_phi.bin").exists())


if __name__ == "__main__":
    unittest.main()
