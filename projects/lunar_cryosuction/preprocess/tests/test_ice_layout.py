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

    def test_wall_points_inside_a_plug_are_ice(self):
        # regression (2026-11-09): wall control points that fell in a raster cell just
        # outside the pore were read as air, leaving isolated phi ~ 0.6 dips on the wall
        rough = dict(FLAT, n_bumps=10, R_range=[2.0e-5, 2.2e-5], height_ratio_range=[0.3, 0.5])
        pore = bl.Pore({"Lx": 3e-4, "mean_gap": 8e-5, "walls": rough})
        rs = bl.Raster(pore, 2.5e-7)
        ice, _ = bl.plug(None, rs, 1.0e-4, 2.0e-4, 60.0)
        surf = bl.build_surface(pore, 300, 80)
        rec = bl.node_field(surf, rs, ice, 1.0e-6)
        x = rec[0, :, 0]
        bulk = (x > 1.3e-4) & (x < 1.7e-4)             # well inside the plug on both walls
        for row in (0, -1):
            self.assertGreater(rec[row, bulk, 2].min(), 0.999)

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


class NewFeatures(unittest.TestCase):
    """rinds, necks, bubbles and smoothed walls (added 2026-11-13 for the type-7 v2 layouts)."""

    def setUp(self):
        self.pore = bl.Pore({"Lx": 3e-4, "mean_gap": 8e-5, "walls": FLAT})
        self.rs = bl.Raster(self.pore, 1.25e-7)

    def test_rind_on_a_flat_wall_is_a_strip_with_round_ends(self):
        t, x0, x1 = 5e-6, 1.0e-4, 2.0e-4
        ice = bl.rind(self.rs, t, "bottom", x0, x1)
        want = t * (x1 - x0) + 2 * (math.pi * t * t / 4)   # strip + a quarter disc at each end
        self.assertAlmostEqual(ice.sum() * self.rs.dx ** 2 / want, 1.0, delta=2e-2)
        # uniform thickness: the film reaches t above the wall, not more, and only on the bottom wall
        i = int(round(1.5e-4 / self.rs.dx))
        col = np.flatnonzero(ice[i]) * self.rs.dx - self.rs.yb[i]
        self.assertAlmostEqual(col.max() / t, 1.0, delta=0.05)
        self.assertFalse(ice[:, self.rs.y > 0.5 * self.pore.H].any())

    def test_rind_on_both_walls(self):
        ice = bl.rind(self.rs, 4e-6, "both", 1.0e-4, 2.0e-4)
        top = ice[:, self.rs.y > 0.5 * self.pore.H].sum()
        bot = ice[:, self.rs.y <= 0.5 * self.pore.H].sum()
        self.assertAlmostEqual(top / bot, 1.0, delta=0.03)

    def test_neck_on_a_flat_wall_is_the_cap(self):
        R, th = 2e-5, 60.0
        ice, angs = bl.neck(self.rs, 1.5e-4, "bottom", R, th)
        for a in angs:
            self.assertAlmostEqual(a, th, delta=1.5)
        t = math.radians(th)
        want = R * R * (t - math.sin(t) * math.cos(t))
        self.assertAlmostEqual(ice.sum() * self.rs.dx ** 2 / want, 1.0, delta=2e-2)
        top, _ = bl.neck(self.rs, 1.5e-4, "top", R, th)
        self.assertAlmostEqual(top.sum() / ice.sum(), 1.0, delta=2e-2)

    def test_neck_in_a_valley_meets_both_flanks_near_theta(self):
        rough = dict(FLAT, n_bumps=7, R_range=[3.0e-5, 3.5e-5], height_ratio_range=[0.4, 0.5],
                     overlap_range=[0.25, 0.3], smooth=4e-6)
        pore = bl.Pore({"Lx": 3e-4, "mean_gap": 8e-5, "walls": rough})
        rs = bl.Raster(pore, 2.5e-7)
        xv = bl.wall_valleys(rs, "bottom", 5e-5, 2.5e-4, 2e-6)
        self.assertGreater(len(xv), 0)
        _, angs = bl.neck(rs, xv[0], "bottom", 1.5e-5, 60.0)
        self.assertAlmostEqual(np.mean(angs), 60.0, delta=1.0)

    def test_smoothed_walls_keep_mean_gap_and_throat_gap(self):
        w = dict(FLAT, n_bumps=12, R_range=[2.0e-5, 3.0e-5], height_ratio_range=[0.5, 0.7],
                 overlap_range=[0.2, 0.3])
        spec = {"Lx": 4e-4, "mean_gap": 1e-4, "walls": w, "throat": {"x": 2e-4, "R": 8e-5, "gap": 4e-5}}
        raw = bl.Pore(spec)
        sm = bl.Pore(dict(spec, walls=dict(w, smooth=5e-6)))
        x = np.linspace(0, 4e-4, 40001)
        for p in (raw, sm):
            yb, yt = p.walls(x)
            win = np.abs(x - 2e-4) <= 8e-5
            self.assertAlmostEqual((yt - yb)[win].min() / 4e-5, 1.0, delta=1e-3)

        def min_valley_radius(p):
            yb, _ = p.walls(x)
            d1 = np.gradient(yb, x)
            k = np.gradient(d1, x) / (1 + d1 ** 2) ** 1.5
            return 1.0 / k.max()
        self.assertGreater(min_valley_radius(sm), 3 * min_valley_radius(raw))
        # the mean gap (without the throat) is kept: compare against a smoothed pore with no throat
        nt = bl.Pore({"Lx": 4e-4, "mean_gap": 1e-4, "walls": dict(w, smooth=5e-6)})
        yb, yt = nt.walls(x)
        self.assertAlmostEqual((yt - yb).mean() / 1e-4, 1.0, delta=1e-3)

    def test_main_with_rind_neck_plug_and_closed_bubble(self):
        spec = {"Lx": 3e-4, "mean_gap": 8e-5, "walls": FLAT, "name": "t2", "theta_ice_init": 60,
                "rinds": [{"t": 8e-6, "wall": "both", "x0": 2e-5, "x1": 3e-4}],
                "plugs": [{"x0": 1.2e-4, "x1": 2.2e-4}],
                "bubbles": [{"x": 1.7e-4, "y": "axis", "R": 2e-5}]}
        with tempfile.TemporaryDirectory() as d:
            j = Path(d) / "l.json"
            j.write_text(json.dumps(spec))
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                rc = bl.main([str(j), "--eps", "1e-6", "--outdir", d, "--no-pdf"])
            out = buf.getvalue()
            self.assertEqual(rc, 0, out)
            self.assertIn("bubble at x 170 um: closed", out)
            rec = bl.read_phi(Path(d) / "t2_phi.bin")
            # the bubble center is gas, the plug next to it is ice
            k = np.argmin(np.abs(rec[0, :, 0] - 1.7e-4))
            mid = rec.shape[0] // 2
            self.assertLess(rec[mid, k, 2], 0.01)
            k2 = np.argmin(np.abs(rec[0, :, 0] - 1.35e-4))
            self.assertGreater(rec[mid, k2, 2], 0.99)

    def test_open_bubble_fails_the_check(self):
        # a bubble in open gas (no ice around it) is reported and the build returns 2
        spec = {"Lx": 3e-4, "mean_gap": 8e-5, "walls": FLAT, "name": "t3", "theta_ice_init": 60,
                "plugs": [{"x0": 2.2e-4, "x1": 3e-4}],
                "bubbles": [{"x": 1.0e-4, "y": "axis", "R": 2e-5}]}
        with tempfile.TemporaryDirectory() as d:
            j = Path(d) / "l.json"
            j.write_text(json.dumps(spec))
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                rc = bl.main([str(j), "--eps", "1e-6", "--outdir", d, "--no-pdf"])
            self.assertEqual(rc, 2)
            self.assertIn("NOT CLOSED", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
