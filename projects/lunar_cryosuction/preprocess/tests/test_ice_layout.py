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


class Periodic(unittest.TestCase):
    """"periodic": true -- the mesh and phi file for -periodic_x 1."""
    LX, W, WT = 1.2e-4, 4.0e-5, 1.6e-5

    def cd(self, **kw):
        w = dict({"profile": "sinhyp", "W": self.W, "wt": self.WT}, **kw)
        return bl.Pore({"Lx": self.LX, "periodic": True, "walls": w})

    def test_knots_are_petiga_periodic_uniform(self):
        # petigaaxis.c IGAAxisInitUniform, periodic branch, N = 8, p = 2, C = 1 on [0, 1]
        U = bl.periodic_uniform_knots(8, 2)
        self.assertEqual(U[0], -0.25)
        self.assertEqual(U[-1], 1.25)
        np.testing.assert_allclose(np.diff(U), 0.125)
        p, n = 2, len(U) - 2 - 2
        for k in range(2):                        # C = p - 1 = 1: k = 0..C
            self.assertAlmostEqual(U[1 - k], U[p] - U[len(U) - 1 - p] + U[n - k])

    def test_surface_wraps_by_one_period(self):
        surf = bl.build_surface(self.cd(), 40, 10)
        cw = surf.control
        self.assertEqual(cw.shape[0], 42)         # N + p columns
        np.testing.assert_allclose(cw[40:, :, 0] - cw[:2, :, 0], self.LX, rtol=0, atol=1e-18)
        np.testing.assert_array_equal(cw[40:, :, 1], cw[:2, :, 1])
        # the mapping and the wall slope agree at u = 0 and u = 1
        v = np.array([0.0, 1.0])
        np.testing.assert_allclose(surf(1.0, v)[:, 1], surf(0.0, v)[:, 1], rtol=0, atol=1e-18)
        g0, g1 = surf.gradient(u=np.array([0.0]), v=v), surf.gradient(u=np.array([1.0]), v=v)
        np.testing.assert_allclose(g1[0, :, 1, 0] / g1[0, :, 0, 0], g0[0, :, 1, 0] / g0[0, :, 0, 0],
                                   rtol=0, atol=1e-12)

    def test_profiles(self):
        for kind in ("cos", "sinhyp"):
            g = bl.gap_profile(np.array([0.0, 0.5 * self.LX, self.LX]), self.LX, self.W, self.WT, kind)
            np.testing.assert_allclose(g, [self.W, self.WT, self.W], rtol=1e-12)
            h = 1e-9                               # slope continuous (and zero) at the seam
            s0 = (bl.gap_profile(h, self.LX, self.W, self.WT, kind) - bl.gap_profile(-h, self.LX, self.W, self.WT, kind))
            self.assertLess(abs(s0) / (2 * h), 1e-6)
        # sinhyp is the hyperbola sqrt(wt^2 + c^2 xi^2) near the throat
        c = math.pi * math.sqrt(self.W ** 2 - self.WT ** 2) / self.LX
        xi = 2e-6
        g = bl.gap_profile(0.5 * self.LX + xi, self.LX, self.W, self.WT, "sinhyp")
        self.assertAlmostEqual(g / math.sqrt(self.WT ** 2 + (c * xi) ** 2), 1.0, delta=1e-3)
        # x_throat moves the throat
        g = bl.gap_profile(0.0, self.LX, self.W, self.WT, "cos", xt=0.0)
        self.assertAlmostEqual(g, self.WT)

    def test_raster_has_whole_cells_per_period(self):
        rs = bl.Raster(self.cd(), 3.7e-7)
        self.assertEqual(len(rs.x), 3 * rs.M)
        self.assertAlmostEqual(rs.M * rs.dx, self.LX, delta=1e-18)
        self.assertAlmostEqual(rs.x0, -self.LX)

    def test_plug_across_the_seam_is_the_shifted_plug(self):
        pore = self.cd(x_throat=0.0)
        rs = bl.Raster(pore, 3.75e-7)
        a = rs.fold(bl.plug(None, rs, -2.5e-5, -1.0e-5, 60.0)[0])     # crosses x = 0
        pore_b = self.cd()
        rs_b = bl.Raster(pore_b, 3.75e-7)
        b = rs_b.fold(bl.plug(None, rs_b, 3.5e-5, 5.0e-5, 60.0)[0])
        np.testing.assert_array_equal(np.roll(rs.middle(a), rs.M // 2, axis=0), rs_b.middle(b))
        lab, n = bl.label_bodies(a, rs)
        self.assertEqual(n, 1)                   # one body, though it crosses the seam

    def test_wrapped_columns_get_the_same_phi(self):
        pore = bl.Pore({"Lx": self.LX, "periodic": True, "mean_gap": 4e-5, "walls": {"flat": True}})
        rs = bl.Raster(pore, 3.75e-7)
        ice = rs.fold(bl.plug(None, rs, -1e-5, 1e-5, 60.0)[0])
        surf = bl.build_surface(pore, 114, 38)
        rec = bl.node_field(surf, rs, ice, 1.5e-6)
        self.assertEqual(rec.shape[1], 116)
        np.testing.assert_allclose(rec[:, 114:, 2], rec[:, :2, 2], rtol=0, atol=1e-12)
        self.assertGreater(rec[:, 0, 2].min(), 0.99)     # x = -h/2 is inside the seam lens

    def test_main_writes_periodic_opts(self):
        spec = {"name": "tp", "Lx": self.LX, "periodic": True, "theta_ice_init": 60,
                "walls": {"profile": "cos", "W": self.W, "wt": self.WT},
                "plugs": [{"x0": -2.5e-5, "x1": -1.0e-5}]}
        with tempfile.TemporaryDirectory() as d:
            j = Path(d) / "l.json"
            j.write_text(json.dumps(spec))
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                rc = bl.main([str(j), "--eps", "1.5e-6", "--outdir", d, "--no-pdf"])
            self.assertEqual(rc, 0, buf.getvalue())
            text = (Path(d) / "tp.opts").read_text()
            self.assertIn("-periodic_x 1", text)
            nx = math.ceil(self.LX * math.sqrt(2) / 1.5e-6)
            ny = math.ceil(self.W * math.sqrt(2) / 1.5e-6)
            self.assertIn(f"# DOF_GRID: {nx} {ny + 2}", text)
            self.assertEqual(bl.read_phi(Path(d) / "tp_phi.bin").shape, (ny + 2, nx, 3))

    def test_periodic_refuses_bumps(self):
        with self.assertRaises(AssertionError):
            bl.Pore({"Lx": 3e-4, "mean_gap": 8e-5, "periodic": True, "walls": FLAT})


if __name__ == "__main__":
    unittest.main()
