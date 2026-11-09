"""Tests for preprocess/build_geometry_dogbone.py (plan type 6).
Run from the project root:  python -m unittest discover -s preprocess/tests -t . -v"""
import contextlib
import io
import re
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "preprocess"))
import build_geometry_dogbone as db  # noqa: E402

LX, W, WC, X0, X1, LS = 1.2e-3, 1.5e-4, 4.0e-5, 4.5e-4, 7.5e-4, 1.0e-4


class Profile(unittest.TestCase):
    def test_psi_is_an_exact_step(self):
        self.assertTrue(np.all(db.psi(np.array([-1.0, 0.0])) == 0.0))
        self.assertTrue(np.all(db.psi(np.array([1.0, 2.0])) == 1.0))
        self.assertAlmostEqual(float(db.psi(np.array([0.5]))[0]), 0.5, places=15)
        t = np.linspace(0, 1, 1001)
        self.assertTrue(np.all(np.diff(db.psi(t)) >= 0.0))

    def test_neck_straight_and_ends_flat(self):
        yb, yt = db.walls(np.linspace(X0, X1, 501), W, WC, X0, X1, LS)
        np.testing.assert_allclose(yt - yb, WC, rtol=0, atol=1e-18)
        for xs in (np.linspace(0, X0 - LS, 501), np.linspace(X1 + LS, LX, 501)):
            yb, yt = db.walls(xs, W, WC, X0, X1, LS)
            self.assertTrue(np.all(yb == 0.0) and np.all(yt == W))

    def test_symmetric_and_steepest_slope(self):
        x = np.linspace(0, LX, 400001)
        yb, yt = db.walls(x, W, WC, X0, X1, LS)
        np.testing.assert_allclose(yb + yt, W, rtol=0, atol=1e-18)
        s = np.abs(np.gradient(yb, x)).max()
        self.assertAlmostEqual(s, (W - WC) / LS, delta=1e-3 * s)   # psi'(1/2) = 2

    def test_open_area(self):
        x = np.linspace(0, LX, 200001)
        yb, yt = db.walls(x, W, WC, X0, X1, LS)
        area = float(np.sum(0.5 * ((yt - yb)[1:] + (yt - yb)[:-1]) * np.diff(x)))
        # psi(t) + psi(1 - t) = 1, so each shoulder removes exactly half its rectangle
        want = LX * W - (W - WC) * (X1 - X0 + LS)
        self.assertAlmostEqual(area / want, 1.0, places=8)


class Builder(unittest.TestCase):
    def _run(self, *args):
        with tempfile.TemporaryDirectory() as d:
            out, opts = Path(d) / "m.dat", Path(d) / "g.opts"
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = db.main(["--out", str(out), "--opts", str(opts), "--png", str(Path(d) / "p.pdf"), *args])
            return rc, buf.getvalue(), opts.read_text() if opts.exists() else ""

    def test_opts_carries_lens_walls_and_profile(self):
        rc, log, text = self._run("--wc", "4e-5", "--Lx", "4e-4", "--W", "5e-5", "--neck-len", "1e-4",
                                  "--shoulder", "4e-5", "--lens-len", "4e-5", "--eps", "1.5e-6")
        self.assertEqual(rc, 0, log)
        self.assertIn("-wall_faces y0,y1", text)
        self.assertRegex(text, r"# DOGBONE: W 5\.000000e-05 wc 4\.000000e-05 x_n0 1\.500000e-04 x_n1 2\.500000e-04")
        r = float(re.search(r"-bridge_rL\s+(\S+)", text).group(1))
        cxl = float(re.search(r"-bridge_cxL\s+(\S+)", text).group(1))
        cxr = float(re.search(r"-bridge_cxR\s+(\S+)", text).group(1))
        # lens in the middle of the right flat section: 290-400 um, centered 345 um, 40 um long
        self.assertAlmostEqual((cxl - r) * 1e6, 325.0, places=3)
        self.assertAlmostEqual((cxr + r) * 1e6, 365.0, places=3)
        self.assertIn("mesh wall vs profile", log)

    def test_lens_on_a_shoulder_is_refused(self):
        with self.assertRaises(SystemExit):
            self._run("--wc", "4e-5", "--lens-x", "8.0e-4")

    def test_lens_too_close_to_the_end_is_refused(self):
        with self.assertRaises(AssertionError):
            self._run("--wc", "4e-5", "--lens-x", "1.12e-3")


if __name__ == "__main__":
    unittest.main()
