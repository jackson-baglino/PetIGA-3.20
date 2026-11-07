"""Unit tests for the core helpers in postprocess/pplib.py.

Scope: load_ssa, grain_radii, rho_vs, supersaturation, read_vts, step_of,
step_times, auto_time_unit, in_time_unit, crossings, refine_tanh,
contour_points, circle_radius. Every check is against an input whose answer is
known exactly (an analytic tanh profile, exact circle points, the solver's own
saturation formula), or against a second, independent output of the same real
run.

Run from the project root:
    venv_cryosuction/bin/python -m unittest discover -s postprocess/tests -t . -p 'test_pplib_core.py'

Tests marked expectedFailure document a bug; the comment above each starts
with "# BUG:" and names the line.
"""
from __future__ import annotations

import math
import os
import re
import sys
import tempfile
import unittest
import warnings

import numpy as np

from postprocess.tests import fixtures as fx

sys.path.insert(0, fx.POSTPROCESS)
sys.path.insert(0, os.path.join(fx.PROJECT, "preprocess"))
import pplib  # noqa: E402
from mesh_params import physics  # noqa: E402


def tanh_profile(x, x0, eps):
    """phi = (1 + tanh((x - x0)/(2 eps)))/2: crosses 0.5 exactly at x0."""
    return fx.profile(np.asarray(x) - x0, eps)


def read_text(path):
    with open(path) as fh:
        return fh.read()


def load_quiet(path, **kw):
    """load_ssa without numpy's 'line has N columns' warnings in the test log."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pplib.load_ssa(path, **kw)


# ---------------------------------------------------------------------------
# load_ssa
# ---------------------------------------------------------------------------
class TestLoadSSA(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run = fx.make_run(self._tmp.name)
        self.path = os.path.join(self.run, "SSA_evo.dat")

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, text):
        with open(self.path, "w") as fh:
            fh.write(text)

    def test_eight_columns_roundtrip_by_dir_and_file(self):
        t = np.array([0.0, 1.5, 4.0, 10.0])
        ice = np.array([1e-8, 1.1e-8, 1.2e-8, 1.3e-8])
        rhov = np.array([2e-11, 2.1e-11, 2.2e-11, 2.3e-11])
        fx.write_ssa(self.run, t, ice, rhov=rhov, air=np.full(4, 3e-8), ssa=np.full(4, 5e-5))
        for arg in (self.run, self.path):
            d = pplib.load_ssa(arg)
            self.assertEqual(d.shape, (4, 8))
            np.testing.assert_allclose(d[:, pplib.TIME], t, rtol=1e-6)
            np.testing.assert_allclose(d[:, pplib.ICE], ice, rtol=1e-6)
            np.testing.assert_allclose(d[:, pplib.RHOV], rhov, rtol=1e-6)
            np.testing.assert_array_equal(d[:, pplib.STEP], np.arange(4))
            np.testing.assert_allclose(d[:, pplib.DT], [1.5, 2.5, 6.0, 6.0], rtol=1e-6)
            np.testing.assert_allclose(d[:, pplib.MASS], 919.0 * ice + rhov, rtol=1e-6)

    def test_missing_empty_header_only(self):
        self.assertIsNone(pplib.load_ssa(os.path.join(self.run, "nope.dat")))
        self.assertIsNone(pplib.load_ssa(self._tmp.name))     # dir without the file
        self.write("")
        self.assertIsNone(load_quiet(self.path))
        self.write("# sub_interf ice t step\n# nothing yet\n")
        self.assertIsNone(load_quiet(self.path))

    def test_one_row(self):
        fx.write_ssa(self.run, [0.0], [2e-8])
        d = pplib.load_ssa(self.run)
        self.assertEqual(d.shape, (1, 8))
        self.assertEqual(d[0, pplib.ICE], 2e-8)

    def test_legacy_four_columns_padded_and_all_rows_kept(self):
        fx.write_ssa(self.run, [0.0, 1.0, 2.0], [1.0, 2.0, 3.0], ncol=4)
        d = pplib.load_ssa(self.run)
        self.assertEqual(d.shape, (3, 5))
        self.assertTrue(np.isnan(d[:, pplib.DT]).all())
        np.testing.assert_array_equal(d[:, pplib.ICE], [1.0, 2.0, 3.0])

    def test_legacy_one_row_four_columns(self):
        fx.write_ssa(self.run, [5.0], [1.0], ncol=4)
        d = pplib.load_ssa(self.run)
        self.assertEqual(d.shape, (1, 5))
        self.assertTrue(np.isnan(d[0, pplib.DT]))

    def test_legacy_five_columns(self):
        fx.write_ssa(self.run, [0.0, 1.0, 3.0], [1.0, 2.0, 3.0], ncol=5)
        d = pplib.load_ssa(self.run)
        self.assertEqual(d.shape, (3, 5))
        np.testing.assert_allclose(d[:, pplib.DT], [1.0, 2.0, 2.0])

    def test_min_cols_rejects_short_files(self):
        fx.write_ssa(self.run, [0.0, 1.0], [1.0, 2.0], ncol=4)
        self.assertIsNone(pplib.load_ssa(self.run, min_cols=5))
        fx.write_ssa(self.run, [0.0, 1.0], [1.0, 2.0], ncol=3)
        self.assertIsNone(pplib.load_ssa(self.run))

    def test_duplicate_steps_keep_last(self):
        # A retried step written twice, adjacent: the second write wins.
        fx.write_ssa(self.run, [0.0, 1.0, 1.0, 2.0], [1.0, 2.0, 2.5, 3.0], steps=[0, 1, 1, 2])
        d = pplib.load_ssa(self.run)
        np.testing.assert_array_equal(d[:, pplib.STEP], [0, 1, 2])
        np.testing.assert_array_equal(d[:, pplib.ICE], [1.0, 2.5, 3.0])

    def test_restart_overlap_keeps_newer_segment(self):
        # Steps 0-5, then a restart from step 3 that runs on to 7.
        steps = [0, 1, 2, 3, 4, 5, 3, 4, 5, 6, 7]
        ice = [0, 1, 2, 3, 4, 5, 30, 40, 50, 60, 70]
        fx.write_ssa(self.run, np.array(steps, float), np.array(ice, float), steps=steps)
        d = pplib.load_ssa(self.run)
        np.testing.assert_array_equal(d[:, pplib.STEP], np.arange(8))
        np.testing.assert_array_equal(d[:, pplib.ICE], [0, 1, 2, 30, 40, 50, 60, 70])

    # BUG: postprocess/pplib.py:83-85 keeps the LAST occurrence of each step in
    # FILE order, so a restarted segment that ends before the old one did
    # (steps 0..5 then 3..4) comes back as steps 0,1,2,5,3,4: time runs
    # backwards and every caller plots a zigzag. Right: rows ordered by step
    # (and ideally old rows past the restart point dropped). Low practical
    # impact: the solver truncates SSA_evo.dat on open, so this needs a
    # concatenated file.
    @unittest.expectedFailure
    def test_restart_that_rewinds_is_monotonic(self):
        steps = [0, 1, 2, 3, 4, 5, 3, 4]
        fx.write_ssa(self.run, np.array(steps, float) * 10.0, np.ones(8), steps=steps)
        d = pplib.load_ssa(self.run)
        self.assertTrue(np.all(np.diff(d[:, pplib.STEP]) > 0), d[:, pplib.STEP])

    def test_nan_rows(self):
        # NaN (and C's "-nan") in t drops the row; NaN in dt keeps it.
        self.write("1 2 0 0 1 0 0 0\n"
                   "1 2 nan 1 1 0 0 0\n"
                   "1 2 2 2 nan 0 0 0\n"
                   "1 2 3 3 -nan 0 0 0\n"
                   "1 2 -nan 4 1 0 0 0\n")
        d = pplib.load_ssa(self.path)
        np.testing.assert_array_equal(d[:, pplib.STEP], [0, 2, 3])
        self.assertTrue(np.isnan(d[1:, pplib.DT]).all())

    def test_all_rows_nan_returns_none(self):
        self.write("nan nan nan nan nan\n")
        self.assertIsNone(pplib.load_ssa(self.path))

    def test_truncated_last_line_is_dropped(self):
        # A run killed mid-write leaves a partial last row.
        self.write("1 2 0 0 1 0 0 0\n1 2 1 1 1 0 0 0\n1 2 2 2\n")
        d = load_quiet(self.path)
        np.testing.assert_array_equal(d[:, pplib.STEP], [0, 1])

    def test_comment_lines_ignored(self):
        self.write("# header\n1 2 0 0 1 0 0 0\n# mid comment\n1 2 1 1 1 0 0 0\n")
        self.assertEqual(pplib.load_ssa(self.path).shape, (2, 8))


# ---------------------------------------------------------------------------
# grain_radii
# ---------------------------------------------------------------------------
class TestGrainRadii(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def grains(self, text, name="grains.dat", where=None):
        path = os.path.join(where or self.tmp, name)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_staged_copy_with_header_comments_and_both_layouts(self):
        run = fx.make_run(self.tmp, geometry={"-grains_file": "inputs/somewhere/grains.dat"})
        self.grains("# packing\n2e-4 1e-4\n1e-5 2e-5 3e-6\n\n4e-5 5e-5 6e-5 7e-6\nbad row here\n",
                    where=run)
        np.testing.assert_array_equal(pplib.grain_radii(run), [3e-6, 7e-6])

    def test_absolute_path_not_staged(self):
        path = self.grains("0 0 1.5e-5\n1 1 2.5e-5\n")
        run = fx.make_run(os.path.join(self.tmp, "r"), geometry={"-grains_file": path})
        np.testing.assert_array_equal(pplib.grain_radii(run), [1.5e-5, 2.5e-5])

    def test_no_option_or_missing_file_or_no_rows(self):
        run = fx.make_run(self.tmp)
        self.assertIsNone(pplib.grain_radii(run))
        run2 = fx.make_run(self.tmp, geom="g2", geometry={"-grains_file": "/no/such/grains.dat"})
        self.assertIsNone(pplib.grain_radii(run2))
        path = self.grains("2e-4 1e-4\n# only a header\n", name="empty.dat")
        run3 = fx.make_run(self.tmp, geom="g3", geometry={"-grains_file": path})
        self.assertIsNone(pplib.grain_radii(run3))

    # BUG: postprocess/pplib.py:268-269 resolves a relative -grains_file
    # against the CURRENT WORKING DIRECTORY, while the docstring (l. 253) says
    # the path is relative to the PROJECT root. The same run gives radii when
    # plotted from the project root and None (so plot_ssa silently falls back
    # to the t = 0 normalization) from anywhere else. Right: also try the path
    # under the project root (e.g. $PETIGA_DIR/projects/lunar_cryosuction or
    # the directory above pplib.py).
    @unittest.expectedFailure
    def test_project_relative_path_independent_of_cwd(self):
        path = self.grains("0 0 4e-6\n")
        rel = os.path.relpath(path, fx.PROJECT)      # "relative to the PROJECT root"
        run = fx.make_run(os.path.join(self.tmp, "r"), geometry={"-grains_file": rel})
        old = os.getcwd()
        try:
            os.chdir(fx.PROJECT)
            from_root = pplib.grain_radii(run)
            os.chdir(self.tmp)
            from_elsewhere = pplib.grain_radii(run)
        finally:
            os.chdir(old)
        np.testing.assert_array_equal(from_root, [4e-6])     # passes
        self.assertIsNotNone(from_elsewhere)                  # fails


# ---------------------------------------------------------------------------
# rho_vs, supersaturation
# ---------------------------------------------------------------------------
class TestRhoVs(unittest.TestCase):
    TEMPS = (-173.0, -120.0, -100.0, -80.0, -50.0, -20.0, -5.0, -1e-9, 0.0)

    def test_equals_mesh_params_physics(self):
        for T in self.TEMPS:
            with self.subTest(T=T):
                self.assertEqual(float(pplib.rho_vs(T)), physics.rho_vs_sat_kg_per_m3(T))

    def test_vectorized_matches_scalar(self):
        T = np.array(self.TEMPS)
        np.testing.assert_array_equal(pplib.rho_vs(T), [pplib.rho_vs(t) for t in T])

    def test_coefficients_match_c_source(self):
        src = read_text(os.path.join(fx.PROJECT, "src", "material_properties.c"))
        body = src[src.index("void RhoVS_I"):]
        body = body[:body.index("}")]
        got = {k: float(v) for k, v in re.findall(r"\b(K[0-5])\s*=\s*([-0-9.eE+]+)", body)}
        self.assertEqual(tuple(got[f"K{i}"] for i in range(6)), pplib._K)
        self.assertIn("Patm = 101325.0", body)
        self.assertIn("bb = 0.62", body)
        main = read_text(os.path.join(fx.PROJECT, "src", "cryosuction_main.c"))
        self.assertRegex(main, r"user\.rho_air\s*=\s*1\.341;")
        self.assertEqual(pplib.RHO_AIR_DEFAULT, 1.341)

    def test_saturation_pressure_against_murphy_koop(self):
        # Independent check of the coefficients: back out P_vs from rho_vs and
        # compare with Murphy and Koop (2005, QJRMS 131, Eq. 7), valid > 110 K.
        def p_from_rho(T):
            r = float(pplib.rho_vs(T))
            return r * 101325.0 / (1.341 * 0.62 + r)

        def murphy_koop(T_C):
            T = T_C + 273.15
            return math.exp(9.550426 - 5723.265 / T + 3.53068 * math.log(T) - 0.00728332 * T)

        for T in (-80.0, -50.0, -20.0, -5.0, 0.0):
            with self.subTest(T=T):
                self.assertAlmostEqual(p_from_rho(T) / murphy_koop(T), 1.0, delta=2e-3)

    def test_magnitudes(self):
        # -20 C: 8.4869e-4 kg/m^3 with the solver's formula (the module is the
        # reference; 8.9e-4 quoted in the task is 5 % high). Near 0 C: ~5.05e-3.
        self.assertAlmostEqual(float(pplib.rho_vs(-20.0)) / 8.4869e-4, 1.0, delta=1e-4)
        self.assertAlmostEqual(float(pplib.rho_vs(0.0)) / 5.0472e-3, 1.0, delta=1e-4)

    def test_lunar_cold_is_finite_positive_and_monotonic(self):
        T = np.linspace(-173.0, 0.0, 200)
        r = pplib.rho_vs(T)
        self.assertTrue(np.all(np.isfinite(r)) and np.all(r > 0))
        self.assertTrue(np.all(np.diff(r) > 0))
        self.assertLess(float(pplib.rho_vs(-173.0)), 1e-18)

    def test_kelvin_by_mistake_is_nonsense(self):
        # Passing Kelvin is not caught: P_vs > Patm gives a NEGATIVE density,
        # which supersaturation() then maps silently to 0.
        self.assertLess(float(pplib.rho_vs(253.15)), 0.0)
        self.assertEqual(float(pplib.supersaturation(1.0, 253.15)), 0.0)

    def test_supersaturation(self):
        for T in (-173.0, -20.0, 0.0):
            rvs = pplib.rho_vs(T)
            with self.subTest(T=T):
                self.assertAlmostEqual(float(pplib.supersaturation(rvs, T)), 0.0, delta=1e-14)
                self.assertAlmostEqual(float(pplib.supersaturation(1.1 * rvs, T)), 0.1, delta=1e-12)
                self.assertAlmostEqual(float(pplib.supersaturation(0.0, T)), -1.0, delta=1e-14)
        rv = np.array([[1e-4, 8.4869e-4], [1e-3, 0.0]])
        s = pplib.supersaturation(rv, np.full((2, 2), -20.0))
        self.assertEqual(s.shape, (2, 2))
        np.testing.assert_allclose(s, rv / pplib.rho_vs(-20.0) - 1.0, rtol=1e-12, atol=1e-15)

    def test_rho_air_scales_linearly(self):
        self.assertAlmostEqual(float(pplib.rho_vs(-20.0, rho_air=2.682) / pplib.rho_vs(-20.0)), 2.0,
                               places=12)


# ---------------------------------------------------------------------------
# read_vts, step_of
# ---------------------------------------------------------------------------
class TestReadVts(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.fn = os.path.join(self._tmp.name, "solV_00042.vts")
        self.X, self.Y = fx.grid(3e-4, 1e-4, 7, 4)          # nx != ny to catch a transpose
        self.F = {"IcePhase": self.X * 1e3 + self.Y * 1e5,
                  "Temperature": -20.0 + 0.0 * self.X,
                  "VaporDensity": np.arange(28.0).reshape(4, 7)}
        fx.write_vts(self.fn, self.X, self.Y, self.F)

    def tearDown(self):
        self._tmp.cleanup()

    def test_roundtrip_all_fields(self):
        f, X, Y = pplib.read_vts(self.fn)
        self.assertEqual(set(f), set(self.F))
        np.testing.assert_array_equal(X, self.X)
        np.testing.assert_array_equal(Y, self.Y)
        for k in self.F:
            np.testing.assert_array_equal(f[k], self.F[k])
            self.assertEqual(f[k].shape, (4, 7))

    def test_want_subset_and_missing_name(self):
        f, X, _ = pplib.read_vts(self.fn, want=("IcePhase", "NoSuchArray"))
        self.assertEqual(set(f), {"IcePhase"})            # missing name silently absent
        np.testing.assert_array_equal(X, self.X)
        f, _, _ = pplib.read_vts(self.fn, want=("NoSuchArray",))
        self.assertEqual(f, {})

    def test_3d_extent_fails_loudly(self):
        fn = os.path.join(self._tmp.name, "solV_3d.vts")
        fx.write_vts(fn, self.X, self.Y, self.F)
        txt = read_text(fn).replace('WholeExtent="0 6 0 3 0 0"', 'WholeExtent="0 6 0 1 0 1"')
        with open(fn, "w") as fh:
            fh.write(txt)
        # 7 x 2 x 2 = 28 points: read_vts ignores z and reshapes to (2, 7, 3) -> error.
        with self.assertRaises(ValueError):
            pplib.read_vts(fn)

    def test_step_of(self):
        self.assertEqual(pplib.step_of(self.fn), 42)
        self.assertEqual(pplib.step_of("solV_00000.vts"), 0)
        self.assertEqual(pplib.step_of("/a/b/vtkOut/solV_123456.vts"), 123456)
        with self.assertRaises(AttributeError):
            pplib.step_of("sol_00042.dat")


# ---------------------------------------------------------------------------
# step_times
# ---------------------------------------------------------------------------
class TestStepTimes(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run = fx.make_run(self._tmp.name)
        self.path = os.path.join(self.run, "outp.txt")

    def tearDown(self):
        self._tmp.cleanup()

    def test_repeated_headers(self):
        steps = np.arange(50)
        t = 0.37 * steps ** 1.5
        fx.write_outp(self.run, steps, t, every_header=7)
        for arg in (self.run, self.path):
            m = pplib.step_times(arg)
            self.assertEqual(sorted(m), list(steps))
            np.testing.assert_allclose([m[s] for s in steps], t, rtol=1e-5, atol=1e-300)

    def test_missing_file(self):
        self.assertEqual(pplib.step_times(os.path.join(self.run, "nope.txt")), {})
        self.assertEqual(pplib.step_times(self._tmp.name), {})

    def test_other_tables_are_ignored(self):
        lines = [
            "    STEP |     TIME [s] |    DT [s] |    TOT_ICE |    TOT_AIR |      TEMP |  TOT_RHOV "
            "| I-A INTERF | TOTAL_MASS",
            # solver's real data-row format (" %5d | %12.5e | ...")
            "     0 |  0.00000e+00 | 1.000e-04 |  2.156e-08 |  2.531e-08 | -9.375e-07 | 2.148e-11 "
            "|  1.839e-11 |  1.982e-05",
            # percentage row: 8 pipes, no step number
            "       |              |           |    -0.001% |    +0.001% |           |   +0.000% "
            "|            |    -0.000%",
            # SNES table row: 11 pipes, starts with an integer
            "     0 |  5.3957e-04 |   1.108e-15 | 1.000e+00 | 0.000e+00 |   9.564e-13 | 1.000e+00 "
            "| 0.000e+00 |   5.396e-04 | 1.000e+00 | 0.000e+00",
            "     1 |  1.6144e-12 |   1.396e-23 | 1.260e-08 | 1.004e-08 |   1.614e-12 | 1.688e+00 "
            "| 3.334e-09 |   3.447e-17 | 6.388e-14 | 4.488e-04",
            "  BOUNDS: phi_ice [0.0000e+00, 1.0000e+00]  phi_air [-4.4409e-16, 1.0000e+00]",
            "     1 |  1.00000e-04 | 1.334e-04 |  2.156e-08 |  2.531e-08 | -9.375e-07 | 2.148e-11 "
            "|  1.839e-11 |  1.982e-05",
            "  tau_sub  =  1.3049e+02 s",
        ]
        with open(self.path, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        self.assertEqual(pplib.step_times(self.path), {0: 0.0, 1: 1e-4})

    def test_restart_duplicates_last_wins_and_nan_skipped(self):
        fx.write_outp(self.run, [0, 1, 2, 1, 2, 3], [0.0, 1.0, 2.0, 1.5, 2.5, 3.5])
        with open(self.path, "a") as fh:
            fh.write("     4 |          nan | 1.000e-04 |  2.156e-08 |  2.531e-08 | -9.375e-07 "
                     "| 2.148e-11 |  1.839e-11 |  1.982e-05\n")
        self.assertEqual(pplib.step_times(self.path), {0: 0.0, 1: 1.5, 2: 2.5, 3: 3.5})

    def test_header_only(self):
        fx.write_outp(self.run, [], [])
        self.assertEqual(pplib.step_times(self.run), {})


# ---------------------------------------------------------------------------
# auto_time_unit, in_time_unit
# ---------------------------------------------------------------------------
class TestTimeUnits(unittest.TestCase):
    def test_boundaries(self):
        cases = [(0.0, "s"), (-5.0, "s"), (600.0, "s"), (np.nextafter(600.0, 1e9), "min"),
                 (7200.0, "min"), (np.nextafter(7200.0, 1e9), "h"), (3 * 86400.0, "h"),
                 (np.nextafter(3 * 86400.0, 1e9), "d"), (1e12, "d"), (math.inf, "d")]
        for t, u in cases:
            with self.subTest(t=t):
                self.assertEqual(pplib.auto_time_unit(t), u)

    def test_nan_falls_through_to_days(self):
        # Silent: a NaN t_max (all comparisons False) labels the axis in days.
        self.assertEqual(pplib.auto_time_unit(float("nan")), "d")

    def test_in_time_unit(self):
        t = np.array([0.0, 60.0, 3600.0, 86400.0])
        np.testing.assert_array_equal(pplib.in_time_unit(t, "s"), t)
        np.testing.assert_array_equal(pplib.in_time_unit(t, "min"), [0, 1, 60, 1440])
        np.testing.assert_array_equal(pplib.in_time_unit(t, "h"), [0, 1 / 60, 1, 24])
        np.testing.assert_array_equal(pplib.in_time_unit(t, "d"), [0, 1 / 1440, 1 / 24, 1])
        self.assertEqual(float(pplib.in_time_unit(7200, "h")), 2.0)
        with self.assertRaises(KeyError):
            pplib.in_time_unit(t, "days")


# ---------------------------------------------------------------------------
# crossings
# ---------------------------------------------------------------------------
class TestCrossings(unittest.TestCase):
    eps = 1.0
    x = np.linspace(-20.0, 20.0, 401) + 0.0371            # h = 0.1, nodes off 0

    def test_midcell_crossing_is_exact(self):
        # Symmetric about the crossing -> linear interpolation is exact.
        x0 = 0.0371 + 0.05
        c = pplib.crossings(self.x, tanh_profile(self.x, x0, self.eps))
        self.assertEqual(c.size, 1)
        self.assertAlmostEqual(c[0], x0, delta=1e-13)

    def test_general_crossing_error_is_small(self):
        # At the inflection point the profile is cubic, so the error is O(h^3/eps^2).
        h = 0.1
        for x0 in (0.0, 0.0123, 0.03, 0.0789):
            with self.subTest(x0=x0):
                c = pplib.crossings(self.x, tanh_profile(self.x, x0, self.eps))
                self.assertEqual(c.size, 1)
                self.assertLess(abs(c[0] - x0), h ** 3 / self.eps ** 2 / 10)

    def test_decreasing_profile_and_other_levels(self):
        x0 = 1.234
        phi = tanh_profile(-self.x, -x0, self.eps)         # ice on the left
        self.assertAlmostEqual(pplib.crossings(self.x, phi)[0], x0, delta=1e-4)
        lvl = 0.9                                        # 0.9 crossing at x0 + 2 eps atanh(0.8)
        c = pplib.crossings(self.x, tanh_profile(self.x, x0, self.eps), level=lvl)
        self.assertAlmostEqual(c[0], x0 + 2 * self.eps * math.atanh(0.8), delta=1e-3)

    def test_several_crossings(self):
        # eps 0.2 so the tails of one lens are < 1e-15 at the other's edges
        phi = fx.union(fx.lens(self.x, -8.0, -3.0, 0.2), fx.lens(self.x, 4.0, 11.0, 0.2))
        c = pplib.crossings(self.x, phi)
        np.testing.assert_allclose(c, [-8.0, -3.0, 4.0, 11.0], atol=2e-4)

    def test_no_crossing_flat_and_touching(self):
        self.assertEqual(pplib.crossings(self.x, np.full(self.x.size, 0.5)).size, 0)
        self.assertEqual(pplib.crossings(self.x, np.full(self.x.size, 0.2)).size, 0)
        self.assertEqual(pplib.crossings(np.array([0.0]), np.array([0.7])).size, 0)

    def test_nan_interval_is_skipped(self):
        phi = tanh_profile(self.x, 0.0, self.eps)
        k = np.searchsorted(self.x, 0.0)
        phi[k] = np.nan
        self.assertEqual(pplib.crossings(self.x, phi).size, 0)       # silent: no warning

    # BUG: postprocess/pplib.py:428 tests (y[i]-level)*(y[i+1]-level) < 0, so a
    # profile that equals the level EXACTLY on a grid node (any IC centered on a
    # node: tanh(0) = 0 gives phi = 0.5 exactly) has no crossing at all; the
    # interface silently disappears. Right: return that node once.
    @unittest.expectedFailure
    def test_crossing_exactly_on_a_node(self):
        x = np.linspace(0.0, 1.0, 11)
        phi = tanh_profile(x, 0.5, 0.1)
        self.assertEqual(phi[5], 0.5)
        np.testing.assert_allclose(pplib.crossings(x, phi), [0.5])


# ---------------------------------------------------------------------------
# refine_tanh
# ---------------------------------------------------------------------------
class TestRefineTanh(unittest.TestCase):
    eps = 1.0
    x = np.linspace(-20.0, 20.0, 401) + 0.0371

    def refine(self, x, phi, eps):
        c = pplib.crossings(x, phi)
        return c, [pplib.refine_tanh(x, phi, ci, eps) for ci in c]

    def test_exact_on_the_model_profile(self):
        for x0 in (0.0, 0.0123, 0.05, 1.7):
            phi = tanh_profile(self.x, x0, self.eps)
            c, r = self.refine(self.x, phi, self.eps)
            with self.subTest(x0=x0):
                self.assertAlmostEqual(r[0], x0, delta=1e-12)
                self.assertLessEqual(abs(r[0] - x0), abs(c[0] - x0))

    def test_decreasing_profile(self):
        phi = tanh_profile(-self.x, -0.3, self.eps)
        self.assertAlmostEqual(self.refine(self.x, phi, self.eps)[1][0], 0.3, delta=1e-12)

    def test_eps_argument_only_sets_the_window(self):
        phi = tanh_profile(self.x, 0.3, self.eps)
        for e in (0.5, 2.0):
            self.assertAlmostEqual(self.refine(self.x, phi, e)[1][0], 0.3, delta=1e-12)

    def test_band_thinner_than_five_samples_returns_guess(self):
        e = 0.02                                          # h = 0.1 = 5 eps
        phi = tanh_profile(self.x, 0.3, e)
        c, r = self.refine(self.x, phi, e)
        self.assertEqual(r[0], c[0])                     # documented fallback, bit-identical

    def test_near_domain_edge(self):
        xe = np.linspace(0.0, 10.0, 101)
        phi = tanh_profile(xe, 0.55, self.eps)            # band cut by x = 0
        self.assertAlmostEqual(self.refine(xe, phi, self.eps)[1][0], 0.55, delta=1e-12)

    def test_overshoot_outside_band_is_excluded(self):
        phi = tanh_profile(self.x, 0.3, self.eps)
        phi = np.where(phi > 0.995, 1.02, np.where(phi < 0.005, -0.02, phi))
        self.assertAlmostEqual(self.refine(self.x, phi, self.eps)[1][0], 0.3, delta=1e-12)

    def test_overshoot_just_outside_unit_interval_near_band(self):
        # Values in (0.98, 1] and [0, 0.02) are cut; values > 1 are cut too,
        # so ringing of +-0.03 just past the band edges does not bias the fit.
        phi = tanh_profile(self.x, 0.3, self.eps)
        ring = 0.03 * np.sin(5 * self.x) * (np.abs(self.x - 0.3) > 9.0)
        self.assertAlmostEqual(self.refine(self.x, phi + ring, self.eps)[1][0], 0.3, delta=1e-12)

    def test_lens_wider_than_window_is_exact(self):
        phi = fx.lens(self.x, -6.0, 6.0, self.eps)
        _, r = self.refine(self.x, phi, self.eps)
        np.testing.assert_allclose(r, [-6.0, 6.0], atol=1e-12)

    def test_thin_lens_bias_is_measured(self):
        # Two interfaces 8 eps apart both fall in the 6-eps window: the fit is
        # biased by a quarter eps. Silent (no flag). Recorded, not a bug.
        phi = fx.profile(np.minimum(self.x + 4.0, 4.0 - self.x), self.eps)
        _, r = self.refine(self.x, phi, self.eps)
        self.assertGreater(abs(r[0] + 4.0), 0.2)
        self.assertLess(abs(r[0] + 4.0), 0.3)

    # BUG: postprocess/pplib.py:451-454 returns x0 - intercept/slope with no
    # check that the root lies in the fitted window. On a lens 2-3 eps wide
    # (phi never reaches 0.98, both flanks in the window) the line is nearly
    # flat and the "interface" lands 8 to 370 eps away, outside the domain.
    # Right: fall back to x0 when |root - x0| > half_width*eps (or when the
    # band is not monotone).
    @unittest.expectedFailure
    def test_root_never_leaves_the_window(self):
        for w in (2.0, 3.0, 4.0):
            phi = fx.profile(np.minimum(self.x + w / 2, w / 2 - self.x), self.eps)
            c, r = self.refine(self.x, phi, self.eps)
            for ci, ri in zip(c, r):
                with self.subTest(w=w, x0=ci):
                    self.assertLessEqual(abs(ri - ci), 6.0 * self.eps)


# ---------------------------------------------------------------------------
# contour_points
# ---------------------------------------------------------------------------
class TestContourPoints(unittest.TestCase):
    def test_disc_points_lie_on_the_circle(self):
        x = np.linspace(0.0, 10.0, 201) + 0.013
        y = np.linspace(0.0, 8.0, 161) + 0.007
        X, Y = np.meshgrid(x, y)
        R, eps = 3.0, 0.1
        phi = fx.disc(X, Y, 5.0, 4.0, R, eps)
        xl, yl, xr, yr = pplib.contour_points(x, Y, phi)
        self.assertGreater(xl.size, 100)
        self.assertTrue(np.all(xl < xr))
        np.testing.assert_array_equal(yl, yr)
        r = np.hypot(np.r_[xl, xr] - 5.0, np.r_[yl, yr] - 4.0)
        np.testing.assert_allclose(r, R, atol=2e-3)

    def test_rows_with_one_or_three_crossings_are_dropped(self):
        x = np.linspace(0.0, 10.0, 101) + 0.0031
        lvl_rows = {
            0: fx.lens(x, 2.0, 5.0, 0.2),                                   # 2 crossings
            1: fx.profile(x - 4.0, 0.2),                                    # 1 crossing
            2: fx.union(fx.lens(x, 1.0, 3.0, 0.2), fx.profile(x - 7.0, 0.2)),  # 3
            3: fx.union(fx.lens(x, 1.0, 3.0, 0.2), fx.lens(x, 6.0, 8.0, 0.2)),  # 4
            4: np.zeros_like(x),                                            # 0
            5: fx.lens(x, 6.5, 8.5, 0.2),                                   # 2
        }
        phi = np.vstack([lvl_rows[k] for k in range(6)])
        Y = np.repeat(np.arange(6.0)[:, None], x.size, axis=1)
        xl, yl, xr, yr = pplib.contour_points(x, Y, phi)
        np.testing.assert_array_equal(yl, [0.0, 5.0])
        np.testing.assert_allclose(xl, [2.0, 6.5], atol=1e-3)
        np.testing.assert_allclose(xr, [5.0, 8.5], atol=1e-3)

    def test_no_crossings(self):
        x = np.linspace(0, 1, 5)
        out = pplib.contour_points(x, np.zeros((3, 5)), np.ones((3, 5)))
        self.assertTrue(all(a.size == 0 for a in out))
        out = pplib.contour_points(x, np.zeros((1, 5)), fx.profile(x - 0.53, 0.05)[None, :])
        self.assertTrue(all(a.size == 0 for a in out))       # one crossing only

    def test_node_hit_turns_a_two_crossing_row_into_one(self):
        # Consequence of the crossings() node bug: a row whose left meniscus
        # sits exactly on a node is dropped, silently.
        x = np.linspace(0.0, 10.0, 101)
        phi = fx.lens(x, 2.0, 5.03, 0.2)[None, :]
        self.assertEqual(phi[0, 20], 0.5)
        out = pplib.contour_points(x, np.zeros_like(phi), phi)
        self.assertEqual(out[0].size, 0)


# ---------------------------------------------------------------------------
# circle_radius
# ---------------------------------------------------------------------------
class TestCircleRadius(unittest.TestCase):
    def test_exact_arc(self):
        th = np.linspace(0.2, 2.5, 40)
        xc, yc, r = pplib.circle_radius(3 + 2 * np.cos(th), -1 + 2 * np.sin(th))
        self.assertAlmostEqual(xc, 3.0, delta=1e-12)
        self.assertAlmostEqual(yc, -1.0, delta=1e-12)
        self.assertAlmostEqual(r, 2.0, delta=1e-12)

    def test_micron_scale_and_large_radius(self):
        th = np.linspace(-0.5, 0.5, 30)
        R = 9.375e-5
        xc, yc, r = pplib.circle_radius(1e-4 + R * np.cos(th), 6e-5 + R * np.sin(th))
        self.assertAlmostEqual(r / R, 1.0, delta=1e-9)
        th = np.linspace(-1e-3, 1e-3, 20)
        _, _, r = pplib.circle_radius(1e6 * np.sin(th), 1e6 * np.cos(th) - 1e6)
        self.assertAlmostEqual(r / 1e6, 1.0, delta=1e-6)

    def test_noisy_arc(self):
        rng = np.random.default_rng(0)
        th = np.linspace(0.2, 2.5, 200)
        sig = 0.01
        px = 3 + 2 * np.cos(th) + rng.normal(0, sig, th.size)
        py = -1 + 2 * np.sin(th) + rng.normal(0, sig, th.size)
        xc, yc, r = pplib.circle_radius(px, py)
        self.assertAlmostEqual(r, 2.0, delta=5 * sig / math.sqrt(th.size) * 3)
        self.assertAlmostEqual(xc, 3.0, delta=0.02)
        self.assertAlmostEqual(yc, -1.0, delta=0.02)

    def test_point_count_limit(self):
        th = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        self.assertAlmostEqual(pplib.circle_radius(np.cos(th), np.sin(th))[2], 1.0, delta=1e-12)
        self.assertIsNone(pplib.circle_radius(np.cos(th[:4]), np.sin(th[:4])))
        self.assertIsNone(pplib.circle_radius(np.empty(0), np.empty(0)))

    # BUG: postprocess/pplib.py:492-502 never checks the rank or conditioning
    # of A. For collinear points A = [2x, 2y, 1] has rank 2, lstsq returns its
    # minimum-norm solution, and a finite "circle" comes back (radius 0.32 for
    # ten points on y = 0, 2.2 for ten on y = 2x + 1). Right: None (no finite
    # circle). Fix: return None when lstsq's rank < 3 (or when the singular
    # values show A is ill-conditioned).
    @unittest.expectedFailure
    def test_collinear_points_return_none(self):
        px = np.linspace(0.0, 1.0, 10)
        for py in (np.zeros(10), 2 * px + 1):
            with self.subTest(slope=py[1] - py[0]):
                self.assertIsNone(pplib.circle_radius(px, py))


# ---------------------------------------------------------------------------
# The real local run (read only): two independent outputs must agree.
# ---------------------------------------------------------------------------
@unittest.skipUnless(fx.real_run_available(), "real run absent: " + fx.REAL_RUN)
class TestRealRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ssa = pplib.load_ssa(fx.REAL_RUN)

    def test_ssa_steps_contiguous(self):
        steps = self.ssa[:, pplib.STEP]
        self.assertEqual(steps[0], 0)
        np.testing.assert_array_equal(np.diff(steps), 1)
        self.assertTrue(np.all(np.diff(self.ssa[:, pplib.TIME]) > 0))
        # dt column is the step to the next row's t (to the 7 digits of %e)
        t = self.ssa[:, pplib.TIME]
        np.testing.assert_allclose(np.diff(t), self.ssa[:-1, pplib.DT], rtol=1e-5, atol=2e-6 * t.max())

    def test_step_times_agree_with_ssa(self):
        tm = pplib.step_times(fx.REAL_RUN)
        steps = self.ssa[:, pplib.STEP].astype(int)
        self.assertEqual(sorted(tm), list(steps))
        np.testing.assert_allclose([tm[s] for s in steps], self.ssa[:, pplib.TIME], rtol=1e-5)

    def test_snapshot_time_and_saturation(self):
        fn = os.path.join(fx.REAL_RUN, "vtkOut", "solV_00070.vts")
        if not os.path.isfile(fn):
            self.skipTest("solV_00070.vts absent")
        self.assertEqual(pplib.step_of(fn), 70)
        pvd = read_text(os.path.join(fx.REAL_RUN, "pf.pvd"))
        t70 = float(re.search(r'timestep="([^"]+)"[^>]*solV_00070', pvd).group(1))
        self.assertAlmostEqual(t70 / self.ssa[70, pplib.TIME], 1.0, delta=1e-6)
        f, X, Y = pplib.read_vts(fn, want=("IcePhase", "Temperature", "VaporDensity"))
        self.assertEqual(f["IcePhase"].shape, X.shape)
        self.assertAlmostEqual(float(X.max()), 3.75e-4, delta=1e-9)
        self.assertAlmostEqual(float(Y.max() - Y.min()), 1.25e-4, delta=1e-9)
        # Isothermal -20 C, h = 0.99 reservoir, pore started saturated: the
        # supersaturation in the air must be O(1e-2), not ~1e10 (old rho_vs bug).
        air = f["IcePhase"] < 0.01
        s = pplib.supersaturation(f["VaporDensity"][air], f["Temperature"][air])
        self.assertLess(abs(float(np.median(s))), 0.02)


if __name__ == "__main__":
    unittest.main()
