"""Tests for the figure scripts and run_postprocess.sh on synthetic run folders.

Scripts: plot_mass.py, plot_ssa.py, plot_porosity.py, plot_phi_bounds.py,
plot_timestep.py, plot_fields.py (sol_*.dat -> .vts), vapor_movie.py (needs
ffmpeg), and run_postprocess.sh. Each script runs in a subprocess with the
project venv, from the project root, on a run folder in a TemporaryDirectory
built by postprocess/tests/fixtures.py; the numbers it prints or writes are
checked against values known exactly from the synthetic input. A few internal
functions are also called in-process to check the numbers behind a figure.

The real local run (fixtures.REAL_RUN) is only READ, never written.

Run from the project root:
    venv_cryosuction/bin/python -m unittest discover -s postprocess/tests -t . -p 'test_plot_scripts.py'

Tests marked expectedFailure document a bug; the comment above each starts
with "# BUG:" and names the line.
"""
from __future__ import annotations

import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import matplotlib
matplotlib.use("Agg")

from postprocess.tests import fixtures as fx  # noqa: E402

sys.path.insert(0, fx.POSTPROCESS)
import pplib  # noqa: E402

PETIGA_DIR = os.path.dirname(os.path.dirname(fx.PROJECT))


def run_script(name, *args, env=None):
    """Run postprocess/<name> with the venv python from the project root."""
    return subprocess.run([fx.PY, os.path.join(fx.POSTPROCESS, name), *map(str, args)],
                          cwd=fx.PROJECT, capture_output=True, text=True, env=env, timeout=120)


def append_bounds(run, rows):
    """Append solver-format BOUNDS lines (monitoring.c:146) to outp.txt."""
    with open(os.path.join(run, "outp.txt"), "a") as fh:
        for imin, imax in rows:
            fh.write(f"  BOUNDS: phi_ice [{imin:.4e}, {imax:.4e}]  "
                     f"phi_air [{1 - imax:.4e}, {1 - imin:.4e}]\n")


def disc_ssa_column(R, eps, L, n):
    """sub_interf/eps for a diffuse disc, by midpoint quadrature on an n x n grid."""
    x = (np.arange(n) + 0.5) * L / n - L / 2
    X, Y = np.meshgrid(x, x)
    phi = fx.disc(X, Y, 0.0, 0.0, R, eps)
    return float(np.sum(phi ** 2 * (1 - phi) ** 2) * (L / n) ** 2 / eps)


def uniform_fn(phi, rv):
    """Uniform phi, T = -20 C, uniform rho_v."""
    return lambda X, Y: (phi + 0 * X, -20.0 + 0 * X, rv + 0 * X)


def pvd_entries(run):
    """[(t, file)] in pf.pvd order, parsed here (not through pplib)."""
    with open(os.path.join(run, "pf.pvd")) as fh:
        txt = fh.read()
    return [(float(t), f) for t, f in re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', txt)]


class TmpRun(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def assertPng(self, path):
        self.assertTrue(os.path.isfile(path), path)
        with open(path, "rb") as fh:
            self.assertEqual(fh.read(8), b"\x89PNG\r\n\x1a\n")


# ---------------------------------------------------------------------------
# plot_mass.py
# ---------------------------------------------------------------------------
class TestPlotMass(TmpRun):
    def setUp(self):
        super().setUp()
        import plot_mass
        self.pm = plot_mass

    def test_preferred_path_numbers(self):
        run = fx.make_run(self.root)
        t = np.array([0.0, 10.0, 100.0, 1000.0])
        ice = np.array([2e-8, 1.9e-8, 1.8e-8, 1.7e-8])
        rhov = np.array([2e-11, 3e-11, 4e-11, 5e-11])
        fx.write_ssa(run, t, ice, rhov=rhov, air=3e-8 - ice)
        times, mi, mv, dim = self.pm.compute_masses(run)
        np.testing.assert_allclose(times, t, rtol=1e-6)
        np.testing.assert_allclose(mi, 919.0 * ice, rtol=1e-6)
        np.testing.assert_allclose(mv, rhov, rtol=1e-6)
        self.assertEqual(dim, 2)

    def test_script_writes_all_figures(self):
        run = fx.make_run(self.root)
        fx.write_ssa(run, np.linspace(0, 5000, 30), np.linspace(2e-8, 1.8e-8, 30),
                     rhov=np.linspace(2e-11, 3e-11, 30))
        r = run_script("plot_mass.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertPng(os.path.join(run, "mass.png"))
        for f in ("total", "ice", "vapor", "change_loglog"):
            self.assertPng(os.path.join(run, "mass_plots", f + ".png"))

    def test_truncated_ssa_falls_back_to_outp(self):
        run = fx.make_run(self.root)
        fx.write_ssa(run, np.arange(5.0), np.full(5, 1e-8))
        fx.write_outp(run, np.arange(40), np.arange(40) * 2.0)
        times, mi, mv, _ = self.pm.compute_masses(run)
        self.assertEqual(times.size, 40)
        np.testing.assert_allclose(times, np.arange(40) * 2.0, rtol=1e-5)
        np.testing.assert_allclose(mi, 919.0 * 2.156e-08)          # TOT_ICE column
        np.testing.assert_allclose(mv, 2.148e-11)                   # TOT_RHOV column
        r = run_script("plot_mass.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("SSA_evo.dat truncated", r.stderr)

    def _fallback_run(self, nx, ny, phi_vals, rv=1e-3, ssa_steps=10):
        run = fx.make_run(self.root, geom=f"fb{nx}")
        nrb = fx.iga_rectangle(2e-4, 1e-4, nx, ny)
        snaps = [(s, uniform_fn(p, rv)) for s, p in phi_vals]
        fx.write_iga_run(run, nrb, snaps)
        if ssa_steps:
            fx.write_ssa(run, np.arange(ssa_steps) * 3.0, np.ones(ssa_steps), ncol=5)
        return run

    def test_fallback_time_matching_by_step(self):
        # Snapshots at steps 0, 5 and 12; SSA has steps 0..9 at t = 3*step.
        run = self._fallback_run(4, 2, [(0, 1.0), (5, 0.5), (12, 0.25)])
        times, mi, mv, _ = self.pm.compute_masses(run)
        np.testing.assert_array_equal(times, [0.0, 15.0])            # step 12 dropped silently
        self.assertAlmostEqual(mi[1] / mi[0], 0.5, places=12)
        self.assertEqual(mv[0], 0.0)                                  # phi = 1 -> no air

    def test_fallback_without_ssa_uses_snapshot_index_as_time(self):
        # Silent: with no SSA_evo.dat and no outp.txt, the x axis is the
        # snapshot INDEX labelled in seconds (plot_mass.py:245).
        run = self._fallback_run(4, 2, [(0, 1.0), (50, 0.5)], ssa_steps=0)
        times, _, _, _ = self.pm.compute_masses(run)
        np.testing.assert_array_equal(times, [0.0, 1.0])

    # BUG: postprocess/plot_mass.py:209-217 weights each CONTROL POINT by
    # |det| of np.gradient of the control net. On an open-knot B-spline the
    # end spacing is half the interior one, so the weights sum to
    # (Lx + hx/2)(Ly + hy/2): a uniform phi = 1 on a 20 x 10 element rectangle
    # reports 1.07625 x the true ice mass (7.6 % high; (1 + 1/(2n)) per
    # direction). Confirmed by hand with compute_masses. Right: the exact
    # integral of the spline. Fix: weight coefficient i by
    # (U[i+p+1] - U[i])/(p+1) per direction (the integral of N_i), times the
    # metric, or evaluate the spline at Gauss points. Only the fallback path
    # (runs whose SSA_evo.dat lacks the integral columns) is affected.
    @unittest.expectedFailure
    def test_fallback_uniform_field_integrates_exactly(self):
        run = self._fallback_run(20, 10, [(0, 1.0), (3, 0.5)])
        _, mi, mv, _ = self.pm.compute_masses(run)
        area = 2e-4 * 1e-4
        self.assertAlmostEqual(mi[0] / (919.0 * area), 1.0, delta=1e-6)
        self.assertAlmostEqual(mv[1] / (1e-3 * 0.5 * area), 1.0, delta=1e-6)

    def test_fallback_bias_is_first_order_in_h(self):
        # Records the size of the error above: (1 + 1/(2 nx))(1 + 1/(2 ny)).
        for nx, ny in ((20, 10), (40, 20)):
            run = self._fallback_run(nx, ny, [(0, 1.0)])
            _, mi, _, _ = self.pm.compute_masses(run)
            with self.subTest(nx=nx):
                self.assertAlmostEqual(mi[0] / (919.0 * 2e-8),
                                       (1 + 0.5 / nx) * (1 + 0.5 / ny), places=10)

    def test_no_inputs_fails_loudly(self):
        run = fx.make_run(self.root)
        r = run_script("plot_mass.py", "--dir", run)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("IGA geometry file not found", r.stderr)

    # BUG: postprocess/plot_mass.py:439 prints "Log-log change plot saved to"
    # unconditionally, but _change_loglog_plot (l. 323-326, 341-343) returns
    # without saving when no t > 0 exists (a one-row run) or every m(0) == 0.
    # The log says a file exists that does not (or a stale one from an earlier
    # run is left in place). Right: report only when the file was written.
    @unittest.expectedFailure
    def test_one_row_run_does_not_claim_a_missing_file(self):
        run = fx.make_run(self.root)
        fx.write_ssa(run, [0.0], [1e-8], rhov=[1e-11])
        r = run_script("plot_mass.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertPng(os.path.join(run, "mass.png"))
        claimed = "change_loglog.png" in r.stdout
        exists = os.path.isfile(os.path.join(run, "mass_plots", "change_loglog.png"))
        self.assertEqual(claimed, exists)

    # BUG: postprocess/plot_mass.py:376 picks the mass unit from the IGA
    # dimension only. Under -axisym the solver's integrals carry 2*pi*r
    # (src/assembly.c:399-409), so tot_ice is a VOLUME and the mass is in kg,
    # but the y label says kg m^-1. Right: "kg" when -axisym is set (and the
    # sol_*.dat fallback would also need the 2*pi*r weight).
    @unittest.expectedFailure
    def test_axisymmetric_mass_unit(self):
        import matplotlib.pyplot as plt
        run = fx.make_run(self.root, geometry={"-axisym": "1"})
        fx.write_iga_run(run, fx.iga_rectangle(2e-4, 1e-4, 4, 2), [])
        fx.write_ssa(run, [0.0, 1.0, 2.0], [1e-12, 1e-12, 1e-12], rhov=[1e-15] * 3)
        orig = self.pm.plt.close
        self.pm.plt.close = lambda *a, **k: None
        try:
            plt.close("all")
            self.pm.plot_mass(run)
            labels = [ax.get_ylabel() for n in plt.get_fignums() for ax in plt.figure(n).axes]
        finally:
            self.pm.plt.close = orig
            plt.close("all")
        self.assertTrue(labels)
        for lab in labels:
            self.assertEqual(lab, "Mass  [kg]")


# ---------------------------------------------------------------------------
# plot_ssa.py
# ---------------------------------------------------------------------------
class TestPlotSSA(TmpRun):
    def test_interface_factor_recovers_disc_perimeter(self):
        # Independent check of INTERFACE_FACTOR = 6: for a resolved diffuse
        # disc, 6 * sub_interf/eps equals 2 pi R.
        import plot_ssa
        R, eps = 1.0, 0.02
        col0 = disc_ssa_column(R, eps, 2.6, 1300)
        self.assertAlmostEqual(plot_ssa.INTERFACE_FACTOR * col0 / (2 * math.pi * R), 1.0, delta=1e-4)

    def test_power_law_exponent_and_values(self):
        run = fx.make_run(self.root)
        t = np.r_[0.0, np.geomspace(1.0, 1e5, 60)]
        col0 = np.r_[1e-3, 1e-3 * np.maximum(t[1:], 1.0) ** -0.25]
        col0[0] = 1e-3
        fx.write_ssa(run, t, np.full(t.size, 1e-8), ssa=col0)
        out = os.path.join(self.root, "figs", "ssa.png")
        r = run_script("plot_ssa.py", "--dir", run, "--save", out)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertPng(out)
        n = float(re.search(r"exponent n = ([-+0-9.]+)", r.stdout).group(1))
        self.assertAlmostEqual(n, -0.25, places=4)
        first, last = map(float, re.search(r"surface: (\S+) → (\S+)", r.stdout).groups())
        self.assertAlmostEqual(first / 6e-3, 1.0, places=5)
        self.assertAlmostEqual(last / (6e-3 * 1e5 ** -0.25), 1.0, places=5)

    def test_fit_power_law_direct(self):
        import plot_ssa
        t = np.geomspace(10.0, 1e4, 40)
        fit = plot_ssa.fit_power_law(t, 3.0 * t ** 0.4, 0.0, 1e9)
        self.assertAlmostEqual(fit["n"], 0.4, places=12)
        self.assertAlmostEqual(fit["C"], 3.0, places=10)
        self.assertAlmostEqual(fit["slope_lo"], 0.4, places=10)
        self.assertAlmostEqual(fit["slope_hi"], 0.4, places=10)
        self.assertIsNone(plot_ssa.fit_power_law(t[:2], t[:2], 0.0, 1e9))

    def test_analytic_reference_from_grains_2d_and_3d(self):
        R = np.array([1e-5, 2e-5, 3e-5])
        for dim, a0 in ((2, 2 * math.pi * R.sum()), (3, 4 * math.pi * (R ** 2).sum())):
            run = fx.make_run(os.path.join(self.root, f"d{dim}"),
                              geometry={"-grains_file": "inputs/g/grains.dat", "-dim": str(dim)})
            with open(os.path.join(run, "grains.dat"), "w") as fh:
                fh.write("2e-4 1e-4\n" + "".join(f"0 0 {r}\n" for r in R))
            fx.write_ssa(run, [0.0, 1.0, 2.0], [1e-8] * 3, ssa=[a0 / 6.0] * 3)
            r = run_script("plot_ssa.py", "--dir", run)
            with self.subTest(dim=dim):
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertIn("(3 grain rows); measured/analytic = 1.000", r.stdout)
                self.assertIn("m$^2$" if dim == 3 else " m\n", r.stdout)

    def test_legacy_four_columns_and_one_row(self):
        run = fx.make_run(self.root)
        fx.write_ssa(run, [0.0, 5.0, 10.0, 20.0, 40.0], [1e-8] * 5, ssa=[1e-4, 9e-5, 8e-5, 7e-5, 6e-5],
                     ncol=4)
        r = run_script("plot_ssa.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertPng(os.path.join(run, "ssa.png"))
        run1 = fx.make_run(os.path.join(self.root, "one"))
        fx.write_ssa(run1, [0.0], [1e-8], ssa=[1e-4])
        r = run_script("plot_ssa.py", "--dir", run1)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("power-law fit", r.stdout)

    def test_missing_or_zero_fails_with_message(self):
        run = fx.make_run(self.root)
        r = run_script("plot_ssa.py", "--dir", run)
        self.assertEqual(r.returncode, 1)
        self.assertIn("No usable SSA_evo.dat", r.stderr)
        fx.write_ssa(run, [0.0, 1.0], [1e-8, 1e-8])               # ssa column all 0
        r = run_script("plot_ssa.py", "--dir", run)
        self.assertEqual(r.returncode, 1)
        self.assertIn("zero throughout", r.stderr)

    # BUG: postprocess/plot_ssa.py:188-191 normalizes the right axis by
    # area[0] when there is no packing file; when the first interface value is
    # 0 and later ones are positive (ice that appears after t = 0) the limits
    # become inf and matplotlib raises "Axis limits cannot be NaN or Inf"
    # (traceback, exit 1, no figure). The title (l. 200) divides by area[0]
    # too. Right: draw the figure, skipping the normalized axis.
    @unittest.expectedFailure
    def test_first_value_zero_still_plots(self):
        run = fx.make_run(self.root)
        fx.write_ssa(run, np.arange(6.0), [1e-8] * 6, ssa=np.arange(6.0) * 1e-5)
        r = run_script("plot_ssa.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr[-400:])


# ---------------------------------------------------------------------------
# plot_porosity.py
# ---------------------------------------------------------------------------
class TestPlotPorosity(TmpRun):
    def porosity(self, r):
        return tuple(map(float, re.search(r"porosity: (\S+) → (\S+)", r.stdout).groups()))

    def test_exact_from_ice_and_air(self):
        run = fx.make_run(self.root)
        ice = np.array([1.0e-8, 1.2e-8, 1.5e-8])
        air = 4e-8 - ice
        fx.write_ssa(run, [0.0, 1e3, 1e4], ice, air=air)
        r = run_script("plot_porosity.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertPng(os.path.join(run, "porosity.png"))
        p0, p1 = self.porosity(r)
        self.assertAlmostEqual(p0, 0.75, places=6)
        self.assertAlmostEqual(p1, 0.625, places=6)
        self.assertIn("tot_air / (tot_ice + tot_air)", r.stdout)

    def test_axisymmetric_new_format_uses_the_ratio(self):
        run = fx.make_run(self.root, geometry={"-axisym": "1"})
        fx.write_ssa(run, [0.0, 1.0], [1e-12, 2e-12], air=[3e-12, 2e-12])
        r = run_script("plot_porosity.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.porosity(r), (0.75, 0.5))

    def test_legacy_uses_domain_measure(self):
        run = fx.make_run(self.root, geometry={"-Lx": "2e-4", "-Ly": "1e-4"})
        fx.write_ssa(run, [0.0, 1.0], [5e-9, 1e-8], ncol=5)
        r = run_script("plot_porosity.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.porosity(r), (0.75, 0.5))
        self.assertIn("(legacy)", r.stdout)

    def test_legacy_3d_needs_lz(self):
        run = fx.make_run(self.root, geometry={"-Lx": "1e-4", "-Ly": "1e-4", "-Lz": "1e-4",
                                               "-dim": "3"})
        fx.write_ssa(run, [0.0, 1.0], [2.5e-13, 5e-13], ncol=4)
        r = run_script("plot_porosity.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.porosity(r), (0.75, 0.5))

    def test_unrecoverable_cases_fail_with_message(self):
        run = fx.make_run(self.root)                           # no -Lx/-Ly
        fx.write_ssa(run, [0.0, 1.0], [5e-9, 1e-8], ncol=5)
        r = run_script("plot_porosity.py", "--dir", run)
        self.assertEqual(r.returncode, 1)
        self.assertIn("Cannot compute porosity", r.stderr)
        run2 = fx.make_run(os.path.join(self.root, "ax"),
                           geometry={"-axisym": "1", "-Lx": "1e-4", "-Ly": "1e-4"})
        fx.write_ssa(run2, [0.0, 1.0], [5e-9, 1e-8], ncol=5)
        r = run_script("plot_porosity.py", "--dir", run2)
        self.assertEqual(r.returncode, 1)
        self.assertIn("axisymmetric", r.stderr)
        r = run_script("plot_porosity.py", "--dir", os.path.join(self.root, "nothing"))
        self.assertEqual(r.returncode, 1)


# ---------------------------------------------------------------------------
# plot_phi_bounds.py
# ---------------------------------------------------------------------------
class TestPlotPhiBounds(TmpRun):
    def read_csv(self, run):
        with open(os.path.join(run, "phi_bounds.csv")) as fh:
            return list(csv.DictReader(fh))

    def test_values_csv_and_time_axis(self):
        run = fx.make_run(self.root, solver={"-phase_lo": "-0.1"})
        t = np.array([0.0, 86400.0, 2 * 86400.0])
        fx.write_ssa(run, t, [1e-8] * 3)
        fx.write_outp(run, [0, 1, 2], t)
        append_bounds(run, [(0.0, 1.0), (-2.5e-23, 1.0), (-3.2e-5, 1.0 + 4.0e-7)])
        png = os.path.join(self.root, "pb.png")
        r = run_script("plot_phi_bounds.py", "--dir", run, "--save", png)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertPng(png)
        rows = self.read_csv(run)
        self.assertEqual(len(rows), 3)
        np.testing.assert_allclose([float(x["x"]) for x in rows], [0, 1, 2])     # days
        self.assertEqual(float(rows[2]["phi_ice_min"]), -3.2e-5)
        self.assertEqual(float(rows[2]["undershoot"]), 3.2e-5)
        self.assertIn("worst excursion    3.200e-05", r.stdout)
        self.assertIn("REAL, inspect the run", r.stdout)
        self.assertIn("(guard -0.100)", r.stdout)

    # BUG: postprocess/plot_phi_bounds.py:83 takes the overshoot from phi_ice
    # max only. The solver prints it as %.4e (monitoring.c:146), so any
    # phi_ice max below 1 + 5e-5 prints as "1.0000e+00" and the overshoot reads
    # 0; the same excursion is printed at full precision as phi_air min
    # (= 1 - phi_ice max, here -4.0000e-07), which the script ignores. The
    # real run shows it: "phi_air [-4.4409e-16, ...]" next to phi_ice max
    # "1.0000e+00". Right: overshoot = max(phi_ice_max - 1, -phi_air_min, 0).
    @unittest.expectedFailure
    def test_small_overshoot_is_resolved(self):
        run = fx.make_run(self.root)
        fx.write_outp(run, [0, 1], [0.0, 1.0])
        append_bounds(run, [(0.0, 1.0), (0.0, 1.0 + 4.0e-7)])
        with open(os.path.join(run, "outp.txt")) as fh:
            self.assertIn("phi_ice [0.0000e+00, 1.0000e+00]  phi_air [-4.0000e-07", fh.read())
        r = run_script("plot_phi_bounds.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertAlmostEqual(float(self.read_csv(run)[1]["overshoot"]), 4.0e-7, delta=1e-12)

    def test_rounding_only(self):
        run = fx.make_run(self.root)
        fx.write_outp(run, [0, 1], [0.0, 1.0])
        append_bounds(run, [(0.0, 1.0), (-2.6e-23, 1.0)])
        r = run_script("plot_phi_bounds.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("rounding only", r.stdout)

    def test_rejected_attempt_switches_to_step_index(self):
        # A step that leaves [phase_lo, phase_hi] prints BOUNDS, writes no SSA
        # row and is retried (monitoring.c:146-191), so BOUNDS has one line more
        # than SSA_evo.dat. The script then uses the line index as "step" and
        # reports the REJECTED attempt's excursion as the run's worst.
        run = fx.make_run(self.root)
        fx.write_ssa(run, [0.0, 1.0, 2.0], [1e-8] * 3)
        append_bounds(run, [(0.0, 1.0), (-0.2, 1.0), (-1e-9, 1.0), (0.0, 1.0)])
        r = run_script("plot_phi_bounds.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = self.read_csv(run)
        self.assertEqual([int(x["step"]) for x in rows], [0, 1, 2, 3])
        self.assertEqual([float(x["x"]) for x in rows], [0, 1, 2, 3])
        self.assertIn("worst excursion    2.000e-01", r.stdout)

    def test_missing_outp_or_no_bounds_fails(self):
        run = fx.make_run(self.root)
        r = run_script("plot_phi_bounds.py", "--dir", run)
        self.assertEqual(r.returncode, 1)
        self.assertIn("no outp.txt", r.stderr)
        fx.write_outp(run, [0], [0.0])
        r = run_script("plot_phi_bounds.py", "--dir", run)
        self.assertEqual(r.returncode, 1)
        self.assertIn("no BOUNDS lines", r.stderr)

    # BUG: postprocess/plot_phi_bounds.py:82-84 cannot see a NaN field. The
    # solver's min/max loop skips NaN (every comparison is false), so an
    # all-NaN phi prints "phi_ice [1.0000e+30, -1.0000e+30]"
    # (monitoring.c:107-108, 146); a printed "nan" would not match _RE (l. 37)
    # and the line is dropped. Either way under = over = 0 and the script
    # reports "rounding only". Right: flag min > max, |phi| >= 1e29, or a
    # BOUNDS line that fails to parse.
    @unittest.expectedFailure
    def test_nan_field_is_not_reported_healthy(self):
        run = fx.make_run(self.root)
        fx.write_outp(run, [0, 1], [0.0, 1.0])
        append_bounds(run, [(0.0, 1.0), (1.0e30, -1.0e30)])
        r = run_script("plot_phi_bounds.py", "--dir", run)
        self.assertNotIn("rounding only", r.stdout)


# ---------------------------------------------------------------------------
# plot_timestep.py
# ---------------------------------------------------------------------------
class TestPlotTimestep(TmpRun):
    def test_load_monitor_ignores_other_tables(self):
        import plot_timestep
        run = fx.make_run(self.root)
        fx.write_outp(run, np.arange(30), np.arange(30) * 0.5, every_header=8)
        with open(os.path.join(run, "outp.txt"), "a") as fh:
            fh.write("     0 |  5.3957e-04 |   1.108e-15 | 1.000e+00 | 0.000e+00 |   9.564e-13 "
                     "| 1.000e+00 | 0.000e+00 |   5.396e-04 | 1.000e+00 | 0.000e+00\n")
            fh.write("       |              |           |    -0.001% |    +0.001% |           "
                     "|   +0.000% |            |    -0.000%\n")
        d = plot_timestep.load_monitor(os.path.join(run, "outp.txt"))
        self.assertEqual(d.shape, (30, 3))
        np.testing.assert_array_equal(d[:, 0], np.arange(30))
        np.testing.assert_allclose(d[:, 1], np.arange(30) * 0.5)
        np.testing.assert_allclose(d[:, 2], 1e-4)

    def test_script_saves_figure(self):
        run = fx.make_run(self.root)
        fx.write_outp(run, np.arange(10), np.geomspace(1e-4, 1e3, 10))
        png = os.path.join(self.root, "ts.png")
        r = run_script("plot_timestep.py", "--dir", run, "--save", png)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertPng(png)

    def test_failures_are_loud(self):
        run = fx.make_run(self.root)
        r = run_script("plot_timestep.py", "--dir", run, "--save", os.path.join(self.root, "x.png"))
        self.assertEqual(r.returncode, 1)
        self.assertIn("monitor file not found", r.stderr)
        fx.write_outp(run, [0], [0.0])
        r = run_script("plot_timestep.py", "--dir", run, "--save", os.path.join(self.root, "x.png"))
        self.assertEqual(r.returncode, 1)
        self.assertIn("all time values are zero", r.stderr)

    @unittest.skipUnless(fx.real_run_available(), "real run absent")
    def test_real_run_matches_ssa(self):
        import plot_timestep
        d = plot_timestep.load_monitor(os.path.join(fx.REAL_RUN, "outp.txt"))
        ssa = pplib.load_ssa(fx.REAL_RUN)
        np.testing.assert_array_equal(d[:, 0], ssa[:, pplib.STEP])
        np.testing.assert_allclose(d[:, 2], ssa[:, pplib.DT], rtol=6e-4)       # %9.3e


# ---------------------------------------------------------------------------
# plot_fields.py
# ---------------------------------------------------------------------------
def field_fn(a):
    """A nonlinear test field set (control values are written as given)."""
    return lambda X, Y: (np.sin(X * 1e4 + a) * np.cos(Y * 3e4), -20.0 + 1e5 * X - 3e5 * Y,
                         1e-3 * (1 + a) + 1e-1 * X * Y / 1e-8)


class TestPlotFields(TmpRun):
    def setUp(self):
        super().setUp()
        self.run = fx.make_run(self.root)
        self.nrb = fx.iga_rectangle(3e-4, 1e-4, 6, 3)            # 8 x 5 control points

    def pvd(self):
        return pvd_entries(self.run)

    def test_roundtrip_values_coordinates_and_times(self):
        fx.write_iga_run(self.run, self.nrb, [(0, field_fn(0.0)), (7, field_fn(1.0))])
        fx.write_outp(self.run, np.arange(10), np.arange(10) * 12.5)
        r = run_script("plot_fields.py", "--dir", self.run)
        self.assertEqual(r.returncode, 0, r.stderr)
        GX, GY = self.nrb.points[..., 0], self.nrb.points[..., 1]          # (nx, ny)
        for step, a in ((0, 0.0), (7, 1.0)):
            f, X, Y = pplib.read_vts(os.path.join(self.run, "vtkOut", f"solV_{step:05d}.vts"))
            self.assertEqual(X.shape, (5, 8))
            np.testing.assert_array_equal(X, GX.T)
            np.testing.assert_array_equal(Y, GY.T)
            phi, T, rv = field_fn(a)(GX, GY)
            np.testing.assert_array_equal(f["IcePhase"], phi.T)
            np.testing.assert_array_equal(f["Temperature"], T.T)
            np.testing.assert_array_equal(f["VaporDensity"], rv.T)
            self.assertNotIn("AirPhase", f)                             # 3 DOFs only
        times = [t for t, _ in self.pvd()]
        self.assertEqual(times, [0.0, 87.5])

    def test_one_snapshot_without_outp_uses_step_as_time(self):
        fx.write_iga_run(self.run, self.nrb, [(40, field_fn(0.0))])
        r = run_script("plot_fields.py", "--dir", self.run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("step index as time proxy", r.stdout)
        self.assertEqual([t for t, _ in self.pvd()], [40.0])

    # BUG: postprocess/plot_fields.py:275 uses t = float(step) for a snapshot
    # whose step is missing from outp.txt (truncated log), with no warning
    # when the map is non-empty; _write_pvd (l. 230) then SORTS by that mixed
    # seconds/step value, so ParaView plays the snapshots out of order (here
    # step 12 at "t = 12 s" lands between steps 0 and 5 at 0 and 500 s).
    # Right: keep step order, and fall back to SSA_evo.dat or warn.
    @unittest.expectedFailure
    def test_partial_outp_keeps_step_order(self):
        fx.write_iga_run(self.run, self.nrb, [(0, field_fn(0)), (5, field_fn(0)), (12, field_fn(0))])
        fx.write_outp(self.run, np.arange(10), np.arange(10) * 100.0)
        r = run_script("plot_fields.py", "--dir", self.run)
        self.assertEqual(r.returncode, 0, r.stderr)
        steps = [pplib.step_of(f) for _, f in self.pvd()]
        self.assertEqual(steps, [0, 5, 12])

    def test_existing_vts_is_not_refreshed_without_force(self):
        fx.write_iga_run(self.run, self.nrb, [(3, field_fn(0.0))])
        run_script("plot_fields.py", "--dir", self.run)
        fx.write_iga_run(self.run, self.nrb, [(3, field_fn(2.0))])        # sol rewritten
        run_script("plot_fields.py", "--dir", self.run)
        fn = os.path.join(self.run, "vtkOut", "solV_00003.vts")
        GX, GY = self.nrb.points[..., 0], self.nrb.points[..., 1]
        stale = field_fn(0.0)(GX, GY)[2].T
        np.testing.assert_array_equal(pplib.read_vts(fn)[0]["VaporDensity"], stale)
        run_script("plot_fields.py", "--dir", self.run, "--force")
        np.testing.assert_array_equal(pplib.read_vts(fn)[0]["VaporDensity"],
                                      field_fn(2.0)(GX, GY)[2].T)

    def test_no_snapshots_and_no_geometry(self):
        fx.write_iga_run(self.run, self.nrb, [])
        fx.write_solver_evo(self.run, [0, 1], [0.0, 1.0])                  # not a snapshot
        r = run_script("plot_fields.py", "--dir", self.run)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("No sol_*.dat files", r.stdout)
        self.assertFalse(os.path.exists(os.path.join(self.run, "pf.pvd")))
        r = run_script("plot_fields.py", "--dir", os.path.join(self.root, "empty"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("IGA geometry file not found", r.stderr)

    @unittest.skipUnless(fx.real_run_available(), "real run absent")
    def test_real_run_vts_equals_conversion_of_sol(self):
        # Read-only: convert sol_00070.dat in memory and compare with the
        # shipped vtkOut/solV_00070.vts.
        import plot_fields
        from igakit.io import PetIGA
        sol_fn = os.path.join(fx.REAL_RUN, "sol_00070.dat")
        vts_fn = os.path.join(fx.REAL_RUN, "vtkOut", "solV_00070.vts")
        if not (os.path.isfile(sol_fn) and os.path.isfile(vts_fn)):
            self.skipTest("snapshot 70 absent")
        nrb = PetIGA().read(os.path.join(fx.REAL_RUN, "igasol.dat"))
        sol = PetIGA().read_vec(sol_fn, nrb)
        f, X, Y = pplib.read_vts(vts_fn)
        ny, nx = X.shape
        np.testing.assert_array_equal(f["IcePhase"].ravel(), plot_fields._vtk_scalar(sol[..., 0]))
        np.testing.assert_array_equal(f["VaporDensity"].ravel(), plot_fields._vtk_scalar(sol[..., 2]))
        np.testing.assert_array_equal(X, nrb.points[..., 0].T)
        self.assertEqual(sol.shape[:2], (nx, ny))


# ---------------------------------------------------------------------------
# vapor_movie.py
# ---------------------------------------------------------------------------
@unittest.skipUnless(shutil.which("ffmpeg"), "ffmpeg not on PATH")
class TestVaporMovie(TmpRun):
    def make(self, radii, ny=31):
        run = fx.make_run(self.root)
        X, Y = fx.grid(60e-6, 30e-6, 61, ny)
        rvs = float(pplib.rho_vs(-20.0))
        snaps = []
        for k, R in enumerate(radii):
            phi = fx.disc(X, Y, 30e-6, 15e-6, R, 1e-6)
            rv = rvs * (1.0 + 0.01 * (X / 60e-6 - 0.5) * (k + 1))
            snaps.append((k * 10, k * 3600.0, phi, np.full(X.shape, -20.0), rv))
        fx.write_vts_run(run, X, Y, snaps)
        return run, X, Y, snaps, rvs

    def test_movie_color_limit_and_ice_area(self):
        run, X, Y, snaps, rvs = self.make([10e-6, 9e-6, 8e-6])
        r = run_script("vapor_movie.py", "--dir", run)
        self.assertEqual(r.returncode, 0, r.stderr)
        for f in ("vapor_movie.mp4", "vapor_movie_first.png", "vapor_movie_mid.png"):
            self.assertGreater(os.path.getsize(os.path.join(run, "plots", f)), 1000)
        # Color limit: 99.5th percentile of |sigma| over the air of all snapshots
        # (every sample is used when there are fewer than 50000).
        s = np.concatenate([np.abs(rv / rvs - 1.0)[phi < 0.5] for _, _, phi, _, rv in snaps])
        lim = float(re.search(r"sigma color limit \+/-(\S+)", r.stdout).group(1))
        self.assertAlmostEqual(lim / np.percentile(s, 99.5), 1.0, delta=1e-3)
        # Ice area at the end vs the exact diffuse-disc integral pi (R^2 + pi^2 eps^2 / 3).
        pct = float(re.search(r"ice area ([0-9.]+) % of initial at the end", r.stdout).group(1))
        exact = (8.0 ** 2 + math.pi ** 2 / 3) / (10.0 ** 2 + math.pi ** 2 / 3) * 100
        self.assertAlmostEqual(pct, exact, delta=0.2)

    def test_single_snapshot_with_tracers(self):
        run, *_ = self.make([10e-6])
        out = os.path.join(self.root, "m.mp4")
        r = run_script("vapor_movie.py", "--dir", run, "--tracers", "--n-particles", "40",
                       "--subframes", "1", "--grid", "2", "--save", out)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertGreater(os.path.getsize(out), 1000)
        self.assertTrue(os.path.isfile(os.path.join(self.root, "m_first.png")))
        self.assertIn("1 frames", r.stdout)

    def test_refuses_1d_and_missing_pvd(self):
        run, *_ = self.make([10e-6], ny=1)
        r = run_script("vapor_movie.py", "--dir", run)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("needs a 2D run", r.stderr)
        os.remove(os.path.join(run, "pf.pvd"))
        r = run_script("vapor_movie.py", "--dir", run)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("run plot_fields.py first", r.stderr)


# ---------------------------------------------------------------------------
# run_postprocess.sh
# ---------------------------------------------------------------------------
class TestRunPostprocess(TmpRun):
    SCRIPTS = ("pplib.py", "plot_fields.py", "plot_porosity.py", "plot_ssa.py",
               "plot_timestep.py", "plot_mass.py", "plot_phi_bounds.py", "run_postprocess.sh")

    def stage(self, outp=True, ssa=True, iga=True):
        run = fx.make_run(self.root)
        pp = os.path.join(run, "postprocess")
        os.makedirs(pp)
        for f in self.SCRIPTS:
            shutil.copy(os.path.join(fx.POSTPROCESS, f), pp)
        steps = np.arange(6)
        t = np.r_[0.0, np.geomspace(1.0, 1e4, 5)]
        if ssa:
            fx.write_ssa(run, t, np.linspace(2e-8, 1.9e-8, 6), rhov=2e-11, air=np.full(6, 1e-8),
                         ssa=np.linspace(5e-5, 4e-5, 6))
        if outp:
            fx.write_outp(run, steps, t)
            append_bounds(run, [(0.0, 1.0)] * 6)
        if iga:
            fx.write_iga_run(run, fx.iga_rectangle(2e-4, 1e-4, 4, 2), [(0, field_fn(0)), (5, field_fn(1))])
        return run

    def bash(self, run):
        env = dict(os.environ, PETIGA_DIR=PETIGA_DIR)
        return subprocess.run(["bash", os.path.join(run, "postprocess", "run_postprocess.sh"), run],
                              capture_output=True, text=True, env=env, timeout=180)

    def test_full_run_exits_zero_and_writes_figures(self):
        run = self.stage()
        r = self.bash(run)
        self.assertEqual(r.returncode, 0, r.stdout[-2000:])
        for f in ("porosity", "ssa", "timestep", "mass", "phi_bounds"):
            self.assertPng(os.path.join(run, "plots", f + ".png"))
        for f in ("total", "ice", "vapor", "change_loglog"):
            self.assertPng(os.path.join(run, "plots", "mass", f + ".png"))
        self.assertEqual(sorted(os.listdir(os.path.join(run, "vtkOut"))),
                         ["solV_00000.vts", "solV_00005.vts"])
        self.assertTrue(os.path.isfile(os.path.join(run, "pf.pvd")))
        self.assertIn("Post-processing complete", r.stdout)

    def test_missing_outp_reports_error_but_makes_other_figures(self):
        # Without outp.txt the bounds step always fails, so the whole script
        # exits 1 even though every figure it could make was made.
        run = self.stage(outp=False)
        r = self.bash(run)
        self.assertEqual(r.returncode, 1, r.stdout[-2000:])
        self.assertIn("plot_phi_bounds.py exited 1", r.stdout)
        for f in ("porosity", "ssa", "mass"):
            self.assertPng(os.path.join(run, "plots", f + ".png"))
        self.assertFalse(os.path.exists(os.path.join(run, "plots", "timestep.png")))

    def test_only_outp(self):
        run = self.stage(ssa=False, iga=False)
        r = self.bash(run)
        self.assertEqual(r.returncode, 0, r.stdout[-2000:])
        self.assertIn("igasol.dat not found", r.stdout)
        self.assertIn("SSA_evo.dat not found", r.stdout)
        self.assertEqual(sorted(os.listdir(os.path.join(run, "plots"))),
                         ["phi_bounds.png", "timestep.png"])

    def test_no_postprocess_dir(self):
        run = fx.make_run(self.root)
        r = subprocess.run(["bash", os.path.join(fx.POSTPROCESS, "run_postprocess.sh"), run],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 1)
        self.assertIn("postprocess/ directory not found", r.stdout)


if __name__ == "__main__":
    unittest.main()
