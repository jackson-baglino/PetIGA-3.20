"""Tests for the plan tools T1-T8, figcheck, and the pplib functions added
with them (read_opts precedence, quadrature, gradients, sampler, snapshot
loader).

Each test builds a run whose answer is known exactly, or compares with an
independent quantity the solver wrote itself (SSA_evo.dat, bflux.dat). Tests
named test_edge_* probe the places where a tool is expected to break down;
their docstrings say what the tool does there and whether that is acceptable.

    venv_cryosuction/bin/python -m unittest discover -s postprocess/tests -t . -p 'test_new_tools.py'
"""
from __future__ import annotations

import csv
import io
import math
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

import numpy as np

from postprocess.tests import fixtures as fx

sys.path.insert(0, fx.POSTPROCESS)
import figcheck  # noqa: E402
import ice_bodies  # noqa: E402
import pplib  # noqa: E402
import section_flux  # noqa: E402
import type_summary  # noqa: E402
import water_budget  # noqa: E402

STUDY = os.path.join(fx.PROJECT, "studies", "cryosuction", "scripts")
VEL = os.path.join(fx.PROJECT, "studies", "contact_angle", "velocity_plan_2026-10-01")
RHO_ICE = 919.0


def quiet(fn, *a, **k):
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        return fn(*a, **k)


def disc_area(R, eps):
    """Exact integral of the model's tanh disc profile (R >> eps)."""
    return math.pi * (R * R + math.pi ** 2 * eps ** 2 / 3.0)


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


# ===========================================================================
# pplib: options
# ===========================================================================
class TestReadOpts(unittest.TestCase):
    def test_local_layout_experiment_beats_solver(self):
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d, geom="channel_A", exp="cryo_B", layout="local", tag="v2check")
            o = pplib.read_opts(run)
            self.assertEqual(o["-flag_BC_rhovfix"], "1")    # experiment, not solver.opts' 0
            self.assertEqual(o["-periodic"], "0")           # geometry, not solver.opts' 1

    def test_hpc_layout(self):
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d, geom="g1", exp="e1", layout="hpc")
            names = [os.path.basename(p) for p in pplib.opts_order(run)]
            self.assertEqual(names, ["solver.opts", "g1.opts", "e1.opts"])

    def test_extra_flags_beat_files(self):
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d)
            fx.write_outp(run, [0], [0.0], extra="-dtmax 1.0e2 -phi_open_bc 0 -snes_monitor -temp -30")
            o = pplib.read_opts(run)
            self.assertEqual(o["-dtmax"], "1.0e2")
            self.assertEqual(o["-phi_open_bc"], "0")
            self.assertEqual(o["-snes_monitor"], "")
            self.assertEqual(o["-temp"], "-30")              # a negative value is a value

    def test_edge_unrecognized_names_fall_back_to_t_final(self):
        """Folder renamed by hand: the file that sets -t_final is the experiment."""
        with tempfile.TemporaryDirectory() as d:
            run = os.path.join(d, "renamed")
            os.makedirs(run)
            fx.write_opts(os.path.join(run, "solver.opts"), {"-humidity": "1.0"})
            fx.write_opts(os.path.join(run, "zz_geom.opts"), {"-Nx": "10"})
            fx.write_opts(os.path.join(run, "aa_exp.opts"), {"-humidity": "0.9", "-t_final": "1"})
            self.assertEqual(pplib.read_opts(run)["-humidity"], "0.9")

    def test_edge_same_flag_twice_in_one_file_last_wins(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.opts")
            with open(p, "w") as fh:
                fh.write("-eps 1e-6\n-eps 2e-6   # override\n")
            self.assertEqual(pplib._parse_opts_file(p)["-eps"], "2e-6")

    def test_opt_bool_and_face_list(self):
        o = {"-a": "", "-b": "0", "-c": "FALSE", "-d": "yes", "-w": '"y0, Y1,x1"'}
        self.assertTrue(pplib.opt_bool(o, "-a"))
        self.assertFalse(pplib.opt_bool(o, "-b"))
        self.assertFalse(pplib.opt_bool(o, "-c"))
        self.assertTrue(pplib.opt_bool(o, "-d"))
        self.assertTrue(pplib.opt_bool(o, "-missing", True))
        self.assertEqual(pplib.face_list(o, "-w"), [(1, 0), (1, 1), (0, 1)])
        self.assertEqual(pplib.face_list({}, "-w"), [])


# ===========================================================================
# pplib: quadrature, gradients, sampling
# ===========================================================================
def mapped(GX, GY, Lx=100e-6, Ly=40e-6):
    """A smooth curvilinear map of the rectangle (rows not at constant y)."""
    return GX + 0.08 * Lx * np.sin(np.pi * GY / Ly) * np.sin(np.pi * GX / Lx), \
        GY + 0.10 * Ly * np.sin(np.pi * GX / Lx)


class TestQuadrature(unittest.TestCase):
    def test_node_weights_sum_to_area_rectangle(self):
        X, Y = fx.grid(3e-4, 1e-4, 37, 13)
        self.assertAlmostEqual(pplib.node_weights(X, Y).sum() / 3e-8, 1.0, places=12)

    def test_node_weights_trapezoid_exact_for_bilinear(self):
        X, Y = fx.grid(2.0, 1.0, 9, 5)
        f = 1 + 2 * X + 3 * Y + 4 * X * Y
        exact = 2 * 1 + 2 ** 2 * 1 + 3 * 2 * 0.5 + 4 * 2 * 0.5   # integral over [0,2]x[0,1]
        self.assertAlmostEqual((pplib.node_weights(X, Y) * f).sum(), exact, places=12)

    def test_node_weights_mapped_grid_converges(self):
        """Annulus sector r in [1, 2], theta in [0, pi/2]: area 3 pi / 4.
        Straight-sided cells cut the curved edges, an O(h^2) error."""
        errs = []
        for n in (11, 21, 41, 81):
            r, th = np.meshgrid(np.linspace(1, 2, n), np.linspace(0, np.pi / 2, n))
            errs.append(abs(pplib.node_weights(r * np.cos(th), r * np.sin(th)).sum() / (0.75 * np.pi) - 1))
        self.assertLess(errs[-1], 1e-3)
        rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
        np.testing.assert_allclose(rates[-2:], 2.0, atol=0.1)

    def test_gradient_exact_for_linear_on_mapped_grid(self):
        u, v = np.meshgrid(np.linspace(0, 100e-6, 41), np.linspace(0, 40e-6, 31))
        X, Y = mapped(u, v)
        F = 3.0 + 5.0 * X - 7.0 * Y
        gx, gy = pplib.grid_gradient(F, X, Y)
        np.testing.assert_allclose(gx, 5.0, rtol=1e-8)
        np.testing.assert_allclose(gy, -7.0, rtol=1e-8)

    def test_sampler_rect_and_mapped_and_outside(self):
        X, Y = fx.grid(1.0, 1.0, 11, 11)
        S = pplib.GridSampler(X, Y)
        self.assertTrue(S.rect)
        v = S(2 * X + Y, np.array([0.33, 1.5]), np.array([0.71, 0.5]))
        self.assertAlmostEqual(v[0], 2 * 0.33 + 0.71, places=12)
        self.assertTrue(np.isnan(v[1]))
        u, w = np.meshgrid(np.linspace(0, 100e-6, 31), np.linspace(0, 40e-6, 21))
        XM, YM = mapped(u, w)
        SM = pplib.GridSampler(XM, YM)
        self.assertFalse(SM.rect)
        px, py = np.array([50e-6]), np.array([20e-6])
        self.assertAlmostEqual(SM(XM - 2 * YM, px, py)[0], 50e-6 - 40e-6, delta=1e-12)


# ===========================================================================
# pplib: snapshot loader
# ===========================================================================
def lens_run(d, a_of_t, b_of_t, times, Lx=200e-6, Ly=40e-6, nx=200, ny=20, eps=0.5e-6, rhov=None,
             T=None, walls="y0,y1", every_snapshot_outp=True):
    """Igakit run of a plug a(t) < x < b(t) spanning a channel."""
    run = fx.make_run(d, geometry={"-eps": f"{eps:g}", "-Lx": f"{Lx:g}", "-Ly": f"{Ly:g}",
                                   "-wall_faces": walls})
    nrb = fx.iga_rectangle(Lx, Ly, nx, ny)
    snaps = []
    for k, t in enumerate(times):
        a, b = a_of_t(t), b_of_t(t)
        snaps.append((k, lambda X, Y, a=a, b=b: (fx.lens(X, a, b, eps),
                                                  (-20.0 + 0 * X) if T is None else T(X, Y),
                                                  (8e-4 + 0 * X) if rhov is None else rhov(X, Y))))
    fx.write_iga_run(run, nrb, snaps)
    if every_snapshot_outp:
        fx.write_outp(run, range(len(times)), times)
    return run


class TestIterSnapshots(unittest.TestCase):
    def test_time_from_ssa_when_outp_missing(self):
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(d, lambda t: 50e-6, lambda t: 100e-6, [0.0, 10.0], every_snapshot_outp=False)
            fx.write_ssa(run, [0.0, 10.0], [1.0, 1.0])
            ts = [s.t for s in pplib.iter_snapshots(run)]
            self.assertEqual(ts, [0.0, 10.0])

    def test_edge_unknown_time_is_skipped_not_guessed(self):
        """No outp.txt row and no SSA row: the snapshot is dropped with a note.
        plot_fields.py instead writes the STEP NUMBER as the time into pf.pvd
        in this case (plot_fields.py:271), so the .vts path can carry a wrong
        time silently; the nurbs path cannot."""
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(d, lambda t: 50e-6, lambda t: 100e-6, [0.0, 10.0], every_snapshot_outp=False)
            fx.write_ssa(run, [0.0], [1.0])           # step 1 has no time anywhere
            err = io.StringIO()
            with redirect_stderr(err):
                steps = [s.step for s in pplib.iter_snapshots(run)]
            self.assertEqual(steps, [0])
            self.assertIn("no time for step 1", err.getvalue())

    def test_restart_chain_later_leg_wins(self):
        with tempfile.TemporaryDirectory() as d:
            a = lens_run(os.path.join(d, "A"), lambda t: 50e-6, lambda t: 100e-6, [0.0, 10.0, 20.0])
            b = lens_run(os.path.join(d, "B"), lambda t: 60e-6, lambda t: 100e-6, [0.0, 10.0, 20.0])
            # leg B restarts at step 2 (t = 20) and goes on to 3, 4
            for f in os.listdir(b):
                if f.startswith("sol_00000") or f.startswith("sol_00001"):
                    os.remove(os.path.join(b, f))
            fx.write_iga_run(b, fx.iga_rectangle(200e-6, 40e-6, 200, 20),
                             [(3, lambda X, Y: (fx.lens(X, 61e-6, 100e-6, 0.5e-6), -20 + 0 * X, 8e-4 + 0 * X))])
            fx.write_outp(b, [2, 3], [20.0, 30.0])
            snaps = list(pplib.iter_snapshots([a, b]))
            self.assertEqual([s.step for s in snaps], [0, 1, 2, 3])
            self.assertTrue(snaps[2].path.startswith(b))

    def test_nurbs_and_vts_agree_on_area(self):
        with tempfile.TemporaryDirectory() as d:
            eps = 1e-6
            run = lens_run(d, lambda t: 60e-6, lambda t: 140e-6, [0.0], nx=100, ny=10, eps=eps)
            # Control-net values written as .vts, the way plot_fields.py does.
            nrb = fx.iga_rectangle(200e-6, 40e-6, 100, 10)
            GX, GY = nrb.points[..., 0].T, nrb.points[..., 1].T
            fx.write_vts_run(run, GX, GY, [(0, 0.0, fx.lens(GX, 60e-6, 140e-6, eps), -20 + 0 * GX,
                                            8e-4 + 0 * GX)])
            areas = {}
            for src in ("nurbs", "vts"):
                s = next(pplib.iter_snapshots(run, source=src))
                areas[src] = (pplib.node_weights(s.X, s.Y) * s.phi).sum()
            self.assertAlmostEqual(areas["nurbs"] / (80e-6 * 40e-6), 1.0, delta=1e-3)
            self.assertAlmostEqual(areas["vts"] / areas["nurbs"], 1.0, delta=5e-3)


# ===========================================================================
# T1 ice_bodies
# ===========================================================================
class TestIceBodies(unittest.TestCase):
    def test_moving_lens_faces_speed_and_acceleration(self):
        """Plug whose upper face moves at constant speed and lower face accelerates."""
        v0, acc = 2e-9, 1e-13
        times = [0.0, 1000.0, 2500.0, 3000.0, 4600.0]          # uneven on purpose
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(d, lambda t: 50e-6 + v0 * t, lambda t: 150e-6 + 0.5 * acc * t * t, times,
                           nx=400, ny=8, eps=0.5e-6)
            rows, contacts, events, _ = quiet(ice_bodies.analyze, [run])
            ice_bodies.add_rates(rows)
            self.assertEqual(len({r["body"] for r in rows}), 1)
            self.assertEqual(events, [])
            for r in rows:
                t = r["t_s"]
                self.assertAlmostEqual(r["face_lo"], 50e-6 + v0 * t, delta=2e-9)
                self.assertAlmostEqual(r["face_hi"], 150e-6 + 0.5 * acc * t * t, delta=2e-9)
            mid = rows[2]
            self.assertAlmostEqual(mid["v_face_lo"] / v0, 1.0, delta=0.01)
            # np.gradient on uneven times is exact for quadratics at interior points
            self.assertAlmostEqual(mid["v_face_hi"] / (acc * mid["t_s"]), 1.0, delta=0.01)
            self.assertAlmostEqual(mid["a_face_hi"] / acc, 1.0, delta=0.05)
            self.assertLess(abs(mid["a_face_lo"]), 0.05 * acc)
            # contact lines: each face meets both walls (y0, y1) at the face's x
            c0 = [c for c in contacts if c["step"] == 0]
            self.assertEqual(sorted(round(c["x"] * 1e6, 2) for c in c0), [50.0, 50.0, 150.0, 150.0])
            self.assertEqual({c["wall"] for c in c0}, {"y0", "y1"})

    def test_area_and_com_of_two_discs(self):
        with tempfile.TemporaryDirectory() as d:
            eps, Lx, Ly = 0.5e-6, 120e-6, 60e-6
            run = fx.make_run(d, geometry={"-eps": f"{eps:g}"})
            nrb = fx.iga_rectangle(Lx, Ly, 240, 120)
            fx.write_iga_run(run, nrb, [(0, lambda X, Y: (fx.union(fx.disc(X, Y, 30e-6, 30e-6, 15e-6, eps),
                                                                    fx.disc(X, Y, 85e-6, 25e-6, 10e-6, eps)),
                                                          -20 + 0 * X, 8e-4 + 0 * X))])
            fx.write_outp(run, [0], [0.0])
            rows, _, _, _ = quiet(ice_bodies.analyze, [run])
            self.assertEqual(len(rows), 2)
            big, small = sorted(rows, key=lambda r: -r["area_m2"])
            self.assertAlmostEqual(big["area_m2"] / disc_area(15e-6, eps), 1.0, delta=2e-3)
            self.assertAlmostEqual(small["area_m2"] / disc_area(10e-6, eps), 1.0, delta=2e-3)
            self.assertAlmostEqual(big["com_x"], 30e-6, delta=2e-8)
            self.assertAlmostEqual(small["com_y"], 25e-6, delta=2e-8)
            self.assertAlmostEqual(big["area_core_m2"] / (math.pi * 15e-6 ** 2), 1.0, delta=5e-3)

    def test_merge_and_split_events_and_new_ids(self):
        eps = 0.5e-6
        shapes = [  # two discs -> touching (merged) -> two again
            lambda X, Y: fx.union(fx.disc(X, Y, 35e-6, 30e-6, 12e-6, eps), fx.disc(X, Y, 85e-6, 30e-6, 12e-6, eps)),
            lambda X, Y: fx.union(fx.disc(X, Y, 45e-6, 30e-6, 12e-6, eps), fx.disc(X, Y, 69e-6, 30e-6, 12e-6, eps),
                                  fx.lens(X, 45e-6, 69e-6, eps) * (np.abs(Y - 30e-6) < 4e-6)),
            lambda X, Y: fx.union(fx.disc(X, Y, 35e-6, 30e-6, 12e-6, eps), fx.disc(X, Y, 85e-6, 30e-6, 12e-6, eps)),
        ]
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d, geometry={"-eps": f"{eps:g}"})
            nrb = fx.iga_rectangle(120e-6, 60e-6, 120, 60)
            fx.write_iga_run(run, nrb, [(k, lambda X, Y, f=f: (f(X, Y), -20 + 0 * X, 8e-4 + 0 * X))
                                        for k, f in enumerate(shapes)])
            fx.write_outp(run, [0, 1, 2], [0.0, 1.0, 2.0])
            rows, _, events, _ = quiet(ice_bodies.analyze, [run])
            kinds = [e[2] for e in events]
            self.assertIn("merge", kinds)
            self.assertIn("split", kinds)
            ids0 = {r["body"] for r in rows if r["step"] == 0}
            ids1 = {r["body"] for r in rows if r["step"] == 1}
            self.assertEqual(len(ids1), 1)
            self.assertFalse(ids1 & ids0)          # a merged body never inherits an old id

    def test_edge_corner_touching_cores_are_two_bodies(self):
        """4-connectivity: cores that share only a grid corner stay separate."""
        phi = np.zeros((6, 6))
        phi[1:3, 1:3] = 1.0
        phi[3:5, 3:5] = 1.0
        lab, n = ice_bodies.label_cores(phi)
        self.assertEqual(n, 2)

    def test_edge_tiny_core_below_min_area_is_not_a_body(self):
        with tempfile.TemporaryDirectory() as d:
            eps = 0.5e-6
            run = fx.make_run(d, geometry={"-eps": f"{eps:g}"})
            nrb = fx.iga_rectangle(100e-6, 50e-6, 200, 100)
            fx.write_iga_run(run, nrb, [(0, lambda X, Y: (fx.union(fx.disc(X, Y, 30e-6, 25e-6, 12e-6, eps),
                                                                    fx.disc(X, Y, 80e-6, 25e-6, 0.6e-6, eps)),
                                                          -20 + 0 * X, 8e-4 + 0 * X))])
            fx.write_outp(run, [0], [0.0])
            rows, _, _, _ = quiet(ice_bodies.analyze, [run])
            self.assertEqual(len(rows), 1)
            # its diffuse ice is given to the remaining body, so the total is kept
            s = next(pplib.iter_snapshots(run))
            total = (pplib.node_weights(s.X, s.Y) * s.phi).sum()
            self.assertAlmostEqual(rows[0]["area_m2"] / total, 1.0, places=10)

    def test_edge_body_against_domain_end_has_one_face(self):
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(d, lambda t: -10e-6, lambda t: 80e-6, [0.0])
            rows, _, _, _ = quiet(ice_bodies.analyze, [run])
            self.assertTrue(np.isnan(rows[0]["face_lo"]))
            self.assertAlmostEqual(rows[0]["face_hi"], 80e-6, delta=2e-9)

    def test_edge_single_snapshot_rates_are_nan(self):
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(d, lambda t: 50e-6, lambda t: 100e-6, [0.0])
            rows, _, _, _ = quiet(ice_bodies.analyze, [run])
            ice_bodies.add_rates(rows)
            self.assertTrue(np.isnan(rows[0]["v_face_lo"]))
            self.assertTrue(np.isnan(rows[0]["a_face_lo"]))

    def test_edge_rates_never_bridge_a_gap(self):
        """A face that is missing in one snapshot (NaN) splits the series; no
        velocity is formed across the gap."""
        rows = [{"body": 0, "t_s": float(t), "face_lo": f, "com_x": 0.0, "com_y": 0.0, "face_hi": 0.0,
                 "area_m2": 1.0} for t, f in enumerate([1.0, 2.0, np.nan, 4.0, 5.0, 6.0])]
        ice_bodies.add_rates(rows)
        v = [r["v_face_lo"] for r in rows]
        self.assertTrue(np.isnan(v[2]))
        self.assertAlmostEqual(v[0], 1.0)
        self.assertAlmostEqual(v[4], 1.0)
        self.assertTrue(np.isnan(rows[0]["a_face_lo"]))           # only two points before the gap
        self.assertAlmostEqual(rows[4]["a_face_lo"], 0.0)

    def test_sample_at_interpolates_and_reports_never(self):
        rows = [{"body": 0, "step": k, "t_s": 10.0 * k, "face_lo": 1e-5 * k, "gap": np.nan}
                for k in range(5)]
        rows += [{"body": 1, "step": k, "t_s": 10.0 * k, "face_lo": 5e-5, "gap": np.nan} for k in range(5)]
        out = ice_bodies.sample_at(rows, "face_lo", 2.5e-5)
        self.assertEqual(len(out), 1)                 # body 1 never crosses: no row
        self.assertAlmostEqual(out[0]["t_s"], 25.0)
        # exactly on a sample: one row, not two
        self.assertEqual(len(ice_bodies.sample_at(rows, "face_lo", 2e-5)), 1)

    def test_mapped_mesh_axis_line(self):
        """Rows of a mapped mesh are not at constant y; the axis line is
        interpolated column by column, so the faces are still exact."""
        def mp(GX, GY):
            return GX, GY * (1.0 + 0.3 * np.sin(np.pi * GX / 200e-6))
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d, geometry={"-eps": "0.5e-6", "-wall_faces": "y0,y1"})
            nrb = fx.iga_rectangle(200e-6, 40e-6, 400, 16, mapping=mp)
            fx.write_iga_run(run, nrb, [(0, lambda X, Y: (fx.lens(X, 60e-6, 130e-6, 0.5e-6), -20 + 0 * X,
                                                          8e-4 + 0 * X))])
            fx.write_outp(run, [0], [0.0])
            rows, contacts, _, _ = quiet(ice_bodies.analyze, [run])
            self.assertAlmostEqual(rows[0]["face_lo"], 60e-6, delta=5e-9)
            self.assertAlmostEqual(rows[0]["face_hi"], 130e-6, delta=5e-9)
            top = sorted(c["x"] for c in contacts if c["wall"] == "y1")
            np.testing.assert_allclose(top, [60e-6, 130e-6], atol=5e-9)

    def test_cli_writes_outputs_and_figure_is_clean(self):
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(d, lambda t: 50e-6 + 1e-9 * t, lambda t: 150e-6, [0.0, 1e3, 2e3, 3e3], nx=200, ny=8)
            out = os.path.join(d, "out")
            res = subprocess.run([fx.PY, os.path.join(fx.POSTPROCESS, "ice_bodies.py"), "--dir", run,
                                  "--out", out, "--plot", "--vertical", "--at", "face_lo=5.15e-5"],
                                 capture_output=True, text=True, cwd=fx.PROJECT)
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertNotIn("figure check", res.stderr)
            at = read_rows(os.path.join(out, "ice_bodies_at.csv"))
            self.assertAlmostEqual(float(at[0]["t_s"]), 1500.0, delta=5.0)
            self.assertTrue(os.path.isfile(os.path.join(out, "plots", "ice_bodies.png")))


# ===========================================================================
# T2 water_budget
# ===========================================================================
class TestWaterBudget(unittest.TestCase):
    def test_exact_balance_closes(self):
        t = np.linspace(0, 1e4, 101)
        ice = 2e-8 - 1e-13 * t                        # m^2
        vap = 2e-11 + 1e-16 * t                       # kg/m
        J = RHO_ICE * 1e-13 - 1e-16                   # outflow that balances exactly
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d)
            fx.write_ssa(run, t, ice, rhov=vap)
            fx.write_bflux(run, t, {"x0": (J * np.ones_like(t), -2.8e6 * J * np.ones_like(t))})
            b = water_budget.budget([run])
            self.assertLess(np.abs(b["closure"]).max() / b["lost"][-1], 1e-6)
            self.assertAlmostEqual(b["heat"][-1] / b["out"][-1], 2.8e6, delta=1.0)

    def test_closed_box_without_bflux(self):
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d)
            t = np.linspace(0, 100, 11)
            fx.write_ssa(run, t, 1e-8 * np.ones_like(t), rhov=1e-12)
            b = water_budget.budget([run])
            self.assertFalse(b["has_bflux"])
            self.assertTrue(np.all(b["out"] == 0))
            self.assertLess(np.abs(b["closure"]).max(), 1e-20)

    def test_edge_leak_is_reported_as_closure(self):
        """Ice lost with no outflow and no vapor gain: the whole loss is closure error."""
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d)
            t = np.linspace(0, 100, 11)
            fx.write_ssa(run, t, 1e-8 - 1e-11 * t, rhov=0.0)
            fx.write_bflux(run, t, {"x0": (0.0, 0.0)})
            b = water_budget.budget([run])
            self.assertAlmostEqual(abs(b["closure"][-1]) / b["lost"][-1], 1.0, places=6)

    def test_edge_legacy_ssa_without_vapor_column_fails_loudly(self):
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d)
            fx.write_ssa(run, [0, 1, 2], [1e-8, 1e-8, 1e-8], ncol=5)
            with self.assertRaises(SystemExit):
                water_budget.budget([run])

    def test_edge_bflux_starting_late_warns(self):
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d)
            t = np.linspace(0, 100, 11)
            fx.write_ssa(run, t, 1e-8 * np.ones_like(t))
            fx.write_bflux(run, t[3:], {"x0": (0.0, 0.0)})
            err = io.StringIO()
            with redirect_stderr(err):
                water_budget.budget([run])
            self.assertIn("starts at", err.getvalue())

    def test_edge_restart_duplicate_steps_keep_last(self):
        with tempfile.TemporaryDirectory() as d:
            run = fx.make_run(d)
            fx.write_ssa(run, [0, 1, 2], [1e-8] * 3)
            fx.write_bflux(run, [0.0, 1.0, 2.0, 1.0, 2.0], {"x0": ([9.0, 9.0, 9.0, 1.0, 1.0], 0.0)},
                           steps=[0, 1, 2, 1, 2])
            t, J, Q, faces = water_budget.load_bflux(run)
            np.testing.assert_array_equal(t, [0.0, 1.0, 2.0])
            np.testing.assert_array_equal(J, [9.0, 1.0, 1.0])

    def test_restart_chain_concatenates(self):
        with tempfile.TemporaryDirectory() as d:
            a = fx.make_run(os.path.join(d, "a"))
            b = fx.make_run(os.path.join(d, "b"))
            fx.write_ssa(a, [0, 1, 2], [3e-8, 2e-8, 1e-8], rhov=0.0)
            fx.write_ssa(b, [2, 3], [1e-8, 0.0], rhov=0.0, steps=[2, 3])
            J = RHO_ICE * 1e-8
            fx.write_bflux(a, [0, 1, 2], {"x0": (J, 0.0)})
            fx.write_bflux(b, [2, 3], {"x0": (J, 0.0)}, steps=[2, 3])
            bud = water_budget.budget([a, b])
            np.testing.assert_array_equal(bud["t"], [0, 1, 2, 3])
            self.assertLess(np.abs(bud["closure"]).max() / bud["lost"][-1], 1e-6)

    def test_recapture_translating_lens_is_one_and_shrinking_is_zero(self):
        times = [0.0, 1.0, 2.0, 3.0]
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(os.path.join(d, "move"), lambda t: 50e-6 + 5e-6 * t, lambda t: 100e-6 + 5e-6 * t, times)
            sub, dep, n = water_budget.recapture([run], "auto", 1)
            self.assertEqual(n, 4)
            self.assertAlmostEqual(dep / sub, 1.0, delta=1e-3)
            self.assertAlmostEqual(sub / RHO_ICE, 3 * 5e-6 * 40e-6, delta=0.01 * 3 * 5e-6 * 40e-6)
            run2 = lens_run(os.path.join(d, "shrink"), lambda t: 50e-6 + 5e-6 * t, lambda t: 150e-6, times)
            sub2, dep2, _ = water_budget.recapture([run2], "auto", 1)
            self.assertLess(dep2 / sub2, 1e-3)

    def test_edge_recapture_misses_reversal_between_snapshots(self):
        """A face that recedes and returns between two snapshots is invisible:
        recapture is a lower bound that depends on cadence (documented)."""
        times = [0.0, 1.0, 2.0]
        with tempfile.TemporaryDirectory() as d:
            run = lens_run(d, lambda t: 50e-6 + 5e-6 * min(t, 2 - t), lambda t: 100e-6, times)
            sub1, dep1, _ = water_budget.recapture([run], "auto", 1)
            sub2, dep2, _ = water_budget.recapture([run], "auto", 2)
            self.assertAlmostEqual(dep1 / sub1, 1.0, delta=1e-3)
            self.assertLess(sub2, 1e-3 * sub1)           # at half cadence nothing happened


# ===========================================================================
# T3 section_flux
# ===========================================================================
class TestSectionFlux(unittest.TestCase):
    def run_linear(self, d, phi_fn=None, mapping=None, T_fn=None):
        Lx, Ly, g = 100e-6, 40e-6, 2.0        # rho_v = 8e-4 + g x  (kg/m^3)
        run = fx.make_run(d)
        nrb = fx.iga_rectangle(Lx, Ly, 50, 20, mapping=mapping)
        fx.write_iga_run(run, nrb, [(0, lambda X, Y: ((0 * X) if phi_fn is None else phi_fn(X, Y),
                                                      (-20 + 0 * X) if T_fn is None else T_fn(X, Y),
                                                      8e-4 + g * X))])
        fx.write_outp(run, [0], [0.0])
        return run, g, Ly

    def test_linear_vapor_field_exact(self):
        with tempfile.TemporaryDirectory() as d:
            run, g, Ly = self.run_linear(d)
            s = next(pplib.iter_snapshots(run))
            J, Q, li, lt = section_flux.line_flux(s, (50e-6, 0.0), (50e-6, Ly))
            exact = -pplib.dif_vap(-20.0) * 1.0 * g * Ly       # normal +x, flow toward -x
            self.assertAlmostEqual(J / exact, 1.0, delta=1e-9)
            self.assertAlmostEqual(li, lt)
            # reversing the line flips the normal and the sign
            J2, *_ = section_flux.line_flux(s, (50e-6, Ly), (50e-6, 0.0))
            self.assertAlmostEqual(J2 / J, -1.0, delta=1e-12)

    def test_heat_flux_linear_temperature(self):
        with tempfile.TemporaryDirectory() as d:
            G = 100.0
            run, g, Ly = self.run_linear(d, T_fn=lambda X, Y: -20 + G * (X - 50e-6))
            s = next(pplib.iter_snapshots(run))
            J, Q, *_ = section_flux.line_flux(s, (30e-6, 0.0), (30e-6, Ly))
            self.assertAlmostEqual(Q / (-pplib.K_AIR * G * Ly), 1.0, delta=1e-9)
            # D_v(T) is evaluated at the local T, not at temp0
            self.assertAlmostEqual(J / (-pplib.dif_vap(-20 + G * (30e-6 - 50e-6)) * g * Ly), 1.0, delta=1e-9)

    def test_ice_blocks_vapor_flux(self):
        with tempfile.TemporaryDirectory() as d:
            run, g, Ly = self.run_linear(d, phi_fn=lambda X, Y: 1.0 + 0 * X)
            s = next(pplib.iter_snapshots(run))
            J, Q, *_ = section_flux.line_flux(s, (50e-6, 0.0), (50e-6, Ly))
            self.assertAlmostEqual(J / (-pplib.dif_vap(-20.0) * pplib.AIR_LIM * g * Ly), 1.0, delta=1e-6)

    def test_mapped_mesh_linear_field(self):
        with tempfile.TemporaryDirectory() as d:
            run, g, Ly = self.run_linear(d, mapping=lambda GX, GY: (GX + 0.1 * GY, GY))
            s = next(pplib.iter_snapshots(run))
            # a vertical line through the sheared domain, fully inside it
            J, _, li, lt = section_flux.line_flux(s, (60e-6, 1e-6), (60e-6, 39e-6))
            self.assertAlmostEqual(J / (-pplib.dif_vap(-20.0) * g * 38e-6), 1.0, delta=1e-6)

    def test_edge_line_partly_outside_reports_length(self):
        with tempfile.TemporaryDirectory() as d:
            run, g, Ly = self.run_linear(d)
            s = next(pplib.iter_snapshots(run))
            J, _, li, lt = section_flux.line_flux(s, (50e-6, -20e-6), (50e-6, Ly))
            self.assertAlmostEqual(li / lt, Ly / (Ly + 20e-6), delta=0.01)
            self.assertAlmostEqual(J / (-pplib.dif_vap(-20.0) * g * Ly), 1.0, delta=0.01)

    def test_edge_zero_length_line_raises(self):
        with tempfile.TemporaryDirectory() as d:
            run, g, Ly = self.run_linear(d)
            s = next(pplib.iter_snapshots(run))
            with self.assertRaises(ValueError):
                section_flux.line_flux(s, (1e-6, 1e-6), (1e-6, 1e-6))

    def test_storage_between_two_lines(self):
        """Conservation in the air: what enters between two lines minus what
        leaves equals the vapor stored there (zero in steady linear profile)."""
        with tempfile.TemporaryDirectory() as d:
            run, g, Ly = self.run_linear(d)
            s = next(pplib.iter_snapshots(run))
            J1, *_ = section_flux.line_flux(s, (20e-6, 0.0), (20e-6, Ly))
            J2, *_ = section_flux.line_flux(s, (80e-6, 0.0), (80e-6, Ly))
            self.assertAlmostEqual(J1 / J2, 1.0, delta=1e-9)


# ===========================================================================
# T4 rhov_bc and T5 check_open_bc (scripts)
# ===========================================================================
def rho_vs_indep(T_C):
    """ASHRAE ice saturation density, written out again from
    src/material_properties.c RhoVS_I (independent of preprocess/ and pplib)."""
    T = T_C + 273.15
    k = (-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918)
    P = math.exp(k[0] / T + k[1] + k[2] * T + k[3] * T * T + k[4] * T ** 3 + k[5] * math.log(T))
    return 1.341 * 0.62 * P / (101325.0 - P)


def run_script(path, *args):
    cmd = ([fx.PY, path] if path.endswith(".py") else ["bash", path]) + [str(a) for a in args]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=fx.PROJECT)


def opt_value(out, flag):
    for line in out.splitlines():
        if line.startswith(flag + " "):
            return float(line.split()[1])
    raise AssertionError(f"{flag} not printed:\n{out}")


class TestRhovBC(unittest.TestCase):
    S = os.path.join(STUDY, "rhov_bc.py")

    def test_local_saturation_pair(self):
        r = run_script(self.S, "--T0", "-20", "--G", "100", "--Lx", "2e-3", "--local-sat")
        self.assertEqual(r.returncode, 0, r.stderr)
        lo, hi = opt_value(r.stdout, "-rhovfix_lo"), opt_value(r.stdout, "-rhovfix_hi")
        self.assertAlmostEqual(lo, rho_vs_indep(-20.1) / rho_vs_indep(-20.0), places=9)
        self.assertAlmostEqual(hi, rho_vs_indep(-19.9) / rho_vs_indep(-20.0), places=9)

    def test_gibbs_thomson_factor(self):
        r = run_script(self.S, "--T0", "-20", "--h-open", "0.99", "--chi", "1e5")
        d0 = 0.109 * (2.99e-26 / 919.0) / (1.38e-23 * 253.15)   # KP Eq. 13, rough constants
        self.assertAlmostEqual(opt_value(r.stdout, "-rhovfix_lo") / 0.99 - 1, d0 * 1e5, delta=0.02 * d0 * 1e5)

    def test_no_gradient_gives_one(self):
        r = run_script(self.S, "--T0", "-20", "--local-sat")
        self.assertEqual(opt_value(r.stdout, "-rhovfix_lo"), 1.0)

    def test_edge_inputs_rejected(self):
        for args in (["--T0", "-20", "--G", "100", "--local-sat"],     # G without Lx
                     ["--T0", "-300", "--local-sat"],                  # below 0 K
                     ["--T0", "-20", "--h-face", "-1"],
                     ["--T0", "-20", "--local-sat", "--h-open", "0.9"]):
            r = run_script(self.S, *args)
            self.assertEqual(r.returncode, 2, args)
            self.assertNotIn("-rhovfix_lo", r.stdout)

    def test_edge_lunar_temperature_warns(self):
        r = run_script(self.S, "--T0", "-173", "--local-sat")
        self.assertIn("WARNING", r.stdout)


class TestCheckOpenBC(unittest.TestCase):
    S = os.path.join(STUDY, "check_open_bc.sh")
    G = "channel_2D_H125um_eps0.39um_deadend"
    E = "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-3"     # sets -phi_open_bc 0 and pins x0

    def check(self, *extra):
        args = [self.G, self.E] + (["--", *extra] if extra else [])
        return run_script(self.S, *args)

    def test_file_as_written_passes(self):
        r = self.check()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_open_face_with_open_phi_fails(self):
        self.assertEqual(self.check("-phi_open_bc", "1").returncode, 1)
        self.assertEqual(self.check("-phi_open_bc", "true").returncode, 1)

    def test_edge_bare_flag_means_true(self):
        """PETSc reads a bare boolean flag as true."""
        self.assertEqual(self.check("-phi_open_bc").returncode, 1)

    def test_no_open_face_cases_pass(self):
        for extra in (["-phi_open_bc", "1", "-flag_BC_rhovfix", "0"],    # no reservoir
                      ["-phi_open_bc", "1", "-periodic", "1"],           # solver drops rhovfix
                      ["-phi_open_bc", "1", "-wall_faces", "x0,y0,y1,x1"]):  # reservoir is a wall
            r = self.check(*extra)
            self.assertEqual(r.returncode, 0, (extra, r.stdout))

    def test_edge_axisym_axis_is_not_a_reservoir(self):
        # (the channel geometry lists y0 as a wall; make x1 the only wall here)
        r = self.check("-phi_open_bc", "1", "-rhovfix_faces", "y0", "-wall_faces", "x1", "-axisym", "1")
        self.assertEqual(r.returncode, 0, r.stdout)
        r = self.check("-phi_open_bc", "1", "-rhovfix_faces", "y0", "-wall_faces", "x1")
        self.assertEqual(r.returncode, 1, r.stdout)

    def test_edge_tests_file_bad_lines_and_missing_names(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.txt")
            with open(p, "w") as fh:
                fh.write(f"# comment\n\n{self.G}:{self.E}   # ok\nnot_a_pair\n{self.G}:no_such_experiment\n")
            r = run_script(self.S, p)
            self.assertEqual(r.returncode, 1)
            self.assertIn("not a geometry:experiment line", r.stdout)
            self.assertIn("experiment not found", r.stdout)
            self.assertEqual(r.stdout.count("OK "), 1)

    def test_edge_usage(self):
        self.assertEqual(run_script(self.S).returncode, 2)
        self.assertEqual(run_script(self.S, "/no/such/file").returncode, 2)


# ===========================================================================
# T6 type_summary, T8 compare_velocity, figcheck
# ===========================================================================
class TestTypeSummary(unittest.TestCase):
    def test_time_to_lose(self):
        t = np.array([0.0, 10.0, 20.0])
        lost = np.array([0.0, 4.0, 8.0])
        self.assertAlmostEqual(type_summary.time_to_lose(t, lost, 10.0, 0.5), 12.5)
        self.assertTrue(np.isnan(type_summary.time_to_lose(t, lost, 10.0, 0.9)))

    def test_edge_manifest_errors(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.txt")
            for text in ("4 only_two\n", "4 A /no/such/folder\n", "# nothing\n"):
                with open(p, "w") as fh:
                    fh.write(text)
                with self.assertRaises(SystemExit):
                    type_summary.read_manifest(p)

    def test_end_to_end_synthetic(self):
        times = [0.0, 1e3, 2e3, 3e3]
        with tempfile.TemporaryDirectory() as d:
            runs = []
            for name, v in (("A", 2e-9), ("B", 4e-9)):
                run = lens_run(os.path.join(d, name), lambda t, v=v: 50e-6 + v * t, lambda t: 150e-6, times,
                               nx=200, ny=8)
                ice = (150e-6 - (50e-6 + v * np.array(times))) * 40e-6
                fx.write_ssa(run, times, ice, rhov=0.0)
                fx.write_bflux(run, times, {"x0": (RHO_ICE * v * 40e-6, 0.0)})
                runs.append((name, run))
            m = os.path.join(d, "m.txt")
            with open(m, "w") as fh:
                fh.writelines(f"7 {n} {r}\n" for n, r in runs)
            out = os.path.join(d, "out")
            res = subprocess.run([fx.PY, os.path.join(fx.POSTPROCESS, "type_summary.py"), "--manifest", m,
                                  "--out", out, "--vertical"], capture_output=True, text=True, cwd=fx.PROJECT)
            self.assertEqual(res.returncode, 0, res.stdout + res.stderr)
            rows = read_rows(os.path.join(out, "type7_summary.csv"))
            self.assertEqual([r["name"] for r in rows], ["A", "B"])
            for r, v in zip(rows, (2e-9, 4e-9)):
                self.assertAlmostEqual(float(r["v_face_lo_end_m_per_s"]) / v, 1.0, delta=0.01)
                self.assertLess(float(r["closure_rel"]), 1e-5)
            # the run folders were not written into
            self.assertFalse(os.path.exists(os.path.join(runs[0][1], "ice_bodies.csv")))


class TestCompareVelocityAtPosition(unittest.TestCase):
    def setUp(self):
        sys.path.insert(0, VEL)
        import compare_velocity
        self.cv = compare_velocity

    def tearDown(self):
        sys.path.remove(VEL)

    def run_dict(self, ell, v):
        t = np.linspace(0, 10, 11) * 86400.0
        return {"geom": "channel", "t": t, "ell": ell(t), "v": {"mean": v(t)}, "v_th": {"mean": 2 * v(t)},
                "v_thin": {"mean": 3 * v(t)}}

    def test_interpolates_at_first_crossing(self):
        r = self.run_dict(lambda t: 1e-4 + 1e-6 * t / 86400.0, lambda t: t / 86400.0)
        tc, vm, vt, vn = self.cv.at_position(r, "mean", 1.045e-4)
        self.assertAlmostEqual(tc / 86400.0, 4.5)
        self.assertAlmostEqual(vm, 4.5)
        self.assertAlmostEqual(vt, 9.0)

    def test_edge_never_reached_is_nan(self):
        r = self.run_dict(lambda t: 1e-4 - 1e-6 * t / 86400.0, lambda t: 0 * t + 1.0)
        self.assertTrue(np.isnan(self.cv.at_position(r, "mean", 2e-4)[1]))

    def test_edge_crossing_inside_first_day_is_ignored(self):
        r = self.run_dict(lambda t: 1e-4 + 1e-6 * t / 86400.0, lambda t: 0 * t + 1.0)
        self.assertTrue(np.isnan(self.cv.at_position(r, "mean", 1.005e-4)[0]))


class TestFigcheck(unittest.TestCase):
    def fig(self, **kw):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt.subplots(figsize=(4, 3), **kw)

    def test_clean_figure_passes(self):
        fig, ax = self.fig()
        ax.plot([0, 1], [0, 1], label="a")
        ax.set_xlabel("x")
        ax.legend(loc="lower right")
        self.assertEqual(figcheck.check_figure(fig), [])

    def test_overlaps_and_small_font_are_caught(self):
        fig, ax = self.fig()
        ax.plot([0, 1], [0, 1])
        ax.text(0.5, 0.5, "on the line")
        ax.text(0.52, 0.51, "on the text")
        ax.set_xlabel("tiny", fontsize=5)
        probs = figcheck.check_figure(fig)
        self.assertTrue(any("text/text" in p for p in probs))
        self.assertTrue(any("line" in p for p in probs))
        self.assertTrue(any("smallest printed font" in p for p in probs))

    def test_placed_width_scales_font(self):
        fig, ax = self.fig()
        ax.set_xlabel("x", fontsize=10)
        ax.tick_params(labelsize=10)
        self.assertEqual(figcheck.check_figure(fig, placed_width_in=4.0), [])
        self.assertTrue(figcheck.check_figure(fig, placed_width_in=2.0))    # 10 pt -> 5 pt

    def test_colorbar_and_overlap_ok(self):
        fig, ax = self.fig()
        im = ax.imshow(np.zeros((4, 4)))
        fig.colorbar(im, ax=ax)
        t = ax.text(1.5, 1.5, "inside image")
        self.assertTrue(any("image" in p for p in figcheck.check_figure(fig)))
        t.set_gid("overlap_ok")
        self.assertEqual(figcheck.check_figure(fig), [])

    def test_edge_crowded_panels_are_caught(self):
        """Three panels squeezed into 4 inches: neighboring tick labels collide."""
        fig, axs = self.fig(ncols=3)
        for a in axs:
            a.plot([0, 5.5], [0, -0.0009])
        fig.tight_layout()
        self.assertTrue(any("text/text" in p for p in figcheck.check_figure(fig)))

    def test_edge_hidden_tick_labels_ignored(self):
        """Ticks outside the view keep invisible labels; they must not count."""
        fig, axs = self.fig(ncols=3)
        fig.set_size_inches(9, 3)
        for a in axs:
            a.plot([0, 5.5], [0, -0.0009])
        fig.tight_layout()
        self.assertEqual(figcheck.check_figure(fig), [])


# ===========================================================================
# Integration on the real stage-0c run (skipped when it is not on this machine)
# ===========================================================================
@unittest.skipUnless(fx.real_run_available(), "real run not present")
class TestRealRun(unittest.TestCase):
    R = fx.REAL_RUN

    def test_body_area_matches_solver_tot_ice(self):
        a = pplib.load_ssa(self.R)
        rows, _, _, _ = quiet(ice_bodies.analyze, [self.R], every=10)
        by_step = {int(r[pplib.STEP]): r[pplib.ICE] for r in a}
        for r in rows:
            if r["step"] in by_step:
                # SSA_evo.dat carries 7 significant digits
                self.assertAlmostEqual(r["area_m2"] / by_step[r["step"]], 1.0, delta=2e-6)

    def test_face_flux_matches_bflux(self):
        rows = section_flux.run([self.R], [lambda s: section_flux.face_line(s, "x0")], every=5)
        ref = section_flux.bflux_at(self.R, "x0", [r[0] for r in rows])
        n = 0
        for st, _t, _k, J, Q, _li, _lt in rows:
            jb, qb = ref[st]
            if jb != 0:
                self.assertAlmostEqual(J / jb, 1.0, delta=1e-6)
                self.assertAlmostEqual(Q / qb, 1.0, delta=1e-4)
                n += 1
        self.assertGreater(n, 2)

    def test_water_budget_closes(self):
        b = water_budget.budget([self.R])
        self.assertLess(np.abs(b["closure"]).max() / b["lost"][-1], 2e-3)
        # isothermal run: the heat drawn in is the latent heat of what left
        self.assertAlmostEqual(b["heat"][-1] / b["out"][-1] / 2.83e6, 1.0, delta=0.01)

    def test_read_opts_on_real_run(self):
        o = pplib.read_opts(self.R)
        self.assertEqual(o["-flag_BC_rhovfix"], "1")
        self.assertEqual(o["-periodic"], "0")
        self.assertEqual(o["-phi_open_bc"], "0")


if __name__ == "__main__":
    unittest.main()
