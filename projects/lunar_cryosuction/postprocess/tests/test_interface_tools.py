"""Known-answer and edge-case tests for the interface-measurement tools.

Tools under test (postprocess/):
    contact_angle.py       macroscopic contact angle at a regolith wall
    meniscus_velocity.py   normal speed of the two menisci of a channel bridge
    wedge_gt_velocity.py   Gibbs-Thomson check on the wedge-band centerline
    gt_balance.py          s = d0*chi + beta*v_n on the phi = 0.5 contour (pvpython)

Every case is built from an analytic shape with the model profile
phi = (1 + tanh(d/(2 eps)))/2, so the right answer is known before the tool
runs. Where a test exposes a bug it asserts the CORRECT behavior and is marked
expectedFailure, with a "# BUG:" comment naming file:line.

Run from the project root:
    venv_cryosuction/bin/python -m unittest discover -s postprocess/tests -t . \
        -p 'test_interface_tools.py'
"""
from __future__ import annotations

import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from postprocess.tests import fixtures as fx

sys.path.insert(0, fx.POSTPROCESS)
import contact_angle as ca            # noqa: E402
import gt_balance as gtb              # noqa: E402
import meniscus_velocity as mv        # noqa: E402
import wedge_gt_velocity as wg        # noqa: E402
from pplib import rho_vs              # noqa: E402

PVPYTHON = "/Applications/ParaView-6.1.1.app/Contents/bin/pvpython"
T0 = -20.0
D0 = 1.0152e-9          # d0_sub0 the solver prints at -20 C (K&P Eq. 13)
BETA0 = 7.9408e5        # -beta_sub0 at -20 C, alpha_c = 1e-2
RHO_ICE = 919.0
SQRT2 = float(np.sqrt(2.0))   # eps/h of the production meshes (h = eps/sqrt 2)


# ---------------------------------------------------------------------------
# Helpers (kept here; fixtures.py is shared and not modified)
# ---------------------------------------------------------------------------
def bridge_distance(X, Y, a, b, H, theta_deg):
    """Signed distance into the ice of a bridge whose two menisci are circular
    arcs meeting the flat walls y = 0 and y = H at contact angle theta.

    The angle is the solver's: cos(theta) = m . n, with m pointing into the ice
    and n the outward wall normal. theta < 90: ice outside circles centered in
    the vapor (concave, wetting); theta > 90: ice inside circles (convex).
    a and b are the mid-plane positions of the left and right menisci.
    R = H / (2 |cos theta|).
    """
    c = np.cos(np.radians(theta_deg))
    if abs(c) < 1e-12:
        return np.minimum(X - a, b - X)
    R = H / (2.0 * abs(c))
    yc = 0.5 * H
    if c > 0:
        xl, xr = a - R, b + R
        dl = np.minimum(np.hypot(X - xl, Y - yc) - R, X - xl)
        dr = np.minimum(np.hypot(X - xr, Y - yc) - R, xr - X)
    else:
        dl = R - np.hypot(X - (a + R), Y - yc)
        dr = R - np.hypot(X - (b - R), Y - yc)
    return np.minimum(dl, dr)


def prepend_banner(run, eps=None, d0=None, mob=None, alph=None, rho_ice=RHO_ICE):
    """Put the startup-banner lines the tools grep for at the top of outp.txt
    (format copied from a real run's outp.txt). Call after fx.write_outp."""
    path = os.path.join(run, "outp.txt")
    body = open(path).read() if os.path.isfile(path) else ""
    lines = []
    if eps is not None:
        lines.append(f"   eps      =  {eps:.4e} m\n")
    lines.append(f"   rho_ice  =  {rho_ice:.4e} kg/m³,  rho_air = 1.3410e+00 kg/m³\n")
    if d0 is not None:
        lines.append(f"   d0_sub0          = {d0:.4e} m   (capillary length, K&P Eq. 13 at T0)\n")
    if mob is not None:
        lines.append(f"   mob_sub  =  {mob:.4e} m/s   [M&F: 4.33e-7]\n")
    if alph is not None:
        lines.append(f"   alph_sub =  {alph:.4e} 1/s\n")
    with open(path, "w") as fh:
        fh.write("".join(lines) + body)


def quiet(fn, *a, **k):
    """Call fn with stdout and stderr captured; return (result, stderr text)."""
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        out = fn(*a, **k)
    return out, err.getvalue()


def run_tool(script, *args):
    return subprocess.run([fx.PY, os.path.join(fx.POSTPROCESS, script), *args],
                          capture_output=True, text=True, cwd=fx.PROJECT)


def read_csv(path, comment="#"):
    with open(path) as fh:
        lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith(comment)]
    head = lines[0].split(",")
    rows = [dict(zip(head, ln.split(","))) for ln in lines[1:]]
    return head, rows


# Flat channel shared by the contact-angle and meniscus tests.
H, LX, EPS = 10e-6, 40e-6, 0.2e-6
FLAT_WALLS = [((0.0, 0.0), 0), ((H, 0.0), 1)]
FLAT_BOUNDS = ((0.0, 0.0), (H, 0.0))


def flat_mesh(eps_over_h=SQRT2, mapping=None):
    h = EPS / eps_over_h
    return fx.iga_rectangle(LX, H, int(round(LX / h)), int(round(H / h)), mapping=mapping)


def measure_angle(theta, eps_over_h=SQRT2, a=14e-6, b=26e-6, mapping=None,
                  transform=None, nu=601, nv=301):
    """Write one bridge snapshot and return contact_angle.measure() on it."""
    with tempfile.TemporaryDirectory() as tmp:
        run = fx.make_run(tmp)
        nrb = flat_mesh(eps_over_h, mapping)

        def f(X, Y):
            phi = fx.profile(bridge_distance(X, Y, a, b, H, theta), EPS)
            return (transform(phi) if transform else phi), T0, 1e-3
        fx.write_iga_run(run, nrb, [(0, f)])
        (_, X, Y, phi), = list(ca.read_snapshots(run, nu, nv))
        return ca.measure(X, Y, phi, EPS, FLAT_WALLS, FLAT_BOUNDS, 5.0)


# Wedge (symmetric about Ly/2, apex off the left edge), mesh parametric [0,1]^2
# like build_geometry_wedge.py / build_geometry_multi_grain.py produce.
W_LX, W_LY, W_WT, W_WM = 60e-6, 40e-6, 10e-6, 40e-6
W_M = (W_WM - W_WT) / (2.0 * W_LX)
W_EPS = 0.4e-6


def wedge_walls(mirror=False):
    """(bot_y0, bot_slope, top_y0, top_slope, apex_x)."""
    bot0, top0 = 0.5 * (W_LY - W_WT), 0.5 * (W_LY + W_WT)
    apex = -W_WT / (2.0 * W_M)
    if mirror:
        return bot0 - W_M * W_LX, +W_M, top0 + W_M * W_LX, -W_M, W_LX - apex
    return bot0, -W_M, top0, +W_M, apex


def wedge_mapping(walls, asym=0.0):
    b0, bs, t0, ts, _ = walls

    def m(U, V):
        X = U * W_LX
        yb, yt = b0 + bs * X, t0 + ts * X + asym * X
        return X, yb + V * (yt - yb)
    return m


def wedge_run(tmp, r1_of_t, r2_of_t, times, vL, vR, mirror=False, T_of_xy=None,
              sigma_override=None, geometry_extra=None, experiment_extra=None,
              banner=True, outp=True, asym=0.0, eps_over_h=SQRT2,
              extra_line=None, vanish_steps=(), beta_true=BETA0):
    """A wedge-band run whose vapor field satisfies the Gibbs-Thomson relation
    sigma*rho_vs(T)/rho_vs(T0) = d0*chi + beta*v_n exactly at both menisci,
    with chi = -1/r1 (inner, concave) and +1/r2 (outer, convex).

    rho_v is built as rho_vs(T) + rho_vs(T0)*(A + B*rho), rho = |X - apex|, so
    the (scaled) supersaturation is linear in radius and the tool's vapor-side
    linear extrapolation is exact."""
    walls = wedge_walls(mirror)
    b0, bs, t0_, ts, apex_x = walls
    apex_y = 0.5 * W_LY
    geom = {"-p": 2, "-C": 1, "-Lx": W_LX, "-Ly": W_LY, "-eps": W_EPS,
            "-wall_bot_y0": b0, "-wall_bot_slope": bs, "-wall_top_y0": t0_,
            "-wall_top_slope": ts, "-wedge_apex_x": apex_x, "-wedge_apex_y": apex_y,
            "-wall_faces": "y0,y1", **(geometry_extra or {})}
    exp = {"-temp": T0, "-beta_sub0": BETA0, **(experiment_extra or {})}
    run = fx.make_run(tmp, geom="wedge_2D", exp="wedge_exp", geometry=geom, experiment=exp)
    h = W_EPS / eps_over_h
    nrb = fx.iga_rectangle(1.0, 1.0, int(round(W_LX / h)), int(round(W_WM / h)),
                           mapping=wedge_mapping(walls, asym))
    steps = [10 * k for k in range(len(times))]
    snaps = []
    for k, (s, t) in enumerate(zip(steps, times)):
        r1, r2 = r1_of_t(t), r2_of_t(t)
        sL = D0 * (-1.0 / r1) + beta_true * vL
        sR = D0 * (+1.0 / r2) + beta_true * vR
        B = (sR - sL) / (r2 - r1)
        A = sL - B * r1
        gone = k in vanish_steps

        def f(X, Y, r1=r1, r2=r2, A=A, B=B, gone=gone):
            rho = np.hypot(X - apex_x, Y - apex_y)
            T = T0 if T_of_xy is None else T_of_xy(X, Y)
            phi = np.zeros_like(X) if gone else fx.profile(np.minimum(rho - r1, r2 - rho), W_EPS)
            return phi, T, rho_vs(T) + rho_vs(T0) * (A + B * rho)
        snaps.append((s, f))
    fx.write_iga_run(run, nrb, snaps)
    if outp:
        fx.write_outp(run, steps, times, extra=extra_line)
        if banner:
            prepend_banner(run, eps=W_EPS, d0=D0)
    return run


# ===========================================================================
# contact_angle.py
# ===========================================================================
class TestContactAngle(unittest.TestCase):

    def assertAngles(self, res, theta, tol, n=4):
        self.assertEqual(len(res), n, f"theta={theta}: {len(res)} estimates")
        for r in res:
            self.assertAlmostEqual(r["theta"], theta, delta=tol,
                                   msg=f"{r['meniscus']}@{r['wall']}: {r['theta']:.4f} vs {theta}")

    def test_theta_scan_flat_channel(self):
        """Circular-arc bridge, eps/h = 1.41 (production resolution). Measured
        error <= 0.013 deg over 30..150."""
        for theta in (30, 60, 90, 120, 150):
            with self.subTest(theta=theta):
                self.assertAngles(measure_angle(theta), theta, 0.05)

    def test_theta_near_0_and_180(self):
        """Nearly tangent arcs: 5 deg -> -0.07, 175 deg -> +0.10 deg error."""
        for theta in (5, 175):
            with self.subTest(theta=theta):
                self.assertAngles(measure_angle(theta), theta, 0.3)

    def test_under_resolved_eps(self):
        """eps/h = 0.7 and 0.4 (fewer than one control point per eps): error
        grows to 0.09 and 0.56 deg (175 deg: 1.4 deg) but stays bounded."""
        for eoh, tol in ((0.7, 0.15), (0.4, 0.8)):
            for theta in (30, 60, 90, 120, 150):
                with self.subTest(eps_over_h=eoh, theta=theta):
                    self.assertAngles(measure_angle(theta, eps_over_h=eoh), theta, tol)

    def test_under_resolved_near_zero_gives_no_estimate_not_a_wrong_one(self):
        """theta = 5 at eps/h = 0.4: the tool returns no estimate (silently, per
        snapshot), never a wrong number."""
        res = measure_angle(5, eps_over_h=0.4)
        for r in res:
            self.assertAlmostEqual(r["theta"], 5.0, delta=2.0)

    def test_contact_line_outside_domain_edge(self):
        """Left meniscus so close to x = 0 that its contact lines lie outside
        the domain (theta = 30, sagitta 2.9 um, vertex at 2 um). Rows near the
        walls lose the left crossing and are dropped; the fitted arc still
        extrapolates to the right angle."""
        res = measure_angle(30, a=2.0e-6)
        self.assertAngles(res, 30, 0.1)

    def test_phi_overshoot(self):
        """phi slightly above 1 and below 0 (a common solver artifact) does not
        move the angle."""
        res = measure_angle(60, transform=lambda p: 1.04 * p - 0.01)
        self.assertAngles(res, 60, 0.05)

    def test_graded_mesh_rows_not_uniform_in_y(self):
        """Mesh graded toward the walls (x ruled, y nonlinear in v). sample_phi's
        'exact' ruled inverse assumes y linear in v, so its orientation probe is
        approximate here; the angle must still be right."""
        def grade(GX, GY):
            s = GY / H
            return GX, H * (s - 0.08 * np.sin(2 * np.pi * s))
        for theta in (45, 135):
            with self.subTest(theta=theta):
                self.assertAngles(measure_angle(theta, mapping=grade), theta, 0.05)

    def test_non_ruled_mesh_fails_loudly(self):
        """x varying along a column breaks the row-scan contour; the tool must
        refuse rather than measure."""
        def shear(GX, GY):
            return GX + 0.05 * GY, GY
        with self.assertRaises(SystemExit):
            measure_angle(60, mapping=shear)

    def _wedge_angle(self, shape):
        walls = wedge_walls()
        b0, bs, t0_, ts, apex_x = walls
        apex_y = 0.5 * W_LY
        with tempfile.TemporaryDirectory() as tmp:
            run = fx.make_run(tmp)
            h = W_EPS / np.sqrt(2.0)
            nrb = fx.iga_rectangle(1.0, 1.0, int(round(W_LX / h)), int(round(W_WM / h)),
                                   mapping=wedge_mapping(walls))
            fx.write_iga_run(run, nrb, [(0, lambda X, Y: (shape(X, Y, apex_x, apex_y), T0, 1e-3))])
            (_, X, Y, phi), = list(ca.read_snapshots(run, 601, 301))
            wb, wt = (b0, bs), (t0_, ts)
            return ca.measure(X, Y, phi, W_EPS, [(wb, 0), (wt, 1)], (wb, wt), 5.0)

    def test_sloped_walls_apex_band_is_90(self):
        """Apex-centered annulus in a wedge: both menisci are curved (one
        concave, one convex) yet meet both sloped walls at exactly 90 deg."""
        def band(X, Y, ax, ay):
            rho = np.hypot(X - ax, Y - ay)
            return fx.profile(np.minimum(rho - 40e-6, 60e-6 - rho), W_EPS)
        self.assertAnglesWedge(self._wedge_angle(band), [90.0] * 4, 0.05)

    def test_sloped_walls_clipped_disc(self):
        """Disc on the wedge axis, clipped by both sloped walls: ice inside the
        circle, so cos(theta) = -d_perp/R with d_perp the perpendicular
        distance from the center to the wall line."""
        cx, R = 30e-6, 15e-6
        b0, bs, *_ = wedge_walls()
        d_perp = abs(0.5 * W_LY - (b0 + bs * cx)) / np.hypot(1.0, bs)
        theta = np.degrees(np.arccos(-d_perp / R))

        def disc(X, Y, ax, ay):
            return fx.disc(X, Y, cx, 0.5 * W_LY, R, W_EPS)
        self.assertAnglesWedge(self._wedge_angle(disc), [theta] * 4, 0.1)

    def assertAnglesWedge(self, res, expect, tol):
        self.assertEqual(len(res), len(expect))
        for r, e in zip(res, expect):
            self.assertAlmostEqual(r["theta"], e, delta=tol,
                                   msg=f"{r['meniscus']}@{r['wall']}: {r['theta']:.4f} vs {e:.4f}")

    def test_extrapolate_recovers_exponential(self):
        t = np.linspace(0.0, 5e4, 30)
        th = 60.0 + 25.0 * np.exp(-t / 1.2e4)
        ext = ca.extrapolate(t, th)
        self.assertIsNotNone(ext)
        self.assertAlmostEqual(ext[0], 60.0, delta=1e-4)
        self.assertAlmostEqual(ext[2], 1.2e4, delta=1.0)

    def test_extrapolate_without_times_returns_none(self):
        """Missing outp.txt gives NaN times; the extrapolation must decline."""
        self.assertIsNone(ca.extrapolate([np.nan] * 12, np.linspace(70, 60, 12)))

    def test_main_end_to_end_and_opts_precedence(self):
        """CLI on a run with flags in the files a real run uses: -eps, -Ly,
        -wall_faces in the geometry; surface energies in the experiment; one
        override after `--` (echoed in outp.txt) that must win."""
        theta = 60.0
        with tempfile.TemporaryDirectory() as tmp:
            run = fx.make_run(tmp, geom="chan", exp="ca_exp",
                              geometry={"-Lx": LX, "-Ly": H, "-eps": EPS, "-wall_faces": "y0,y1,x1"},
                              experiment={"-gamma_ia": 0.109, "-gamma_as": 0.3, "-gamma_is": 0.2455})
            nrb = flat_mesh()
            fx.write_iga_run(run, nrb, [(s, lambda X, Y: (fx.profile(
                bridge_distance(X, Y, 14e-6, 26e-6, H, theta), EPS), T0, 1e-3)) for s in (0, 5)])
            fx.write_outp(run, [0, 5], [0.0, 50.0], extra="-gamma_is 0.25")
            p = run_tool("contact_angle.py", "--dir", run)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(os.path.join(run, ca.CSV_NAME)) as fh:
                first = fh.readline()
            _, rows = read_csv(os.path.join(run, ca.CSV_NAME))
        young = np.degrees(np.arccos((0.3 - 0.25) / 0.109))
        self.assertIn("0.109,0.25,0.3", first)
        self.assertEqual(len(rows), 8)
        for r in rows:
            self.assertAlmostEqual(float(r["theta_young_deg"]), young, places=3)
            self.assertAlmostEqual(float(r["theta_deg"]), theta, delta=0.05)
        self.assertEqual({float(r["time_s"]) for r in rows}, {0.0, 50.0})

    def test_no_wall_faces_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = fx.make_run(tmp, geometry={"-Ly": H, "-eps": EPS, "-wall_faces": "x1"})
            fx.write_iga_run(run, flat_mesh(), [(0, lambda X, Y: (fx.lens(X, 1e-5, 3e-5, EPS), T0, 0))])
            p = run_tool("contact_angle.py", "--dir", run)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("-wall_faces", p.stderr)

    @unittest.skipUnless(fx.real_run_available(), "real run not present")
    def test_real_run_final_snapshot_matches_csv(self):
        """Re-measure the real run's last snapshot in memory (through symlinks;
        nothing is written into the run) and compare with its contact_angle.csv."""
        csv = os.path.join(fx.REAL_RUN, ca.CSV_NAME)
        if not os.path.isfile(csv):
            self.skipTest("no contact_angle.csv in the real run")
        sols = sorted(f for f in os.listdir(fx.REAL_RUN) if f.startswith("sol_"))
        last = sols[-1]
        step = int(last[4:9])
        with tempfile.TemporaryDirectory() as tmp:
            for f in ("igasol.dat", last):
                os.symlink(os.path.join(fx.REAL_RUN, f), os.path.join(tmp, f))
            (_, X, Y, phi), = list(ca.read_snapshots(tmp, 601, 301))
        from pplib import opt_float, read_opts
        opts = read_opts(fx.REAL_RUN)
        eps, Ly = opt_float(opts, "-eps"), opt_float(opts, "-Ly")
        res = ca.measure(X, Y, phi, eps, [((0.0, 0.0), 0), ((Ly, 0.0), 1)],
                         ((0.0, 0.0), (Ly, 0.0)), 5.0)
        _, rows = read_csv(csv)
        ref = {(r["meniscus"], r["wall"]): float(r["theta_deg"]) for r in rows
               if int(r["step"]) == step}
        self.assertEqual(len(res), len(ref))
        for r in res:
            self.assertAlmostEqual(r["theta"], ref[(r["meniscus"], r["wall"])], delta=1e-3)


# ===========================================================================
# meniscus_velocity.py
# ===========================================================================
def lens_run(tmp, a_of_t, b_of_t, times, steps=None, eps_over_h=SQRT2, outp=True,
             phi_of=None, mapping=None, name="mv_exp", layout_root=None):
    steps = list(range(0, 10 * len(times), 10)) if steps is None else steps
    run = fx.make_run(layout_root or tmp, geom="chan", exp=name,
                      geometry={"-Lx": LX, "-Ly": H, "-eps": EPS, "-wall_faces": "y0,y1"})
    nrb = flat_mesh(eps_over_h, mapping)
    snaps = []
    for s, t in zip(steps, times):
        a, b = a_of_t(t), b_of_t(t)
        if phi_of is None:
            snaps.append((s, lambda X, Y, a=a, b=b: (fx.lens(X, a, b, EPS), T0, 1e-3)))
        else:
            snaps.append((s, lambda X, Y, a=a, b=b, t=t: (phi_of(X, Y, a, b, t), T0, 1e-3)))
    fx.write_iga_run(run, nrb, snaps)
    if outp:
        fx.write_outp(run, steps, times)
    return run


def mv_rows(*dirs, nu="801", nv="101"):
    args = []
    for d in dirs:
        args += ["--dir", d]
    p = run_tool("meniscus_velocity.py", *args, "--nu", nu, "--nv", nv)
    if p.returncode != 0:
        return p, None
    _, rows = read_csv(os.path.join(dirs[-1], mv.CSV_NAME))
    return p, {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


class TestMeniscusVelocity(unittest.TestCase):
    VA, VB = 3e-9, 1e-9      # left face moves left, right face moves right: both grow

    def test_flat_lens_uneven_times(self):
        """Plug a(t) = 15um - VA t, b(t) = 25um + VB t, sampled at uneven times.
        The tool reports the MEAN of the two face speeds, (VA+VB)/2, positive
        for growth; positions match a(t), b(t)."""
        t = np.array([0.0, 100.0, 250.0, 300.0, 700.0, 1000.0])
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6 - self.VA * t, lambda t: 25e-6 + self.VB * t, t)
            p, d = mv_rows(run)
        self.assertEqual(p.returncode, 0, p.stderr)
        v = 0.5 * (self.VA + self.VB)
        for k in ("v_wall_m_s", "v_mid_m_s", "v_mean_m_s"):
            np.testing.assert_allclose(d[k], v, rtol=1e-4, err_msg=k)
        np.testing.assert_allclose(d["x_mid_l"], 15e-6 - self.VA * t, atol=2e-11)
        np.testing.assert_allclose(d["x_mid_r"], 25e-6 + self.VB * t, atol=2e-11)
        np.testing.assert_allclose(d["area_m2"], H * (10e-6 + (self.VA + self.VB) * t), rtol=1e-5)

    def test_shrinking_lens_is_negative(self):
        t = np.linspace(0.0, 1000.0, 5)
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6 + 2e-9 * t, lambda t: 25e-6 - 2e-9 * t, t)
            p, d = mv_rows(run)
        np.testing.assert_allclose(d["v_mid_m_s"], -2e-9, rtol=1e-4)
        np.testing.assert_allclose(d["v_mean_m_s"], -2e-9, rtol=1e-4)

    def test_stationary_interface_is_zero(self):
        t = np.array([0.0, 10.0, 11.0, 500.0])
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15.3e-6, lambda t: 24.9e-6, t)
            p, d = mv_rows(run)
        for k in ("v_wall_m_s", "v_mid_m_s", "v_mean_m_s"):
            self.assertLess(np.abs(d[k]).max(), 1e-18, k)

    def test_quadratic_motion_uneven_times(self):
        """a(t) = a0 - c t^2: np.gradient on uneven times is exact for a
        quadratic at interior points (one-sided at the ends)."""
        c = 1e-12
        t = np.array([0.0, 100.0, 180.0, 400.0, 450.0, 800.0, 1000.0])
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6 - c * t * t, lambda t: 25e-6 + c * t * t, t)
            p, d = mv_rows(run)
        # 0.2 % observed: position ripple (~7e-11 m = 3e-4 h) over short dt
        np.testing.assert_allclose(d["v_mid_m_s"][1:-1], 2 * c * t[1:-1], rtol=5e-3)

    def test_under_resolved_eps(self):
        t = np.linspace(0.0, 1000.0, 5)
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6 - self.VA * t, lambda t: 25e-6 + self.VB * t, t,
                           eps_over_h=0.4)
            p, d = mv_rows(run)
        np.testing.assert_allclose(d["v_mid_m_s"], 0.5 * (self.VA + self.VB), rtol=2e-3)

    def test_phi_overshoot(self):
        t = np.linspace(0.0, 1000.0, 5)
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6 - self.VA * t, lambda t: 25e-6 + self.VB * t, t,
                           phi_of=lambda X, Y, a, b, t: 1.03 * fx.lens(X, a, b, EPS) - 0.01)
            p, d = mv_rows(run)
        np.testing.assert_allclose(d["v_mid_m_s"], 0.5 * (self.VA + self.VB), rtol=1e-3)

    def test_interface_near_domain_edge(self):
        """Left face 2 eps from x = 0: half of its diffuse band is cut off."""
        t = np.linspace(0.0, 100.0, 5)
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 2 * EPS - 1e-9 * t, lambda t: 25e-6 + 1e-9 * t, t)
            p, d = mv_rows(run)
        self.assertEqual(p.returncode, 0, p.stderr)
        np.testing.assert_allclose(d["v_mid_m_s"], 1e-9, rtol=1e-2)

    def test_curved_bridge_mean_speed(self):
        """Concave (theta = 60) bridge translating rigidly to the right at U:
        the left face recedes, the right advances, so v_mid = 0, v_wall = 0,
        dA/dt = 0 and v_mean = 0."""
        U = 1e-9
        t = np.linspace(0.0, 1000.0, 5)
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 14e-6 + U * t, lambda t: 26e-6 + U * t, t,
                           phi_of=lambda X, Y, a, b, t: fx.profile(
                               bridge_distance(X, Y, a, b, H, 60.0), EPS))
            p, d = mv_rows(run)
        for k in ("v_wall_m_s", "v_mid_m_s", "v_mean_m_s"):
            self.assertLess(np.abs(d[k]).max(), 1e-3 * U, k)
        # refine_tanh fits a straight line to atanh(2 phi - 1) across a CURVED
        # band: a constant 2.5e-10 m (1.3e-3 eps) offset, harmless for speeds.
        np.testing.assert_allclose(d["x_mid_l"], 14e-6 + U * t, atol=5e-10)

    def test_two_snapshots_fail_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6, lambda t: 25e-6, [0.0, 10.0])
            p, _ = mv_rows(run)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("at least 3", p.stderr)

    def test_missing_outp_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6, lambda t: 25e-6, [0.0, 10.0, 20.0], outp=False)
            p, _ = mv_rows(run)
        self.assertNotEqual(p.returncode, 0)

    def test_lens_vanishes_in_one_snapshot(self):
        """A snapshot with no ice is skipped (silently); the speed across the
        gap is still right because positions are differentiated against t."""
        t = np.linspace(0.0, 1000.0, 6)

        def phi_of(X, Y, a, b, tt):
            return np.zeros_like(X) if tt == t[2] else fx.lens(X, a, b, EPS)
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 15e-6 - 1e-9 * t, lambda t: 25e-6 + 1e-9 * t, t,
                           phi_of=phi_of)
            p, d = mv_rows(run)
        self.assertEqual(len(d["time_s"]), 5)
        np.testing.assert_allclose(d["v_mid_m_s"], 1e-9, rtol=1e-4)

    def test_restart_chain_without_duplicate(self):
        """Two legs, each with its own outp.txt, step numbers restarting at 1."""
        VA = VB = 1e-9
        tA = np.array([0.0, 200.0, 400.0])
        tB = np.array([600.0, 800.0, 1000.0])
        with tempfile.TemporaryDirectory() as tmp:
            A = lens_run(tmp, lambda t: 15e-6 - VA * t, lambda t: 25e-6 + VB * t, tA,
                         steps=[0, 1, 2], name="legA")
            B = lens_run(tmp, lambda t: 15e-6 - VA * t, lambda t: 25e-6 + VB * t, tB,
                         steps=[1, 2, 3], name="legB")
            p, d = mv_rows(A, B)
        self.assertEqual(p.returncode, 0, p.stderr)
        np.testing.assert_allclose(d["time_s"], np.r_[tA, tB])
        np.testing.assert_allclose(d["v_mid_m_s"], 1e-9, rtol=1e-4)

    # BUG: meniscus_velocity.py:147-171 -- a restart leg writes sol_00000.dat at
    # t = t_start (src/monitoring.c:395, step <= 1 always written; step resets to
    # 0 because the solver only calls TSSetTime), which duplicates the previous
    # leg's last snapshot. The tool keeps both rows; np.gradient then divides by
    # dt = 0 and the CSV carries NaN/inf velocities at the join (verified by
    # hand: v = nan at the two duplicated rows). Right: drop or merge rows with
    # equal t (pplib.iter_snapshots keeps one per step; dedupe on t here).
    @unittest.expectedFailure
    def test_restart_chain_with_duplicate_join(self):
        VA = VB = 1e-9
        tA = np.array([0.0, 200.0, 400.0])
        tB = np.array([400.0, 600.0, 800.0])
        with tempfile.TemporaryDirectory() as tmp:
            A = lens_run(tmp, lambda t: 15e-6 - VA * t, lambda t: 25e-6 + VB * t, tA,
                         steps=[0, 1, 2], name="legA")
            B = lens_run(tmp, lambda t: 15e-6 - VA * t, lambda t: 25e-6 + VB * t, tB,
                         steps=[0, 1, 2], name="legB")
            p, d = mv_rows(A, B)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertTrue(np.isfinite(d["v_mid_m_s"]).all(), d["v_mid_m_s"])
        np.testing.assert_allclose(d["v_mid_m_s"], 1e-9, rtol=1e-4)

    # BUG: meniscus_velocity.py:99-106 and 118-125 -- crossings_in() keeps only the
    # first and last phi = 0.5 crossing of each row, so a bridge that has split
    # into two plugs is read as one: the two inner faces drop out of the arc
    # length L and v_mean = (dA/dt)/L comes out 2x too large, with no warning.
    # Right: count every crossing (L = 4H here) or refuse rows with != 2 crossings.
    @unittest.expectedFailure
    def test_split_lens_mean_speed(self):
        v = 1e-9
        t = np.linspace(0.0, 1000.0, 5)

        def two_plugs(X, Y, a, b, tt):
            return fx.union(fx.lens(X, a, 16e-6 + v * tt, EPS),
                            fx.lens(X, 24e-6 - v * tt, b, EPS))
        with tempfile.TemporaryDirectory() as tmp:
            run = lens_run(tmp, lambda t: 10e-6 - v * t, lambda t: 30e-6 + v * t, t,
                           phi_of=two_plugs)
            p, d = mv_rows(run)
        np.testing.assert_allclose(d["v_mean_m_s"], v, rtol=1e-3)

    # BUG: meniscus_velocity.py:69 and 116 -- the snapshot keeps y1d = C[0, :, 1],
    # the node heights of the FIRST column only, and integrates the area with it
    # as if every column had that y. On a mapped (wedge) mesh the channel width
    # varies with x, so area and dA/dt (hence v_mean) are silently wrong (here
    # the area is 0.4x the truth: column 0 is the 10 um throat). Right:
    # integrate with the 2-D physical coordinates (pplib.node_weights).
    @unittest.expectedFailure
    def test_mapped_mesh_area(self):
        walls = wedge_walls()
        b0, bs, t0_, ts, _ = walls
        a, b = 20e-6, 40e-6
        with tempfile.TemporaryDirectory() as tmp:
            run = fx.make_run(tmp, geometry={"-Ly": W_LY, "-eps": W_EPS})
            h = W_EPS / np.sqrt(2.0)
            nrb = fx.iga_rectangle(1.0, 1.0, int(round(W_LX / h)), int(round(W_WM / h)),
                                   mapping=wedge_mapping(walls))
            fx.write_iga_run(run, nrb, [(0, lambda X, Y: (fx.lens(X, a, b, W_EPS), T0, 0))])
            (_, x1d, y1d, phi), = list(mv.snapshots(run, 1201, 201))
        g = mv.measure(x1d, y1d, phi, W_EPS, W_LY)
        width = lambda x: (t0_ + ts * x) - (b0 + bs * x)   # noqa: E731
        exact = 0.5 * (width(a) + width(b)) * (b - a)       # width is linear in x
        self.assertAlmostEqual(g["area"] / exact, 1.0, delta=2e-3)


# ===========================================================================
# wedge_gt_velocity.py
# ===========================================================================
class TestWedgeGTVelocity(unittest.TestCase):
    VL, VR = 2e-11, 5e-12
    TIMES = np.linspace(0.0, 2e5, 9)

    def analyse(self, run, **kw):
        kw.setdefault("nu", 2001)
        kw.setdefault("nv", 61)
        kw.setdefault("vn_win", 5)
        (data, meta), err = quiet(wg.analyse, run, progress=False, **kw)
        return data, meta, err

    def band(self, tmp, **kw):
        return wedge_run(tmp, lambda t: 40e-6 - self.VL * t, lambda t: 60e-6 + self.VR * t,
                         self.TIMES, self.VL, self.VR, **kw)

    def check(self, data, vL=None, vR=None, rtol=2e-3):
        vL = self.VL if vL is None else vL
        vR = self.VR if vR is None else vR
        np.testing.assert_allclose(data["vn_left_meas"], vL, rtol=rtol, err_msg="left meas")
        np.testing.assert_allclose(data["vn_right_meas"], vR, rtol=rtol, err_msg="right meas")
        np.testing.assert_allclose(data["vn_left_pred"], vL, rtol=rtol, err_msg="left pred")
        np.testing.assert_allclose(data["vn_right_pred"], vR, rtol=rtol, err_msg="right pred")

    def test_known_answer_isothermal(self):
        """Measured and predicted v_n both recover the imposed speeds (inner
        concave front 2e-11, outer convex 5e-12 m/s); chi has the right sign
        and the circle-fit curvature agrees with 1/r."""
        with tempfile.TemporaryDirectory() as tmp:
            data, meta, _ = self.analyse(self.band(tmp))
        self.assertTrue(meta["left_is_inner"])
        self.assertAlmostEqual(meta["d0"], D0, delta=1e-13)
        self.check(data, rtol=2e-4)
        r1 = 40e-6 - self.VL * self.TIMES
        r2 = 60e-6 + self.VR * self.TIMES
        np.testing.assert_allclose(data["r_left"], r1, rtol=1e-5)
        np.testing.assert_allclose(data["r_right"], r2, rtol=1e-5)
        np.testing.assert_allclose(data["chi_left"], -1.0 / r1, rtol=1e-5)
        np.testing.assert_allclose(data["chi_left_fit"], -1.0 / r1, rtol=1e-2)
        np.testing.assert_allclose(data["chi_right_fit"], 1.0 / r2, rtol=1e-2)

    def test_mirrored_apex_right_of_domain(self):
        """Apex off the RIGHT edge: the inner (concave) front is now the right
        one; signs must follow."""
        with tempfile.TemporaryDirectory() as tmp:
            data, meta, _ = self.analyse(self.band(tmp, mirror=True))
        self.assertFalse(meta["left_is_inner"])
        # left meniscus is now the OUTER one (radius r2), right is the inner (r1)
        self.check(data, vL=self.VR, vR=self.VL, rtol=2e-4)
        self.assertTrue((data["chi_right"] < 0).all() and (data["chi_left"] > 0).all())

    def test_nonisothermal_rho_vs_ratio(self):
        """Field at -21 C while -temp is -20 C: the prediction must carry
        rho_vs(T)/rho_vs(T0) = 0.91 (dropping it would move v_pred by ~5 %).

        A uniform offset is used rather than a gradient: with a gradient the
        fixture's O(h^2) spline error in rho_vs(T(x)) (~2e-7 relative) is
        already a 2 % error in sigma ~ 1e-5, i.e. the test would measure the
        fixture, not the tool. The same sensitivity applies to real data:
        sigma is a 1e-5 difference of two 1e-3 numbers."""
        with tempfile.TemporaryDirectory() as tmp:
            run = self.band(tmp, T_of_xy=lambda X, Y: -21.0 + 0.0 * X)
            data, _, _ = self.analyse(run)
        self.check(data, rtol=2e-4)

    def test_unevenly_spaced_times(self):
        t = np.array([0.0, 1e4, 1.5e4, 5e4, 6e4, 1.2e5, 1.3e5, 2e5])
        with tempfile.TemporaryDirectory() as tmp:
            run = wedge_run(tmp, lambda t: 40e-6 - self.VL * t, lambda t: 60e-6 + self.VR * t,
                            t, self.VL, self.VR)
            data, _, _ = self.analyse(run)
        self.check(data, rtol=1e-3)

    def test_stationary_band(self):
        """sigma = d0*chi exactly and nothing moves: both v_n ~ 0, far below
        the d0*chi/beta scale (~3e-14 m/s)."""
        with tempfile.TemporaryDirectory() as tmp:
            run = wedge_run(tmp, lambda t: 40e-6, lambda t: 60e-6, self.TIMES, 0.0, 0.0)
            data, _, _ = self.analyse(run)
        scale = D0 / 40e-6 / BETA0
        for k in ("vn_left_meas", "vn_right_meas", "vn_left_pred", "vn_right_pred"):
            self.assertLess(np.nanmax(np.abs(data[k])), 1e-3 * scale, k)

    def test_under_resolved_eps(self):
        with tempfile.TemporaryDirectory() as tmp:
            data, _, _ = self.analyse(self.band(tmp, eps_over_h=0.5))
        self.check(data, rtol=5e-3)

    def test_extra_opts_override_beta(self):
        """-beta_sub0 passed after `--` (echoed in outp.txt) beats the
        experiment file; the vapor field is built with the overriding value."""
        beta = 2.0 * BETA0
        with tempfile.TemporaryDirectory() as tmp:
            run = self.band(tmp, extra_line=f"-beta_sub0 {beta:.6e}", beta_true=beta)
            data, meta, _ = self.analyse(run)
        self.assertAlmostEqual(meta["beta"], beta, delta=1.0)
        self.check(data, rtol=2e-4)

    def test_band_vanishes_in_one_snapshot_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            data, _, err = self.analyse(self.band(tmp, vanish_steps=(4,)))
        self.assertIn("0 centerline crossings", err)
        self.assertEqual(len(data["time"]), len(self.TIMES) - 1)
        self.check(data, rtol=2e-4)

    def test_vts_source_curvilinear_rows(self):
        """--source vtkOut on a wedge grid whose rows are NOT at constant y; the
        centerline is recovered column by column (centreline_weights)."""
        walls = wedge_walls()
        b0, bs, t0_, ts, apex_x = walls
        with tempfile.TemporaryDirectory() as tmp:
            run = self.band(tmp)
            U, V = np.meshgrid(np.linspace(0, 1, 1201), np.linspace(0, 1, 120))
            X, Y = wedge_mapping(walls)(U, V)
            apex_y = 0.5 * W_LY
            snaps = []
            for k, t in enumerate(self.TIMES):
                r1, r2 = 40e-6 - self.VL * t, 60e-6 + self.VR * t
                sL, sR = -D0 / r1 + BETA0 * self.VL, D0 / r2 + BETA0 * self.VR
                B = (sR - sL) / (r2 - r1)
                rho = np.hypot(X - apex_x, Y - apex_y)
                snaps.append((10 * k, t, fx.profile(np.minimum(rho - r1, r2 - rho), W_EPS),
                              T0 + 0 * X, rho_vs(T0) * (1 + sL + B * (rho - r1))))
            fx.write_vts_run(run, X, Y, snaps)
            self.assertFalse(np.allclose(Y[60, :], Y[60, 0]))
            data, _, _ = self.analyse(run, source="vtkOut")
        self.check(data, rtol=1e-2)

    def test_asymmetric_walls_fail_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.band(tmp, asym=0.05)
            with self.assertRaises(SystemExit):
                self.analyse(run)

    def test_mob_sub_override_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.band(tmp, experiment_extra={"-mob_sub": 1e-8})
            with self.assertRaises(SystemExit):
                self.analyse(run)

    def test_missing_d0_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.band(tmp, banner=False)
            with self.assertRaises(SystemExit) as cm:
                self.analyse(run)
        self.assertIn("-d0_sub0", str(cm.exception))

    def test_missing_outp_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.band(tmp, outp=False, geometry_extra={"-d0_sub0": D0})
            with self.assertRaises(SystemExit):
                self.analyse(run)

    def test_two_snapshots_fail_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = wedge_run(tmp, lambda t: 40e-6, lambda t: 60e-6, [0.0, 1e4], 0.0, 0.0)
            with self.assertRaises(SystemExit):
                self.analyse(run)

    def test_rectangle_mesh_with_physical_knots_fails(self):
        """_iter_igasol samples u, v in [0, 1] (wedge_gt_velocity.py:199-200);
        a solver-built rectangle has knots in metres. The tool must not return
        numbers; it currently dies with a bare AssertionError from igakit."""
        with tempfile.TemporaryDirectory() as tmp:
            run = fx.make_run(tmp)
            fx.write_iga_run(run, fx.iga_rectangle(40e-6, 20e-6, 40, 20),
                             [(0, lambda X, Y: (fx.lens(X, 1e-5, 3e-5, 1e-6), T0, 1e-3))])
            with self.assertRaises((SystemExit, AssertionError, ValueError)):
                list(wg._iter_igasol(run, 10e-6, 101, 11))


# ===========================================================================
# gt_balance.py (needs ParaView's pvpython)
# ===========================================================================
def gt_run(tmp, T, v0, eps=0.4e-6, L=30e-6, R=10e-6, n=151, alph=1.6087e7):
    """Disc (ice inside, chi = +1/R) at t = 0 and a vapor hole in ice
    (chi = -1/R) at t = 100 s, uniform T. The vapor density is set so that
    s = d0*chi + beta*v0 on the interface, with d0 = d0_sub0 and
    beta = 5 rho_ice/(eps alph rho_vs(T0)) -- the relations the solver's
    mob_sub/alph_sub derivation implies (cryosuction_main.c:1097-1166)."""
    mob = alph * D0 * rho_vs(T) / (15.0 * RHO_ICE)
    beta = 5.0 * RHO_ICE / (eps * alph * rho_vs(T))
    run = fx.make_run(tmp, geometry={"-Lx": L, "-Ly": L, "-eps": eps},
                      experiment={"-temp": T})
    X, Y = fx.grid(L, L, n, n)
    disc = fx.disc(X, Y, L / 2, L / 2, R, eps)
    rv1 = rho_vs(T) * (1 + D0 / R + beta * v0) + 0 * X
    rv2 = rho_vs(T) * (1 - D0 / R + beta * v0) + 0 * X
    fx.write_vts_run(run, X, Y, [(0, 0.0, disc, T + 0 * X, rv1),
                                 (10, 100.0, 1.0 - disc, T + 0 * X, rv2)])
    fx.write_outp(run, [0, 10], [0.0, 100.0])
    prepend_banner(run, eps=eps, d0=D0, mob=mob, alph=alph)
    return run, beta


def gt_csv(run):
    p = subprocess.run([PVPYTHON, os.path.join(fx.POSTPROCESS, "gt_balance.py"), run,
                        "--times", "2"], capture_output=True, text=True, cwd=fx.PROJECT)
    if p.returncode != 0:
        raise AssertionError(p.stdout[-1500:] + p.stderr[-1500:])
    a = np.loadtxt(os.path.join(run, "gt_balance.csv"), delimiter=",", skiprows=1)
    return {k: a[:, i] for i, k in enumerate(
        ("t", "x", "y", "chi", "d0_chi", "s", "resid", "v_n"))}


@unittest.skipUnless(os.access(PVPYTHON, os.X_OK), "pvpython not available")
class TestGTBalance(unittest.TestCase):
    R, V0 = 10e-6, 2e-10

    def test_disc_and_hole_known_answer(self):
        """chi = +1/R (disc) and -1/R (hole) to 1.5 %; v_n = v0 to 2 %."""
        with tempfile.TemporaryDirectory() as tmp:
            run, _ = gt_run(tmp, T0, self.V0)
            d = gt_csv(run)
        for t, sign in ((0.0, 1.0), (100.0, -1.0)):
            sel = d["t"] == t
            self.assertGreater(sel.sum(), 100)
            self.assertAlmostEqual(np.median(d["chi"][sel]) * self.R, sign, delta=0.015)
            self.assertAlmostEqual(np.median(d["v_n"][sel]) / self.V0, 1.0, delta=0.02)
        np.testing.assert_allclose(d["d0_chi"] / d["chi"], D0, rtol=1e-3)

    # BUG: gt_balance.py:226 -- rho_vs is evaluated at a hard-coded T0_K = 253.15
    # (-20 C) instead of the run's -temp, so d0 = 15 M rho_ice/(alph rho_vs) and
    # beta are off by rho_vs(T)/rho_vs(-20 C) for every other temperature. At
    # -40 C this gives d0 = 0.124 d0_sub0 and v_n = 8.6x the truth; the tool
    # prints "ratio = 0.124311" next to "<- must agree" but still writes the CSV.
    # Verified by hand with this fixture. Right: take T0 from read_opts -temp (and
    # rho_air from -rho_air) and use pplib.rho_vs.
    @unittest.expectedFailure
    def test_temperature_other_than_minus_20(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, _ = gt_run(tmp, -40.0, self.V0)
            d = gt_csv(run)
        np.testing.assert_allclose(d["d0_chi"] / d["chi"], D0, rtol=1e-3)
        sel = d["t"] == 0.0
        self.assertAlmostEqual(np.median(d["v_n"][sel]) / self.V0, 1.0, delta=0.02)


class TestGTBalanceLog(unittest.TestCase):
    """scan_log needs no ParaView."""

    def test_missing_outp_fails_loudly(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = fx.make_run(tmp)
            with self.assertRaises(SystemExit):
                gtb.scan_log(run)

    def test_scan_log_synthetic_banner(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = fx.make_run(tmp)
            fx.write_outp(run, [0], [0.0])
            prepend_banner(run, eps=3.9358e-7, d0=D0, mob=1.0054e-9, alph=1.6087e7)
            c = gtb.scan_log(run)
        self.assertEqual((c["eps"], c["mob_sub"], c["alph_sub"], c["rho_ice"], c["d0_sub0"]),
                         (3.9358e-7, 1.0054e-9, 1.6087e7, 919.0, D0))

    @unittest.skipUnless(fx.real_run_available(), "real run not present")
    def test_scan_log_real_run(self):
        """Read-only: the real banner yields the values printed in it."""
        c = gtb.scan_log(fx.REAL_RUN)
        self.assertEqual(c["eps"], 3.9358e-7)
        self.assertEqual(c["rho_ice"], 919.0)
        self.assertEqual(c["mob_sub"], 1.0054e-9)
        self.assertEqual(c["alph_sub"], 1.6087e7)
        self.assertEqual(c["d0_sub0"], 1.0152e-9)
        d0 = 15 * c["mob_sub"] * c["rho_ice"] / (c["alph_sub"] * float(rho_vs(T0)))
        self.assertAlmostEqual(d0 / c["d0_sub0"], 1.0, delta=2e-4)


if __name__ == "__main__":
    unittest.main()
