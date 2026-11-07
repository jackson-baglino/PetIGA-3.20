"""Tests for the cryosuction study scripts (studies/cryosuction/scripts/).

Every number a script prints is checked against an independent calculation
written here from the physics: the saturation curve is transcribed from the C
source (src/material_properties.c, RhoVS_I) and checked to still be there, not
taken from preprocess/mesh_params; derivatives are analytic, not central
differences; the harvest cutoff is solved exactly by bisection rather than
from the linearized formula; areas are clipped polygons, not the script's
quadrature. Scripts run as subprocesses with venv_cryosuction/bin/python.

Nothing here runs the solver. preflight.sh is run only as a COPY inside a
temporary tree whose "binary" is a stub shell script that prints a fake
header; the real binary is never started.

Run from the project root:
    venv_cryosuction/bin/python -m unittest discover -s studies/cryosuction/tests -t studies/cryosuction/tests -v

Tests marked expectedFailure document a bug; each carries a "# BUG:" comment
naming file:line and the wrong vs right behavior.
"""
from __future__ import annotations

import math
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PY = str(ROOT / "venv_cryosuction" / "bin" / "python")
SCRIPTS = ROOT / "studies" / "cryosuction" / "scripts"
sys.path.insert(0, str(ROOT))
from postprocess.tests import fixtures as F  # noqa: E402

# ---------------------------------------------------------------------------
# Independent physics (from the C source and the papers, not from mesh_params)
# ---------------------------------------------------------------------------
RHO_ICE = 919.0          # src/cryosuction_main.c user.rho_ice
RHO_AIR = 1.341          # src/cryosuction_main.c user.rho_air
P_ATM = 101325.0         # src/material_properties.c Patm
BB = 0.62                # src/material_properties.c bb
DV0 = 2.178e-5           # src/cryosuction_main.c user.dif_vap
M_HK, KB_HK = 2.99e-26, 1.38e-23           # Hertz-Knudsen constants used by the planning tools
M_D0, KB_D0 = 2.9915e-26, 1.380649e-23     # include/NASA_types.h (d0)
GAMMA = 0.109            # src/cryosuction_main.c gamma_ia_bulk


def _ashrae_from_c_source():
    """K0..K5 read from src/material_properties.c, so a change there is caught."""
    txt = (ROOT / "src" / "material_properties.c").read_text()
    out = []
    for k in range(6):
        m = re.search(rf"\bK{k}\s*=\s*([-+0-9.eE]+)", txt)
        out.append(float(m.group(1)))
    return out


K_ASH = _ashrae_from_c_source()


def p_vs(T_K):
    k = K_ASH
    return math.exp(k[0] / T_K + k[1] + k[2] * T_K + k[3] * T_K**2 + k[4] * T_K**3 + k[5] * math.log(T_K))


def rho_vs(T_C):
    """The solver's saturation density, RhoVS_I: rho_air * 0.62 * p / (Patm - p)."""
    p = p_vs(T_C + 273.15)
    return RHO_AIR * BB * p / (P_ATM - p)


def s_analytic(T_C):
    """d ln(rho_vs)/dT for RhoVS_I, analytically: dlnp/dT * Patm / (Patm - p)."""
    t = T_C + 273.15
    k = K_ASH
    dlnp = -k[0] / t**2 + k[2] + 2 * k[3] * t + 3 * k[4] * t**2 + k[5] / t
    p = p_vs(t)
    return dlnp * P_ATM / (P_ATM - p)


def d_v(T_C):
    return DV0 * ((T_C + 273.15) / 273.15) ** 1.81


def beta0(T_C, alpha):
    """Unscaled K&P beta0 = (1/alpha) sqrt(2 pi m / (kB T)) rho_ice / rho_vs [s/m]."""
    return (1 / alpha) * math.sqrt(2 * math.pi * M_HK / (KB_HK * (T_C + 273.15))) * RHO_ICE / rho_vs(T_C)


def d0(T_C):
    """K&P Eq. 13, d0 = gamma (m / rho_ice) / (kB T)."""
    return GAMMA * (M_D0 / RHO_ICE) / (KB_D0 * (T_C + 273.15))


def k_diff(T_C):
    """Diffusion resistance per length, K = rho_ice / (D_v rho_vs) [s/m^2]."""
    return RHO_ICE / (d_v(T_C) * rho_vs(T_C))


def p_mk(T):
    """Murphy and Koop 2005, Eq. 7 (ice)."""
    return math.exp(9.550426 - 5723.265 / T + 3.53068 * math.log(T) - 0.00728332 * T)


def dln_rho_mk(T):
    """d/dT ln(p_mk(T) M / (R T)), analytic."""
    return 5723.265 / T**2 + 3.53068 / T - 0.00728332 - 1.0 / T


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def run(script, *args, env=None, cwd=ROOT):
    return subprocess.run([PY, str(SCRIPTS / script), *map(str, args)], capture_output=True, text=True,
                          cwd=cwd, env=env, timeout=60)


NUM = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|nan|inf)"


def num(pattern, text, group=1):
    m = re.search(pattern.replace("{N}", NUM), text)
    if m is None:
        raise AssertionError(f"pattern {pattern!r} not found in:\n{text}")
    return float(m.group(group))


def nums(pattern, text):
    return [float(x) for x in re.findall(pattern.replace("{N}", NUM), text)]


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


class AnchorValues(unittest.TestCase):
    """The independent physics above reproduces the values quoted in the project."""

    def test_anchor_values_at_minus_20(self):
        self.assertLess(rel(beta0(-20, 1e-2), 7.94e5), 1e-3)      # K&P beta0 at alpha_c 1e-2
        self.assertLess(rel(d0(-20), 1.015e-9), 1e-3)
        self.assertLess(rel(d_v(-20), 1.898e-5), 1e-3)
        self.assertLess(rel(rho_vs(-20), 8.4869e-4), 1e-4)

    def test_ashrae_pressure_matches_murphy_koop(self):
        # The solver's pressure fit agrees with M&K Eq. 7 to 0.1 % from 173 K up
        # (the density differs because rho_air is held at 1.341 kg/m^3).
        for T in (173.15, 223.15, 253.15, 273.15):
            self.assertLess(rel(p_vs(T), p_mk(T)), 2e-3, T)


# ---------------------------------------------------------------------------
# front_speed.py
# ---------------------------------------------------------------------------
class FrontSpeed(unittest.TestCase):
    def speed(self, *args):
        r = run("front_speed.py", *args)
        self.assertEqual(r.returncode, 0, r.stderr)
        return (num(r"beta0 = {N} s/m", r.stdout), num(r"K = {N} s/m\^2", r.stdout),
                num(r"v_n = {N} m/s", r.stdout), r.stdout)

    def test_series_law_against_independent(self):
        for T, a, h, l in ((-20, 1e-2, 0.99, 9.4e-5), (-5, 1e-3, 0.9, 1e-3), (-60, 1.0, 0.5, 1e-6)):
            with self.subTest(T=T, a=a, h=h, l=l):
                b, K, v, _ = self.speed("--T0", T, "--alpha", a, "--h", h, "--l", l)
                self.assertLess(rel(b, beta0(T, a)), 1e-4)
                self.assertLess(rel(K, k_diff(T)), 1e-4)
                self.assertLess(rel(v, abs(h - 1) / (beta0(T, a) + k_diff(T) * l)), 1e-4)

    def test_l_zero_gives_kinetic_bound(self):
        _, _, v, out = self.speed("--alpha", 1e-2, "--h", 0.99, "--l", 0)
        self.assertLess(rel(v, 0.01 / beta0(-20, 1e-2)), 1e-4)
        self.assertLess(rel(num(r"\|sigma\|/beta0 = {N}", out), v), 1e-4)

    def test_beta0_scales_as_one_over_alpha(self):
        b1, _, _, _ = self.speed("--alpha", 1, "--h", 0.99, "--l", 1e-4)
        b4, _, _, _ = self.speed("--alpha", 1e-4, "--h", 0.99, "--l", 1e-4)
        self.assertLess(rel(b4 / b1, 1e4), 1e-4)

    def test_h_equal_one_gives_zero_speed(self):
        _, _, v, _ = self.speed("--alpha", 1e-2, "--h", 1.0, "--l", 1e-4)
        self.assertEqual(v, 0.0)

    def test_deposition_same_magnitude(self):
        # h > 1: the front grows; the magnitude (what comp_eps needs) is symmetric in |h - 1|.
        _, _, vs, _ = self.speed("--alpha", 1e-2, "--h", 0.99, "--l", 1e-4)
        _, _, vd, _ = self.speed("--alpha", 1e-2, "--h", 1.01, "--l", 1e-4)
        self.assertLess(rel(vd, vs), 1e-4)

    def test_printed_comp_eps_command_runs_and_gives_pi_half(self):
        _, _, v, out = self.speed("--alpha", 1e-2, "--h", 0.99, "--l", 9.4e-5, "--Lx", 3.75e-4, "--Ly", 1.25e-4)
        line = next(x for x in out.splitlines() if x.startswith("comp_eps:"))
        argv = line.split(":", 1)[1].split()
        self.assertEqual(argv[1], "preprocess/comp_eps.py")
        r = subprocess.run([PY, *argv[1:], "--as-json"], capture_output=True, text=True, cwd=ROOT, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        import json
        res = json.loads(r.stdout)
        # K&P Eq. 45 resolution number at the chosen eps; the kinetic bound binds here, so Pi = safety = 0.5.
        pi = beta0(-20, 1e-2) * v * res["eps_m"] / d0(-20)
        self.assertEqual(res["eps_bounds_m"]["binding"], "B-KINETIC")
        self.assertLess(abs(pi - 0.5), 1e-3)

    def test_alpha_zero_fails_without_a_number(self):
        r = run("front_speed.py", "--alpha", 0, "--h", 0.99, "--l", 1e-4)
        self.assertNotEqual(r.returncode, 0)
        self.assertNotIn("v_n =", r.stdout)

    def test_non_numeric_argument_is_a_usage_error(self):
        r = run("front_speed.py", "--alpha", "abc", "--h", 0.99, "--l", 1e-4)
        self.assertEqual(r.returncode, 2)
        self.assertIn("invalid float value", r.stderr)

    # Fixed 2026-11-06 (was a bug): front_speed.py:47-48 accepts alpha_c > 1 and prints v_n (1.75e-9 m/s at
    # alpha 5); right: refuse alpha_c outside (0, 1] like mesh_params does.
    def test_alpha_above_one_is_refused(self):
        r = run("front_speed.py", "--alpha", 5, "--h", 0.99, "--l", 1e-4)
        self.assertNotEqual(r.returncode, 0)

    # Fixed 2026-11-06 (was a bug): front_speed.py:47-48 accepts h < 0 (sigma > 1, unphysical) and prints
    # v_n = 2.3e-7 m/s; right: refuse h < 0.
    def test_negative_humidity_is_refused(self):
        r = run("front_speed.py", "--alpha", 1e-2, "--h=-0.5", "--l", 1e-4)
        self.assertNotEqual(r.returncode, 0)

    # Fixed 2026-11-06 (was a bug): front_speed.py:47-48 accepts a negative path length; with
    # l = -1e-5 m it prints v_n larger than the kinetic upper bound; right: refuse l < 0.
    def test_negative_length_is_refused(self):
        r = run("front_speed.py", "--alpha", 1e-2, "--h", 0.99, "--l=-1e-5")
        self.assertNotEqual(r.returncode, 0)

    # Fixed 2026-11-06 (was a bug): front_speed.py:31-35 has no temperature range. At --T0 -173 (100 K,
    # lunar) it prints v_n = 3.6e-26 m/s from the ASHRAE fit and the AIR
    # diffusivity, and at +20 C a speed for ice above melting. Right: refuse
    # outside mesh_params' (-150, -0.001) C, as comp_eps does for the same input.
    def test_temperature_outside_formula_range_is_refused(self):
        for T in ("-173", "20"):
            r = run("front_speed.py", f"--T0={T}", "--alpha", 1e-2, "--h", 0.99, "--l", 1e-4)
            self.assertNotEqual(r.returncode, 0, T)


# ---------------------------------------------------------------------------
# supersat_compare.py
# ---------------------------------------------------------------------------
class SupersatCompare(unittest.TestCase):
    def test_all_sections_against_independent(self):
        T, R, H, L, gap = -20.0, 50e-6, 125e-6, 375e-6, 1e-4
        r = run("supersat_compare.py", "--temp", T, "--radius", R, "--width", H, "--theta", 30, 150,
                "--humidity", 0.99, "--xi-v", 1e-3, "--length", L, "--grad-T", 100, "--gap", gap)
        self.assertEqual(r.returncode, 0, r.stderr)
        o = r.stdout
        self.assertLess(rel(num(r"d0 = {N} m", o), d0(T)), 1e-4)
        self.assertLess(rel(num(r"rho_vs = {N} kg/m3", o), rho_vs(T)), 1e-4)
        self.assertLess(rel(num(r"D_v = {N} m2/s", o), d_v(T)), 1e-4)
        self.assertLess(rel(num(r"2D disc {N}", o), d0(T) / R), 1e-3)
        self.assertLess(rel(num(r"3D sphere {N}", o), 2 * d0(T) / R), 1e-3)
        m30 = d0(T) * 2 * math.cos(math.radians(30)) / H
        self.assertLess(rel(num(r"theta =  30.0: {N}", o), m30), 1e-3)
        self.assertLess(rel(num(r"theta = 150.0: {N}", o), -m30), 1e-3)
        self.assertLess(rel(num(r"\|sigma_h\|/sigma_GT\(3D\) = +{N}", o), 0.01 / (2 * d0(T) / R)), 1e-3)
        s = s_analytic(T)
        self.assertLess(rel(num(r"s = dln\(rho_vs\)/dT = {N}", o), s), 1e-5)
        self.assertLess(rel(num(r"G = 100 K/m: sigma_T = {N}", o), s * 100 * gap), 1e-3)
        self.assertLess(rel(num(r"G matching meniscus theta =  30.0 \(H = [^)]*\): {N}", o), m30 / (s * gap)), 1e-3)
        self.assertLess(rel(num(r"G matching sigma_h = [-+.e0-9]+: {N}", o), 0.01 / (s * gap)), 1e-3)
        self.assertLess(rel(num(r"rho_vs/\(xi rho_ice\) = {N}", o), rho_vs(T) / (1e-3 * RHO_ICE)), 1e-3)
        self.assertLess(rel(num(r"L\^2/\(xi D_v\) = {N}", o), L**2 / (1e-3 * d_v(T))), 1e-3)

    def test_gap_defaults_to_width(self):
        r = run("supersat_compare.py", "--width", 2e-4, "--grad-T", 10, "--theta", 60)
        self.assertIn("l = 0.0002 m, G = 10 K/m", r.stdout)

    # Fixed 2026-11-06 (was a bug): supersat_compare.py:84 means to skip theta = 90 (zero meniscus drive)
    # with `if sm != 0.0`, but cos(pi/2) = 6.1e-17, so it prints "G matching
    # meniscus theta = 90.0: 8.283e-17 K/m", a meaningless number. Right: skip
    # (test abs(sm) against a tolerance, or compute cos with the angle in degrees exact).
    def test_theta_90_has_no_matching_gradient_line(self):
        r = run("supersat_compare.py", "--theta", 90, "--grad-T", 1)
        self.assertEqual(r.returncode, 0)
        self.assertNotIn("G matching meniscus theta =  90.0", r.stdout)

    def test_zero_radius_and_zero_xi_fail(self):
        for args in (("--radius", 0), ("--xi-v", 0), ("--grad-T", 1, "--gap", 0)):
            with self.subTest(args=args):
                r = run("supersat_compare.py", *args)
                self.assertNotEqual(r.returncode, 0)


# ---------------------------------------------------------------------------
# tgrad_envelope.py
# ---------------------------------------------------------------------------
class TgradEnvelope(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        r = run("tgrad_envelope.py", "--T-lunar", 100, 150, "--G", 10, 1000, "--dT-open", 1, 10,
                "--Lx", 8.5e-4, "--h-open", 0.99)
        assert r.returncode == 0, r.stderr
        cls.o = r.stdout

    def test_murphy_koop_triple_point(self):
        # M&K 2005: 611.657 Pa at 273.16 K.
        self.assertAlmostEqual(num(r"273.16 K: {N} Pa", self.o), 611.66, places=2)
        self.assertLess(rel(p_mk(273.16), 611.657), 1e-5)

    def test_solver_rho_at_minus_20(self):
        self.assertAlmostEqual(num(r"solver rho_vs at -20 C: {N} g/m\^3", self.o), rho_vs(-20) * 1e3, places=4)

    def test_h_star(self):
        for dT in (1, 10):
            got = num(rf"dT +{dT:.1f}   h\* = {{N}}", self.o)
            self.assertAlmostEqual(got, rho_vs(-20) / rho_vs(-20 + dT), places=4)

    def test_cold_trapping_rows(self):
        rows = {}
        for line in self.o.splitlines():
            m = re.match(r"\s+(analog \(solver\)|lunar \(M&K\*?\))\s+([\d.]+)\s+(\S+)\s+(\S+)\s+([\d.]+)%\s+(\S+)\s+(\S+)", line)
            if m:
                rows[float(m.group(2))] = [float(m.group(i)) for i in (3, 4, 5, 6, 7)]
        self.assertEqual(sorted(rows), [100.0, 150.0, 253.2])
        k = s_analytic(-20.0)
        self.assertLess(rel(rows[253.2][0], k), 1e-4)
        for T in (100.0, 150.0):
            k_T = dln_rho_mk(T)
            dln, inv, corr, lt10, lt1000 = rows[T]
            with self.subTest(T=T):
                self.assertLess(rel(dln, k_T), 1e-4)
                self.assertLess(rel(inv, 1 / k_T), 1e-3)
                self.assertAlmostEqual(corr, 100 / (2 * T) / k_T, places=1)
                self.assertLess(rel(lt10, 1 / (10 * k_T)), 1e-3)
                self.assertLess(rel(lt1000, 1 / (1000 * k_T)), 1e-3)

    def test_extrapolation_flagged_below_110_K(self):
        lines = [x for x in self.o.splitlines() if x.strip().startswith("lunar (M&K")]
        self.assertTrue(any("M&K*" in x and " 100.0 " in x for x in lines))
        self.assertTrue(any("M&K)" in x and " 150.0 " in x for x in lines))

    def test_equivalent_drive_and_warming(self):
        k = s_analytic(-20.0)
        sat = num(r"G +1000 K/m: dT across channel [^>]*-> saturation change {N}", self.o)
        self.assertLess(rel(sat, 1000 * 8.5e-4 * k), 1e-3)
        ratio = num(r"G +1000 K/m: .*\({N} x the suction\)", self.o)
        self.assertAlmostEqual(ratio, 1000 * 8.5e-4 * k / 0.01, places=1)
        self.assertAlmostEqual(num(r"analog 253 K: {N}", self.o), rho_vs(-10) / rho_vs(-20), places=2)
        R, M = 8.314462618, 0.01801528
        rho = lambda T: p_mk(T) * M / (R * T)  # noqa: E731
        self.assertLess(rel(num(r"lunar  150 K: {N}", self.o), rho(160) / rho(150)), 2e-3)

    def test_bad_inputs_fail(self):
        for args in (("--h-open", 1.0), ("--T-lunar", 0), ("--G", 0)):
            with self.subTest(args=args):
                self.assertNotEqual(run("tgrad_envelope.py", *args).returncode, 0)


# ---------------------------------------------------------------------------
# plan_estimates.py
# ---------------------------------------------------------------------------
class PlanEstimates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        r = run("plan_estimates.py", "--G", 1, 100, 300, "--cutoff-h", 0.99)
        assert r.returncode == 0, r.stderr
        cls.o = r.stdout
        cls.s = s_analytic(-20.0)

    def test_quoted_planning_numbers(self):
        # Numbers quoted in the planning document.
        self.assertAlmostEqual(num(r"slope with D_v\(T0\) = \S+ s/m2 -> measured/predicted = {N}", self.o), 1.004, places=3)
        self.assertAlmostEqual(num(r"G\* = {N} K/m", self.o), 1.6, delta=0.05)
        self.assertAlmostEqual(num(r"h_open = 0.9900, G = +100.0 K/m: dz = +{N} mm", self.o), 1.05, delta=0.005)
        self.assertAlmostEqual(num(r"h_open = 0.9900, G = +300.0 K/m: dz = +{N} mm", self.o), 0.35, delta=0.005)

    def test_inputs_block(self):
        self.assertLess(rel(num(r"beta0 = {N} s/m", self.o), beta0(-20, 1e-2)), 1e-5)
        self.assertLess(rel(num(r"s = dln\(rho_vs\)/dT = {N}", self.o), self.s), 1e-5)
        self.assertLess(rel(num(r"K = {N} s/m2", self.o), k_diff(-20)), 1e-5)

    def test_b1_slope(self):
        pred = k_diff(-20) / 0.01
        self.assertLess(rel(num(r"slope with D_v\(T0\) = {N}", self.o), pred), 1e-4)
        self.assertLess(rel(num(r"measured/predicted = {N}", self.o), 5.73e12 / pred), 1e-4)
        self.assertLess(rel(num(r"intercept beta0/sigma = {N}", self.o), beta0(-20, 1e-2) / 0.01), 1e-4)

    def test_ripening_pair(self):
        Rs, Rl, L = 20e-6, 50e-6, 200e-6
        dsig = d0(-20) * (1 / Rs - 1 / Rl)
        self.assertLess(rel(num(r"G\* = {N} K/m", self.o), dsig / (self.s * L)), 1e-3)
        # 2D quasi-steady diffusion from a disc Rs to a ring at L: flux per length
        # 2 pi D drho / ln(L/Rs); shrink speed = flux / (2 pi Rs rho_ice).
        v_dif = d_v(-20) * rho_vs(-20) * dsig / (RHO_ICE * Rs * math.log(L / Rs))
        v_kin = dsig / beta0(-20, 1e-2)
        self.assertLess(rel(num(r"diffusion {N}", self.o), v_dif), 1e-3)
        self.assertLess(rel(num(r"kinetics {N}", self.o), v_kin), 1e-3)
        self.assertLess(rel(num(r"series {N} m/s", self.o), 1 / (1 / v_dif + 1 / v_kin)), 1e-3)

    def test_migration_speed_and_pi(self):
        v = d_v(-20) * rho_vs(-20) * self.s * 100 / RHO_ICE
        self.assertLess(rel(num(r"G = +100.0 K/m: v = {N} m/s", self.o), v), 1e-3)
        pi = beta0(-20, 1e-2) * v * 8.584e-7 / d0(-20)
        self.assertAlmostEqual(num(r"G = +100.0 K/m: .*Pi = beta0 v eps/d0 = {N}", self.o), pi, places=3)

    def test_harvest_time_by_integrating_the_front_ode(self):
        # dl/dt = sigma / (f (beta0 + K l)), integrated numerically from l0 to l1.
        b, K = beta0(-20, 1e-2), k_diff(-20)
        for h in (0.99, 0.9):
            l = np.linspace(1e-4, 0.85e-3, 200001)
            t = np.trapezoid((b + K * l) / (1 - h), l) / 86400
            self.assertAlmostEqual(num(rf"h = {h:5.3f}  l1 = 8.50e-04 m: t = +{{N}} days", self.o), t, places=1)

    def test_cutoff_depth_exact_by_bisection(self):
        # Solve rho_vs(T0 - G dz) = h rho_vs(T0) exactly, no linearization.
        for G in (100.0, 300.0):
            target = 0.99 * rho_vs(-20)
            lo, hi = 0.0, 1.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if rho_vs(-20 - G * mid) > target else (lo, mid)
            got = num(rf"h_open = 0.9900, G = +{G:.1f} K/m: dz = +{{N}} mm", self.o) * 1e-3
            self.assertLess(rel(got, lo), 3e-3, G)   # 3 decimals in mm; s varies little over 0.1 K

    def test_cutoff_h_open_one_is_zero(self):
        r = run("plan_estimates.py", "--G", 100, "--cutoff-h", 1.0)
        self.assertAlmostEqual(num(r"h_open = 1.0000, G = +100.0 K/m: dz = +{N} mm", r.stdout), 0.0)

    def test_zero_gradient_or_zero_alpha_fails(self):
        for args in (("--G", 0), ("--alpha", 0), ("--harvest", 1.0)):
            with self.subTest(args=args):
                r = run("plan_estimates.py", *args)
                self.assertNotEqual(r.returncode, 0)

    # Fixed 2026-11-06 (was a bug): plan_estimates.py:104-110 prints negative harvest times for h_open > 1
    # ("t = -24.2 days" at h 1.01; vapor deposits, nothing is harvested) and for
    # l1 < l0. Right: refuse (or skip with a message) when 1 - h <= 0 or l1 <= l0.
    def test_harvest_with_deposition_is_not_a_negative_time(self):
        r = run("plan_estimates.py", "--harvest", 1.01, "--G", 1)
        days = nums(r"h = 1.010 .*t = +{N} days", r.stdout)
        self.assertTrue(r.returncode != 0 or (days and min(days) >= 0), r.stdout)

    # Fixed 2026-11-06 (was a bug): plan_estimates.py:83-94 accepts a ripening pair that overlaps
    # (center distance 10 um < Rs + Rl = 70 um, printed gap -6e-5 m) or with
    # L < Rs (log(L/Rs) < 0), and prints a shrink speed of 835 um/day. Right:
    # refuse gap <= 0 and R_small >= R_large.
    def test_overlapping_pair_is_refused(self):
        r = run("plan_estimates.py", "--pair", 20e-6, 50e-6, 10e-6, "--G", 1)
        self.assertNotEqual(r.returncode, 0)


# ---------------------------------------------------------------------------
# analyze_stage0.py, on synthetic batches with exact answers
# ---------------------------------------------------------------------------
A_EXP = "openface_T-20_curvflow_22h_ac1e-2"
A_INT = "openface_2D_box60um_disc_interior"
A_OFF = "openface_2D_box60um_disc_offcenter"
A_CAP = "openface_2D_box60x30um_cap"
B_GEOM = "channel_2D_H125um_eps0.39um_deadend"
B_EXP = {"B1": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-3", "B2": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-4",
         "B3": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-2", "B4": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-3_xiT1e-3"}
C_GEOM, C_EXP = "throat_2D_L850um_eps0.86um_narrowL_deadend", "cryo_T-20_h0.990_iso_5d_ac1e-2"

D_EFF, R0, D_OFF, EPS = 1.14e-15, 15e-6, 7.5e-6, 6.4e-7
TH, R0C, CX = math.pi / 3, 20e-6, -5e-6


def clipped_circle_area(cx, cy, R, quadrant=False, n=20000):
    """Area of a circle inside x > 0 (and y > 0 if quadrant), by clipping a
    fine polygon (Sutherland-Hodgman) and the shoelace formula."""
    if R <= 0:
        return 0.0
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    poly = list(zip(cx + R * np.cos(a), cy + R * np.sin(a)))

    def clip(pts, axis):
        out = []
        for i in range(len(pts)):
            p, q = pts[i - 1], pts[i]
            pin, qin = p[axis] >= 0, q[axis] >= 0
            if qin:
                if not pin:
                    s = p[axis] / (p[axis] - q[axis])
                    out.append((p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])))
                out.append(q)
            elif pin:
                s = p[axis] / (p[axis] - q[axis])
                out.append((p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])))
        return out
    poly = clip(poly, 0)
    if quadrant and poly:
        poly = clip(poly, 1)
    if len(poly) < 3:
        return 0.0
    x, y = np.array(poly).T
    # Polygon inscribed in the circle: correct the chord deficit to O(1/n^2) by area scaling.
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / (np.sin(2 * np.pi / n) / (2 * np.pi / n))


def a_law_disc(t, D=D_EFF):
    return math.pi * (R0**2 - 2 * D * t)


def a_law_offcenter(t, D=D_EFF):
    R2 = R0**2 - 2 * D * t
    return clipped_circle_area(-D_OFF, 0.0, math.sqrt(max(R2, 0.0)))


def a_law_cap(t, D=D_EFF):
    shape = TH - math.sin(TH) * math.cos(TH)
    R = math.sqrt(max(R0C**2 * shape - 2 * TH * D * t, 0.0) / shape)
    return clipped_circle_area(CX, -R * math.cos(TH), R, quadrant=True)


def analyze(batch, *extra):
    return subprocess.run([PY, str(SCRIPTS / "analyze_stage0.py"), "--batch", str(batch), *map(str, extra)],
                          capture_output=True, text=True, cwd=ROOT, timeout=60)


def verdict(out, tag):
    m = re.search(rf"  {re.escape(tag)}: max rel. error (\S+) over (\d+) rows .*?(PASS|FAIL)", out)
    if m is None:
        raise AssertionError(f"{tag} not reported:\n{out}")
    return float(m.group(1)), int(m.group(2)), m.group(3)


class ClippedAreaSelfCheck(unittest.TestCase):
    def test_clipper_against_closed_forms(self):
        R = 15e-6
        self.assertLess(rel(clipped_circle_area(1.0, 1.0, 0.5), math.pi * 0.25), 1e-6)
        # Segment beyond x = d of a disc centered at -d: 2 * int_d^R sqrt(R^2 - u^2) du.
        d = 7.5e-6
        exact = R * R * (math.pi / 2 - math.asin(d / R)) - d * math.sqrt(R * R - d * d)
        self.assertLess(rel(clipped_circle_area(-d, 0.0, R), exact), 1e-6)
        # Quarter of a disc centered at the origin.
        self.assertLess(rel(clipped_circle_area(0.0, 0.0, R, quadrant=True), math.pi * R * R / 4), 1e-6)


class AnalyzeStage0A(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.b = self.tmp.name
        self.t = np.linspace(0, 7.0e4, 141)

    def tearDown(self):
        self.tmp.cleanup()

    def a3(self, layout="hpc", law=a_law_disc, t=None, tag=""):
        t = self.t if t is None else t
        r = F.make_run(self.b, A_INT, A_EXP, layout=layout, geometry={"-eps": f"{EPS:e}"}, tag=tag)
        F.write_ssa(r, t, [law(x) for x in t])
        F.write_outp(r, range(len(t)), t, tau_sub=EPS**2 / D_EFF)
        return r

    def other(self, geom, exp, law, layout="hpc"):
        r = F.make_run(self.b, geom, exp, layout=layout, geometry={"-eps": f"{EPS:e}"})
        F.write_ssa(r, self.t, [law(x) for x in self.t])
        return r

    def test_exact_curvature_flow_passes_and_neumann_control_fails(self):
        self.a3()
        self.other(A_OFF, A_EXP, a_law_offcenter)
        self.other(A_OFF, A_EXP + "_closed", lambda t: a_law_offcenter(t, D=0.7 * D_EFF))
        self.other(A_CAP, A_EXP, a_law_cap)
        r = analyze(self.b)
        self.assertEqual(r.returncode, 0, r.stderr)
        o = r.stdout
        self.assertLess(rel(num(r"D_eff = {N} m2/s", o), D_EFF), 1e-4)
        self.assertLess(rel(num(r"R0 = {N} m", o), R0), 1e-4)
        self.assertLess(rel(num(r"eps {N}", o), EPS), 1e-6)
        err, rows, v = verdict(o, "A1 open")
        self.assertEqual(v, "PASS")
        self.assertLess(err, 1e-4)
        # Rows with R > d + 4 eps: R^2 = R0^2 - 2 D t > (d + 4 eps)^2.
        self.assertEqual(rows, int(np.sum(R0**2 - 2 * D_EFF * self.t > (D_OFF + 4 * EPS) ** 2)))
        self.assertEqual(verdict(o, "A2 90deg")[2], "FAIL")
        err, _, v = verdict(o, "A4 cap at the corner")
        self.assertEqual(v, "PASS")
        self.assertLess(err, 1e-3)

    def test_perturbed_open_face_fails(self):
        self.a3()
        self.other(A_OFF, A_EXP, lambda t: a_law_offcenter(t) * (1 + 0.03 * t / 7e4))
        err, _, v = verdict(analyze(self.b).stdout, "A1 open")
        self.assertEqual(v, "FAIL")
        self.assertGreater(err, 0.02)

    def test_tolerance_flag(self):
        self.a3()
        self.other(A_OFF, A_EXP, lambda t: a_law_offcenter(t) * 1.01)
        self.assertEqual(verdict(analyze(self.b, "--tol", 0.02).stdout, "A1 open")[2], "PASS")
        self.assertEqual(verdict(analyze(self.b, "--tol", 0.005).stdout, "A1 open")[2], "FAIL")

    def test_local_layout_found_with_stage0_tag(self):
        self.a3(layout="local", tag="stage0_A3")
        self.assertIn("D_eff =", analyze(self.b).stdout)

    def test_local_layout_with_other_tag_is_not_matched(self):
        self.a3(layout="local", tag="v2check")
        self.assertIn("A3 run not found", analyze(self.b).stdout)

    def test_duplicate_local_runs_newest_wins(self):
        old = self.a3(layout="local")
        new = old.replace("2026-01-01__00.00.00", "2026-02-01__00.00.00")
        os.rename(old, new)
        old = self.a3(layout="local", law=lambda t: a_law_disc(t, D=2 * D_EFF))
        o = analyze(self.b).stdout
        self.assertIn("A3 2026-02-01__00.00.00_", o)
        self.assertLess(rel(num(r"D_eff = {N}", o), D_EFF), 1e-4)

    def test_both_layouts_present_hpc_wins_silently(self):
        # Characterization, not a correctness claim: with one run in each
        # layout, the HPC folder sorts last and is used, with no warning that
        # a second candidate exists.
        self.a3(layout="local", law=lambda t: a_law_disc(t, D=2 * D_EFF))
        self.a3(layout="hpc")
        o = analyze(self.b).stdout
        self.assertIn(f"A3 {A_INT}__{A_EXP}:", o)
        self.assertNotIn("duplicate", o.lower())

    def test_crashed_a3_still_gives_d_eff(self):
        self.a3(t=self.t[:30])
        self.assertLess(rel(num(r"D_eff = {N}", analyze(self.b).stdout), D_EFF), 1e-4)

    def test_nan_rows_never_pass(self):
        self.a3()
        r = self.other(A_OFF, A_EXP, a_law_offcenter)
        a = [a_law_offcenter(x) for x in self.t]
        a[10] = float("nan")
        F.write_ssa(r, self.t, a)
        res = analyze(self.b)
        self.assertNotEqual(verdict(res.stdout, "A1 open")[2], "PASS")

    def test_single_row_a3_never_passes(self):
        # A run that crashed after its first row cannot give a slope. Today this
        # raises numpy LinAlgError (exit 1) and stops parts B and C too.
        self.a3(t=self.t[:1])
        self.other(A_OFF, A_EXP, a_law_offcenter)
        res = analyze(self.b)
        self.assertNotIn("PASS", res.stdout)

    def test_missing_eps_gives_fail_not_pass(self):
        r = F.make_run(self.b, A_INT, A_EXP, layout="hpc")
        F.write_ssa(r, self.t, [a_law_disc(x) for x in self.t])
        self.other(A_OFF, A_EXP, a_law_offcenter)
        err, rows, v = verdict(analyze(self.b).stdout, "A1 open")
        self.assertEqual((rows, v), (0, "FAIL"))

    # Fixed 2026-11-06 (was a bug): analyze_stage0.py:43-46 reads -eps from the staged *.opts in
    # ALPHABETICAL order, last one wins. The solver reads solver.opts <
    # geometry < experiment, and "solver.opts" sorts last, so an -eps in
    # solver.opts would override the geometry's. Today no inputs/solver.opts sets
    # -eps, so stage-0 numbers are unaffected. Right: read in the solver's order
    # (as pplib now does, commit 18f4bf66) and honor "Extra opts" from outp.txt.
    def test_eps_follows_solver_precedence(self):
        r = F.make_run(self.b, A_INT, A_EXP, layout="hpc", solver={"-eps": "1.0e-6"},
                       geometry={"-eps": f"{EPS:e}"})
        F.write_ssa(r, self.t, [a_law_disc(x) for x in self.t])
        F.write_outp(r, range(len(self.t)), self.t, tau_sub=1.0)
        self.assertLess(rel(num(r", eps {N},", analyze(self.b).stdout), EPS), 1e-6)

    def test_figure_written(self):
        self.a3()
        self.other(A_OFF, A_EXP, a_law_offcenter)
        out = Path(self.b) / "figA.png"
        r = analyze(self.b, "--fig-a", out)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertGreater(out.stat().st_size, 10000)


class AnalyzeStage0BC(unittest.TestCase):
    T_END, RATE, A0, TAU = 1.0e4, 1.0e-13, 1.0e-8, 100.0

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.b = self.tmp.name
        self.t = np.linspace(0, self.T_END, 101)

    def tearDown(self):
        self.tmp.cleanup()

    def brun(self, key, t=None, rate=None, flux=True, geom=B_GEOM, exp=None):
        t = self.t if t is None else t
        rate = self.RATE if rate is None else rate
        r = F.make_run(self.b, geom, exp or B_EXP[key], layout="hpc")
        F.write_ssa(r, t, self.A0 - rate * t)
        if flux:
            F.write_bflux(r, t, {"x0": (RHO_ICE * rate, 0.0), "x1": (0.0, 0.0)})
        F.write_outp(r, range(len(t)), t, tau_sub=self.TAU)
        return r

    def test_crashed_run_compared_over_common_window(self):
        self.brun("B1")
        self.brun("B2", t=self.t[:40])                        # crashed at 3.9e3 s, same law
        self.brun("B3", rate=1.05 * self.RATE)                # 5 % off
        o = analyze(self.b).stdout
        self.assertLess(rel(num(r"B1: ice change {N} m2", o), -self.RATE * self.T_END), 1e-5)
        self.assertLess(rel(num(r"mean J_x0 {N} kg/\(m s\)", o), RHO_ICE * self.RATE), 1e-6)
        m = re.search(r"B2 xi_v 1e-4: over t <= (\S+) s, ice change differs from B1 by (\S+);.*(PASS|FAIL)", o)
        self.assertAlmostEqual(float(m.group(1)), self.t[39], delta=1.0)
        self.assertLess(float(m.group(2)), 1e-5)
        self.assertEqual(m.group(3), "PASS")
        m = re.search(r"B3 xi_v 1e-2: .* differs from B1 by (\S+);.*\((\S+)\)  (PASS|FAIL)", o)
        self.assertAlmostEqual(float(m.group(1)), 0.05, places=4)
        self.assertAlmostEqual(float(m.group(2)), 0.05, places=4)       # mean J differs by 5 % too
        self.assertEqual(m.group(3), "FAIL")

    def test_ic_transient_on_a_linear_law(self):
        self.brun("B1")
        o = analyze(self.b, "--t-transient", 1000).stdout
        tr = num(r"ice change in the first 1000 s = {N} m2", o)
        late = num(r"would give {N} m2", o)
        self.assertLess(rel(tr, -self.RATE * 1000), 1e-5)
        self.assertLess(rel(late, -self.RATE * 1000), 1e-5)

    def test_solver_cost_totals(self):
        r = self.brun("B1")
        steps = np.arange(51)
        F.write_solver_evo(r, steps, steps * 10.0, newton=2, krylov=30, rej=1)
        o = analyze(self.b).stdout
        self.assertIn("51 steps, 102 Newton its (2.00/step), 1530 Krylov its, 51 rejected", o)

    def test_resolution_window_and_flags(self):
        r = self.brun("B1")
        eps, b0, dd = 3.9358e-7, 7.9408e5, 1.01517e-9
        vn_for = lambda pi: pi * dd / (b0 * eps)  # noqa: E731
        # Pi = 2 inside the 2 tau_sub = 200 s start-up window only: ignored.
        F.write_resolution(r, self.t, np.where(self.t < 200, vn_for(2.0), vn_for(0.3)), eps=eps, beta0=b0, d0=dd)
        o = analyze(self.b).stdout
        self.assertAlmostEqual(num(r"after the {N} s start-up window", o), 200.0)
        self.assertAlmostEqual(num(r"max Pi {N}", o), 0.3, places=3)
        self.assertNotIn("CHECK", o)
        vn = np.full(self.t.size, vn_for(0.3))
        vn[50] = vn_for(1.5)
        F.write_resolution(r, self.t, vn, eps=eps, beta0=b0, d0=dd)
        o = analyze(self.b).stdout
        self.assertIn("(steps > 0.5: 1, > 1: 1)", o)
        self.assertIn("<-- CHECK", o)

    def test_nan_row_in_b2_never_passes(self):
        self.brun("B1")
        r = self.brun("B2")
        a = self.A0 - self.RATE * self.t
        a[-1] = float("nan")
        F.write_ssa(r, self.t, a)
        o = analyze(self.b).stdout
        self.assertNotRegex(o, r"B2 xi_v 1e-4: .*PASS")

    # Fixed 2026-11-06 (was a bug): analyze_stage0.py:252-257: when a B2-B4 run (or B1) has no bflux.dat,
    # bflux() returns (None, None) and mean_j() does `None <= tc`, so the whole
    # analysis dies with TypeError and part C never runs. Line 232 already
    # handles a missing B1 file with nan. Right: report nan for the mean J of
    # that run and carry on.
    def test_missing_bflux_in_b2_does_not_crash(self):
        self.brun("B1")
        self.brun("B2", flux=False)
        self.brun(None, geom=C_GEOM, exp=C_EXP)
        r = analyze(self.b)
        self.assertEqual(r.returncode, 0, r.stderr[-400:])
        self.assertIn("== C: pilot ==", r.stdout)

    def test_mass_balance_ratio_is_one_for_exact_flux(self):
        self.brun(None, geom=C_GEOM, exp=C_EXP)
        o = analyze(self.b).stdout
        self.assertAlmostEqual(num(r"ratio {N}\)", o), 1.0, places=4)
        self.assertAlmostEqual(num(r"ice left {N} of initial", o), 1 - self.RATE * self.T_END / self.A0, places=3)
        self.assertLess(rel(num(r"mass out through x0 {N} kg/m", o), RHO_ICE * self.RATE * self.T_END), 1e-6)

    def test_c_without_bflux_reports_nan(self):
        self.brun(None, geom=C_GEOM, exp=C_EXP, flux=False)
        r = analyze(self.b)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("mass out through x0 nan", r.stdout)


# ---------------------------------------------------------------------------
# plot_setup.py: synthetic IC folders; layout checked with postprocess/figcheck.py
# ---------------------------------------------------------------------------
IC_DIMS = {"chan": (375e-6, 125e-6), "throat": (850e-6, 180e-6), "regpore": (814e-6, 260e-6),
           "off": (60e-6, 45e-6), "int": (60e-6, 45e-6), "cap": (60e-6, 30e-6)}


def make_ic_root(root, dims=IC_DIMS):
    for name, (Lx, Ly) in dims.items():
        run_dir = os.path.join(root, name)
        os.makedirs(run_dir)
        X, Y = F.grid(Lx, Ly, 80, 30)
        phi = F.disc(X, Y, Lx / 2, Ly / 2, Ly / 4, Ly / 60)
        F.write_vts_run(run_dir, X, Y, [(0, 1e-9, phi, 0 * X - 20, 0 * X)])


class PlotSetup(unittest.TestCase):
    def test_cli_writes_both_figures(self):
        with tempfile.TemporaryDirectory() as d:
            make_ic_root(d)
            r = run("plot_setup.py", "--ic-root", d, "--out", os.path.join(d, "out"))
            self.assertEqual(r.returncode, 0, r.stderr)
            for f in ("stage0_setup_cryosuction.png", "stage0_setup_verification.png"):
                self.assertGreater(os.path.getsize(os.path.join(d, "out", f)), 50000)

    def test_missing_ic_folder_fails(self):
        with tempfile.TemporaryDirectory() as d:
            make_ic_root(d, {k: v for k, v in IC_DIMS.items() if k != "regpore"})
            r = run("plot_setup.py", "--ic-root", d, "--out", os.path.join(d, "out"))
            self.assertNotEqual(r.returncode, 0)

    def _layout_problems(self):
        """Run both figure builders in this process and collect figcheck results."""
        sys.path.insert(0, str(SCRIPTS))
        sys.path.insert(0, str(ROOT / "postprocess"))
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.figure import Figure
        from figcheck import check_figure
        import plot_setup
        found = {}
        orig = Figure.savefig

        def capture(fig, fname, *a, **k):
            found[Path(fname).name] = check_figure(fig, min_pt=7.0)
            return orig(fig, fname, *a, **k)
        Figure.savefig = capture
        try:
            with tempfile.TemporaryDirectory() as d:
                make_ic_root(d)
                out = Path(d) / "out"
                out.mkdir()
                plot_setup.main_figure(Path(d), out)
                plot_setup.verification_figure(Path(d), out)
        finally:
            Figure.savefig = orig
        return found

    def test_verification_figure_layout_clean(self):
        self.assertEqual(self._layout_problems()["stage0_setup_verification.png"], [])

    # Fixed 2026-11-06 (was a bug): plot_setup.py:109-113: the three-line subtitle (fig.text at y 0.945)
    # and the legend (anchored at y 0.915) are placed at fixed figure fractions;
    # figcheck finds the subtitle's last line overlapping the legend's first row
    # (in the committed PNG the two touch with no gap). Right: no text/text
    # overlap (move the legend down or derive its y from the subtitle's extent).
    def test_main_figure_layout_clean(self):
        self.assertEqual(self._layout_problems()["stage0_setup_cryosuction.png"], [])


# ---------------------------------------------------------------------------
# preflight.sh: static checks, and a COPY run in a temporary tree with a stub binary
# ---------------------------------------------------------------------------
STUB = """#!/bin/bash
# Stub standing in for the solver: never computes anything.
echo "$@" >> "$STUB_LOG"
echo "   tau_sub  =  1.2000e+02 s"
echo "   dtmax/tau_sub = 2.0002   (ok, <= 5.0)"
echo "   h_x = L/Nx = 1.0e-07 m"
for f in "$@"; do
  [[ -f $f ]] || continue
  grep -q STUB_WARN "$f" && echo "   RESOLUTION WARNING: stub"
  grep -q STUB_CRASH "$f" && { echo "stub crash"; exit 3; }
done
exit 0
"""


class Preflight(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        t = Path(self.tmp.name)
        self.root = t / "tree"
        sd = self.root / "studies" / "cryosuction" / "scripts"
        sd.mkdir(parents=True)
        for name in ("preflight.sh", "check_open_bc.sh"):
            if (SCRIPTS / name).exists():
                shutil.copy2(SCRIPTS / name, sd / name)
        (self.root / "scripts" / "lib").mkdir(parents=True)
        shutil.copy2(ROOT / "scripts" / "lib" / "opts.sh", self.root / "scripts" / "lib" / "opts.sh")
        for sub in ("geometry/fam", "experiment/fam", "experiment/other"):
            (self.root / "inputs" / sub).mkdir(parents=True)
        (self.root / "inputs" / "solver.opts").write_text("-flag_BC_rhovfix 0\n")
        for g in ("g1", "gwarn", "gcrash"):
            (self.root / "inputs" / "geometry" / "fam" / f"{g}.opts").write_text(
                {"g1": "-dim 2\n", "gwarn": "# STUB_WARN\n", "gcrash": "# STUB_CRASH\n"}[g])
        for e in ("e1", "e2", "dup"):
            (self.root / "inputs" / "experiment" / "fam" / f"{e}.opts").write_text("-t_final 1\n")
        (self.root / "inputs" / "experiment" / "other" / "dup.opts").write_text("-t_final 1\n")
        self.log = t / "stub.log"
        self.work = t / "tmpdir"
        self.work.mkdir()
        self.script = sd / "preflight.sh"

    def tearDown(self):
        self.tmp.cleanup()

    def add_stub(self):
        exe = self.root / "lunar_cryosuction"
        exe.write_text(STUB)
        exe.chmod(exe.stat().st_mode | stat.S_IXUSR)

    def preflight(self, text, name="tests.txt"):
        f = Path(self.tmp.name) / name
        if text is not None:
            f.write_text(text)
        env = {**os.environ, "TMPDIR": str(self.work), "STUB_LOG": str(self.log)}
        return subprocess.run(["bash", str(self.script), str(f)], capture_output=True, text=True, env=env, timeout=60)

    def calls(self):
        return self.log.read_text().splitlines() if self.log.exists() else []

    def test_missing_binary_exits_2_with_message(self):
        r = self.preflight("g1:e1\n")
        self.assertEqual(r.returncode, 2)
        self.assertIn("lunar_cryosuction -- run make", r.stderr)
        self.assertEqual(self.calls(), [])

    def test_comments_blanks_and_spaces(self):
        self.add_stub()
        r = self.preflight("# header\n\n  g1 : e1   # trailing comment\n\t\n")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("OK   g1:e1", r.stdout)
        self.assertIn("tau_sub  =  1.2000e+02 s", r.stdout)
        calls = self.calls()
        self.assertEqual(len(calls), 1)
        args = calls[0].split()
        # Solver order: solver.opts < geometry < experiment, then -t_final 0.
        files = [args[i + 1] for i, a in enumerate(args) if a == "-options_file"]
        self.assertEqual(files[0], "inputs/solver.opts")
        self.assertTrue(files[1].endswith("geometry/fam/g1.opts") and files[2].endswith("experiment/fam/e1.opts"))
        self.assertEqual(args[args.index("-t_final") + 1], "0")

    def test_unresolved_and_ambiguous_names_fail(self):
        self.add_stub()
        r = self.preflight("g1:nope\nnogeom:e1\ng1:dup\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("FAIL g1:nope: experiment not found", r.stdout)
        self.assertIn("FAIL nogeom:e1: geometry not found", r.stdout)
        self.assertIn("FAIL g1:dup: experiment not found", r.stdout)     # ambiguous name
        self.assertEqual(self.calls(), [])

    def test_resolution_warning_and_crash_fail(self):
        self.add_stub()
        r = self.preflight("gwarn:e1\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("RESOLUTION WARNING", r.stdout)
        r = self.preflight("gcrash:e1\ng1:e1\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("FAIL gcrash:e1", r.stdout)
        self.assertIn("stub crash", r.stdout)
        self.assertIn("OK   g1:e1", r.stdout)         # a failure does not stop later lines

    def test_missing_tests_file_is_not_ok(self):
        self.add_stub()
        r = self.preflight(None, name="nonexistent.txt")
        self.assertNotEqual(r.returncode, 0)

    def test_work_folder_removed(self):
        self.add_stub()
        self.preflight("g1:e1\n")
        self.assertEqual(list(self.work.iterdir()), [])

    # Fixed 2026-11-06: preflight.sh and check_open_bc.sh dropped a last line
    # with no trailing newline ("g1:e1\\ng1:e2" checked only g1:e1).
    def test_last_line_without_newline_is_checked(self):
        self.add_stub()
        r = self.preflight("g1:e1\ng1:e2")
        self.assertEqual(len(self.calls()), 2, r.stdout)

    def test_empty_tests_file_passes_vacuously(self):
        # Characterization: an empty (or all-comment) tests file exits 0 with
        # nothing checked. Worth knowing before trusting a green preflight.
        self.add_stub()
        r = self.preflight("# nothing\n")
        self.assertEqual(r.returncode, 0)
        self.assertEqual(self.calls(), [])


class CommittedTestsFilesResolve(unittest.TestCase):
    """Every geom:exp line in the committed tests files resolves to exactly one
    opts file in the live tree (the part of preflight that needs no binary)."""

    def test_resolve(self):
        files = sorted((ROOT / "studies" / "cryosuction").glob("*_tests.txt"))
        self.assertTrue(files)
        for f in files:
            for raw in f.read_text().splitlines():
                line = re.sub(r"\s", "", raw.split("#", 1)[0])
                if not line:
                    continue
                g, e = line.split(":", 1)
                for kind, name in (("geometry", g), ("experiment", e)):
                    with self.subTest(file=f.name, line=line, kind=kind):
                        r = subprocess.run(["bash", "-c", 'source "$1"; resolve_opts "$2" "$3"', "_",
                                            str(ROOT / "scripts" / "lib" / "opts.sh"),
                                            str(ROOT / "inputs" / kind), name],
                                           capture_output=True, text=True, timeout=30)
                        if r.returncode != 0 and list((ROOT / "inputs" / "scratch" / kind).rglob(f"{name}.opts")):
                            # Retired to inputs/scratch: the tests file is a record
                            # of a past batch and cannot be resubmitted as is.
                            self.skipTest(f"{name} is retired to inputs/scratch/{kind}")
                        self.assertEqual(r.returncode, 0, r.stderr)
                        self.assertTrue(r.stdout.endswith(f"/{name}.opts"))


if __name__ == "__main__":
    unittest.main()
