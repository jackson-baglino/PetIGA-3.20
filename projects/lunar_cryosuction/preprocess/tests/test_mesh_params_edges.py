"""Edge-case tests for preprocess/mesh_params and preprocess/comp_eps.py.

test_mesh_params.py already pins the golden numbers, band counts, option set,
input-error listing, the alpha_c sweep and the constants against the C source.
This file adds what it does not cover: the eps bounds recomputed independently
from the physics (constants read from the C source, not from mesh_params), the
validity limits, the mesh rounding, the CLI failure modes and --patch.

Run from the project root:
    venv_cryosuction/bin/python -m unittest discover -s preprocess/tests -t . -v

Tests marked expectedFailure document a bug ("# BUG:" comment above each).
"""
import json
import math
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

from preprocess.mesh_params import MeshParamsInput, ToolInputError, compute_mesh_params

ROOT = Path(__file__).resolve().parents[2]
PY = str(ROOT / "venv_cryosuction" / "bin" / "python")

# Independent physics: constants from the C source / K&P, not from mesh_params.
_C = (ROOT / "src" / "material_properties.c").read_text()
K_ASH = [float(re.search(rf"\bK{k}\s*=\s*([-+0-9.eE]+)", _C).group(1)) for k in range(6)]
RHO_ICE, RHO_AIR, P_ATM, BB, DV0 = 919.0, 1.341, 101325.0, 0.62, 2.178e-5
M_HK, KB_HK, M_D0, KB_D0, GAMMA = 2.99e-26, 1.38e-23, 2.9915e-26, 1.380649e-23, 0.109
D_HEAT_MEAN = 0.5 * (2.29 / 1.8e6 + 0.02 / (1.341 * 1004.0))   # (k_i/C_i + k_a/C_a) / 2


def rho_vs(T):
    t = T + 273.15
    k = K_ASH
    p = math.exp(k[0] / t + k[1] + k[2] * t + k[3] * t**2 + k[4] * t**3 + k[5] * math.log(t))
    return RHO_AIR * BB * p / (P_ATM - p)


def beta_unscaled(T, alpha):
    return (1 / alpha) * math.sqrt(2 * math.pi * M_HK / (KB_HK * (T + 273.15))) * RHO_ICE / rho_vs(T)


def d0(T):
    return GAMMA * (M_D0 / RHO_ICE) / (KB_D0 * (T + 273.15))


def dv(T):
    return DV0 * ((T + 273.15) / 273.15) ** 1.81


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def cli(*args):
    return subprocess.run([PY, str(ROOT / "preprocess" / "comp_eps.py"), *map(str, args)],
                          capture_output=True, text=True, cwd=ROOT, timeout=60)


class IndependentBounds(unittest.TestCase):
    CASES = [  # (T, alpha_c, v_n, R, Lx)
        (-20.0, 1e-2, 1.6e-9, 3e-5, 3.75e-4),
        (-20.0, 1.0, 1e-9, 3e-5, 1e-4),          # alpha_c at its upper limit
        (-20.0, 1e-4, 1e-9, 3e-5, 1e-4),
        (-149.0, 1e-2, 1e-12, 1e-5, 1e-4),       # cold end of the accepted range
        (-0.01, 1e-3, 1e-8, 1e-4, 1e-3),         # warm end
    ]

    def test_bounds_eps_and_mesh(self):
        for T, a, v, R, Lx in self.CASES:
            with self.subTest(T=T, alpha=a):
                r = compute_mesh_params(MeshParamsInput(domain_x_m=Lx, grain_radius_m=R, temperature_c=T,
                                                        alpha_c=a, front_velocity_m_per_s=v))
                b_uns = beta_unscaled(T, a)
                b_hk = b_uns * rho_vs(T) / RHO_ICE           # scaled K&P beta'
                exp = {"B-KINETIC": d0(T) / (b_uns * v), "B-VAPOR": dv(T) * b_hk,
                       "B-HEAT": D_HEAT_MEAN * b_hk, "B-CURV": R}
                for k, val in exp.items():
                    self.assertLess(rel(r["eps_bounds_m"][k], val), 1e-9, k)
                self.assertEqual(r["eps_bounds_m"]["binding"], min(exp, key=exp.get))
                eps = 0.5 * min(exp.values())
                self.assertLess(rel(r["eps_m"], eps), 1e-9)
                self.assertEqual(r["mesh"]["nx"], math.ceil(Lx * math.sqrt(2) / eps - 1e-9))

    def test_resolution_number_is_safety_when_kinetic_binds(self):
        # K&P Eq. 45: Pi = beta0 v eps / d0; with eps = safety * d0/(beta0 v), Pi = safety.
        for safety in (0.5, 0.25, 1.0):
            r = compute_mesh_params(MeshParamsInput(domain_x_m=1e-4, alpha_c=1e-2, front_velocity_m_per_s=1e-8,
                                                    safety=safety))
            self.assertEqual(r["eps_bounds_m"]["binding"], "B-KINETIC")
            pi = beta_unscaled(-20, 1e-2) * 1e-8 * r["eps_m"] / d0(-20)
            self.assertAlmostEqual(pi, safety, places=9)

    def test_feature_radius_makes_kinetic_bound_equal_the_radius(self):
        # v_n = d0/(beta0 R_feat) gives B-KINETIC = d0/(beta0 v_n) = R_feat (the circularity noted in memory).
        r = compute_mesh_params(MeshParamsInput(domain_x_m=1e-4, alpha_c=1e-3, feature_radius_m=7e-6))
        self.assertLess(rel(r["eps_bounds_m"]["B-KINETIC"], 7e-6), 1e-12)

    def test_sigma_surf_uses_libbrecht_table(self):
        # Libbrecht sigma0(-20 C) = 3.5e-2 (table entry); alpha = exp(-sigma0/sigma).
        r = compute_mesh_params(MeshParamsInput(domain_x_m=1e-4, sigma_surf=0.01))
        self.assertLess(rel(r["inputs"]["alpha_c_used"], math.exp(-3.5)), 1e-12)


class Limits(unittest.TestCase):
    def ok(self, **kw):
        return compute_mesh_params(MeshParamsInput(**{"domain_x_m": 1e-4, "alpha_c": 1e-2, **kw}))

    def test_temperature_open_interval(self):
        for T in (-149.99, -0.0011):
            self.ok(temperature_c=T)
        for T in (-150.0, -173.15, -0.001, 0.0, 5.0, float("nan")):
            with self.subTest(T=T), self.assertRaises(ToolInputError):
                self.ok(temperature_c=T)

    def test_alpha_limits(self):
        self.ok(alpha_c=1.0)
        self.ok(alpha_c=1e-4)
        for a in (0.0, -1e-3, 1.0000001, float("nan"), float("inf")):
            with self.subTest(a=a), self.assertRaises(ToolInputError):
                self.ok(alpha_c=a)

    def test_velocity_and_lengths(self):
        for kw in ({"front_velocity_m_per_s": 0.0}, {"front_velocity_m_per_s": -1e-9}, {"domain_x_m": 0.0},
                   {"grain_radius_m": 0.0}, {"domain_z_m": 1e-4}, {"safety": 1.01}, {"xi_v": 0.0}):
            with self.subTest(kw=kw), self.assertRaises(ToolInputError):
                self.ok(**kw)

    def test_tiny_and_huge_domains(self):
        self.assertEqual(self.ok(domain_x_m=1e-12)["mesh"]["nx"], 1)
        r = self.ok(domain_x_m=1.0, domain_y_m=1e-12)
        self.assertEqual(r["mesh"]["ny"], 1)
        self.assertEqual(r["mesh"]["nx"], math.ceil(1.0 * math.sqrt(2) / r["eps_m"] - 1e-9))

    def test_exact_multiple_of_h_is_not_rounded_up(self):
        h = self.ok()["mesh"]["h_m"]
        for k in (1, 10, 333, 12345):
            self.assertEqual(self.ok(domain_x_m=k * h)["mesh"]["nx"], k)

    def test_xi_v_quasi_steady_check(self):
        rat = rho_vs(-20) / RHO_ICE
        st = {c["name"]: c["status"] for c in self.ok(xi_v=0.5 * rat)["checks"]}
        self.assertEqual(st["xi_v_quasi_steady"], "ok")
        st = {c["name"]: c["status"] for c in self.ok(xi_v=2 * rat)["checks"]}
        self.assertEqual(st["xi_v_quasi_steady"], "warn")

    # BUG: preprocess/mesh_params/tool.py:92 takes alpha_c from the Libbrecht law
    # without checking it. For sigma_surf below sigma0/69 the law returns its
    # 1e-30 floor; the tool then returns eps = 6.4e-35 m, Nx = 2.2e30 and every
    # check "ok". Right: raise ToolInputError (or a "fail" check) when the
    # derived alpha_c is at the floor.
    @unittest.expectedFailure
    def test_sigma_surf_at_alpha_floor_is_flagged(self):
        try:
            r = compute_mesh_params(MeshParamsInput(domain_x_m=1e-4, sigma_surf=1e-6))
        except ToolInputError:
            return
        self.assertTrue({"fail", "warn"} & {c["status"] for c in r["checks"]}, r["eps_m"])

    # BUG: preprocess/mesh_params/tool.py:60-61 `finite()` accepts bool because
    # bool is an int, so JSON {"alpha_c": true} runs as alpha_c = 1. Right:
    # reject booleans.
    @unittest.expectedFailure
    def test_json_boolean_is_rejected(self):
        r = cli("--json", json.dumps({"domain_x_m": 1e-4, "alpha_c": True}))
        self.assertEqual(r.returncode, 2)


class CommandLine(unittest.TestCase):
    def test_out_of_range_inputs_exit_2_with_message(self):
        for args, word in ((("--T0=-173",), "temperature_c"), (("--vn", 0), "front_velocity"),
                           (("--alpha_range", 1e-3, 1e-2, 1), "n_points")):
            with self.subTest(args=args):
                base = ["--Lx", 1e-4] + ([] if args[0] == "--alpha_range" else ["--alpha", 1e-2])
                r = cli(*base, *args)
                self.assertEqual(r.returncode, 2)
                self.assertIn(word, r.stderr)
                self.assertEqual(r.stdout, "")

    def test_nan_alpha_rejected(self):
        r = cli("--Lx", 1e-4, "--alpha", "nan")
        self.assertEqual(r.returncode, 2)
        self.assertIn("alpha_c", r.stderr)

    def test_two_alpha_sources_is_a_usage_error(self):
        r = cli("--Lx", 1e-4, "--alpha", 1e-2, "--sigma_surf", 0.01)
        self.assertEqual(r.returncode, 2)

    def test_non_numeric_alpha_range_fails_without_output(self):
        # Fails, but with a Python traceback (ValueError is not caught in cli.py:107).
        r = cli("--Lx", 1e-4, "--alpha_range", "a", "b", "c")
        self.assertNotEqual(r.returncode, 0)
        self.assertEqual(r.stdout, "")

    def test_quiet_output(self):
        r = cli("--Lx", 3e-4, "--Ly", 1e-4, "--alpha", 1e-2, "--quiet")
        lines = r.stdout.split()
        res = compute_mesh_params(MeshParamsInput(domain_x_m=3e-4, domain_y_m=1e-4, alpha_c=1e-2))
        self.assertEqual(lines, [f"{res['eps_m']:.4e}", str(res["mesh"]["nx"]), str(res["mesh"]["ny"])])

    def test_patch_touches_only_eps_and_counts(self):
        text = ("# geometry\n-eps 1.0e-6   # old value\n-eps_valid_temp -20\n-Nx 10\n-Nxx 7\n"
                "# -Ny 99 commented\n-Lx 3e-4\n")
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "g.opts"
            f.write_text(text)
            r = cli("--Lx", 3e-4, "--Ly", 1e-4, "--alpha", 1e-2, "--quiet", "--patch", f)
            self.assertEqual(r.returncode, 0, r.stderr)
            eps, nx, ny = r.stdout.split()
            new = f.read_text().splitlines()
        self.assertIn(f"-eps {eps}   # old value", new)
        self.assertIn("-eps_valid_temp -20", new)
        self.assertIn(f"-Nx {nx}", new)
        self.assertIn("-Nxx 7", new)
        self.assertIn("# -Ny 99 commented", new)
        self.assertEqual(new[-1], f"-Ny {ny}")
        self.assertIn("-Lx 3e-4", new)
        self.assertEqual(len(new), len(text.splitlines()) + 1)


if __name__ == "__main__":
    unittest.main()
