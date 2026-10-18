"""Tests for the mesh_params tool. Run from the project root:  python -m unittest discover -s agentic/tests -v"""
import json
import math
import unittest
from pathlib import Path

from agentic.tools.mesh_params import MeshParamsInput, ToolInputError, compute_mesh_params, solver_options, sweep_alpha_c
from agentic.tools.mesh_params import physics
from agentic.tools.mesh_params.constants import CONSTANTS, C_SOURCE_LINKS

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = json.loads((Path(__file__).parent / "golden" / "comp_eps_golden.json").read_text())


def close(a, b, rel=1e-12):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)


class GoldenAgainstOldScript(unittest.TestCase):
    """The redesigned tool must reproduce the numbers of preprocess/comp_eps.py (frozen 2026-10-17)."""

    def test_cases(self):
        for case in GOLDEN["cases"]:
            c, old = case["input"], case["result"]
            new = compute_mesh_params(MeshParamsInput(
                domain_x_m=c["Lx"], domain_y_m=c.get("Ly", 0.0), domain_z_m=c.get("Lz", 0.0), grain_radius_m=c["Rave"],
                temperature_c=c["T0_C"], alpha_c=c["alpha_c"], front_velocity_m_per_s=c.get("v_n", 1e-9),
                safety=c.get("safety", 0.5), thermal_channel=c.get("thermal_channel", "mean"), xi_v=c.get("xi_v")))
            with self.subTest(case=c):
                self.assertTrue(close(new["eps_m"], old["eps"]))
                self.assertEqual((new["mesh"]["nx"], new["mesh"]["ny"], new["mesh"]["nz"]), (old["Nx"], old["Ny"], old["Nz"]))
                for k_new, k_old in (("tau_sub_s", "tau_sub"), ("lam_sub", "lam_sub"), ("M_sub", "M_sub"),
                                     ("alpha_src_per_s", "alpha_src"), ("min_beta_ratio", "min_beta_ratio"),
                                     ("kinetic_fraction", "kinetic_frac")):
                    self.assertTrue(close(new["derived"][k_new], old[k_old]), k_new)
                for k in ("B-HEAT", "B-VAPOR", "B-KINETIC", "B-CURV"):
                    self.assertTrue(close(new["eps_bounds_m"][k], old[k]), k)
                self.assertEqual(new["eps_bounds_m"]["binding"], old["binding"])
                self.assertTrue(close(new["thermo"]["peclet"], old["Pe"]))

    def test_libbrecht(self):
        for r in GOLDEN["sigma0"]:
            self.assertTrue(close(physics.sigma0(r["T0_C"]), r["sigma0"]))
        for r in GOLDEN["alpha_libbrecht"]:
            self.assertTrue(close(physics.alpha_libbrecht(r["T0_C"], r["sigma_surf"]), r["alpha"]))


class Behaviour(unittest.TestCase):
    base = dict(domain_x_m=3e-4, temperature_c=-20.0, alpha_c=1e-3)

    def test_band_element_counts(self):
        r = compute_mesh_params(MeshParamsInput(**self.base))["resolution"]
        self.assertAlmostEqual(r["elements_across_5_95"], 8.33, places=2)
        self.assertAlmostEqual(r["elements_across_1_99"], 13.0, places=1)

    def test_solver_options_only_eps_and_n(self):
        o = compute_mesh_params(MeshParamsInput(domain_x_m=1e-4, domain_y_m=1e-4, temperature_c=-20, alpha_c=1e-3))["solver_options"]
        self.assertEqual(set(o), {"-eps", "-Nx", "-Ny"})
        self.assertEqual(set(solver_options(1e-6, 10, 0, 0, 1)), {"-eps", "-Nx"})

    def test_json_serialisable_and_deterministic(self):
        a = compute_mesh_params(MeshParamsInput(**self.base))
        self.assertEqual(json.dumps(a, sort_keys=True), json.dumps(compute_mesh_params(MeshParamsInput(**self.base)), sort_keys=True))

    def test_input_errors_listed_together(self):
        with self.assertRaises(ToolInputError) as cm:
            compute_mesh_params(MeshParamsInput(domain_x_m=-1.0, temperature_c=+5.0, alpha_c=2.0, safety=0.0))
        msg = str(cm.exception)
        for part in ("domain_x_m", "temperature_c", "alpha_c", "safety"):
            self.assertIn(part, msg)

    def test_exactly_one_of_alpha_and_sigma(self):
        with self.assertRaises(ToolInputError):
            compute_mesh_params(MeshParamsInput(domain_x_m=1e-4))
        with self.assertRaises(ToolInputError):
            compute_mesh_params(MeshParamsInput(domain_x_m=1e-4, alpha_c=1e-3, sigma_surf=1e-3))

    def test_sweep_log_uniform_and_binding_switch(self):
        rows = sweep_alpha_c(MeshParamsInput(domain_x_m=1e-4, temperature_c=-20.0), 1e-4, 1e-1, 4)
        self.assertEqual(len(rows), 4)
        self.assertTrue(close(rows[0][0], 1e-4) and close(rows[-1][0], 1e-1))
        ratios = [rows[k + 1][0] / rows[k][0] for k in range(3)]
        self.assertTrue(all(close(r, ratios[0], 1e-9) for r in ratios))   # log-uniform
        # eps is not monotone in alpha_c: the binding bound switches from the kinetic one (cold, small alpha_c)
        # to the heat one (large alpha_c), which is why a sweep before trusting a mesh is worth running
        self.assertEqual(rows[0][1]["eps_bounds_m"]["binding"], "B-KINETIC")
        self.assertEqual(rows[-1][1]["eps_bounds_m"]["binding"], "B-HEAT")

    def test_two_sided_caveat_is_always_reported(self):
        names = {c["name"] for c in compute_mesh_params(MeshParamsInput(**self.base))["checks"]}
        self.assertIn("thin_interface_two_sided_constants", names)
        self.assertIn("vapor_diffusivity_is_air", names)


class ConstantsMatchSolverSource(unittest.TestCase):
    """Detect drift between this tool's constants and the C solver."""

    def test_links(self):
        for name, value, rel_path, needle in C_SOURCE_LINKS:
            f = ROOT / rel_path
            if not f.exists():
                self.skipTest(f"{rel_path} not present (tool used outside the solver tree)")
            with self.subTest(constant=name):
                self.assertEqual(getattr(CONSTANTS, name), value)
                self.assertIn(needle, f.read_text(), f"{rel_path} no longer contains '{needle}'; update the tool or the table")


if __name__ == "__main__":
    unittest.main()
