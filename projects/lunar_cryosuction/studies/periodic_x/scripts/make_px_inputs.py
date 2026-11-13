#!/usr/bin/env python3
"""make_px_inputs.py: layouts, meshes and the experiment file for the
-periodic_x verification (studies/periodic_x/verification/verify_periodic_x.py).

Five small layouts at eps 1.5 um (h = eps/sqrt(2) = 1.06 um; 114 x 40 nodes,
13.7 k DOF each), written to preprocess/layouts/px_*.json and built with
preprocess/build_ice_layout.py into inputs/geometry/periodic_x/:

  px_flat120x40_lens_mid       flat periodic channel, 40 um gap, Lx 120 um,
                               20 um lens (on the axis) centered at Lx/2
  px_flat120x40_lens_seam      the same lens shifted by Lx/2: across the seam
  px_open_flat120x40_lens_mid  the mid lens in a sealed, non-periodic channel
                               (reference for drift and conservation)
  px_cd120_W40_wt16_lens       converging-diverging cell (sinhyp, W 40 um,
                               throat 16 um at Lx/2), lens 10-25 um left of it
  px_cd120_W40_wt16_shift_lens the cd cell and its lens shifted by -Lx/2
                               (throat on the seam, lens across it)

The experiment px_T-20_iso_theta60_2d_ac1e-2: -20 C, T pinned on the two walls
(isothermal; the x faces are a seam), no vapor reservoir, wall angle 60 deg,
2 days, a snapshot every 6 h. dtmax = 2 tau_sub(1.5 um) = 4.5257e3 s
(studies/cryosuction/scripts/make_type34_inputs.tau_sub).

    venv_cryosuction/bin/python studies/periodic_x/scripts/make_px_inputs.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "studies/cryosuction/scripts"))
from make_type34_inputs import BETA0, G_AS, G_IA, gamma_is, tau_sub  # noqa: E402

EPS = 1.5e-6
LX = 1.2e-4
DTMAX = 2 * tau_sub(EPS)
BASE = {"eps": EPS, "Lx": LX, "theta_ice_init": 60.0, "dtmax": DTMAX, "caps": [], "grains": []}
CD = {"profile": "sinhyp", "W": 4.0e-5, "wt": 1.6e-5}
FLAT = {"flat": True}

LAYOUTS = {
    "px_flat120x40_lens_mid_eps1.5um": dict(
        BASE, periodic=True, mean_gap=4.0e-5, walls=FLAT, plugs=[{"x0": 5.0e-5, "x1": 7.0e-5}],
        comment="periodic_x verification b/c1: flat periodic channel, 20 um lens centered at Lx/2"),
    "px_flat120x40_lens_seam_eps1.5um": dict(
        BASE, periodic=True, mean_gap=4.0e-5, walls=FLAT, plugs=[{"x0": -1.0e-5, "x1": 1.0e-5}],
        comment="periodic_x verification c1: the mid lens shifted by Lx/2, across the seam"),
    "px_open_flat120x40_lens_mid_eps1.5um": dict(
        BASE, periodic=False, mean_gap=4.0e-5, walls=FLAT, plugs=[{"x0": 5.0e-5, "x1": 7.0e-5}],
        comment="periodic_x verification b reference: the mid lens in a sealed, non-periodic channel"),
    "px_cd120_W40_wt16_lens_eps1.5um": dict(
        BASE, periodic=True, walls=CD, plugs=[{"x0": 3.5e-5, "x1": 5.0e-5}],
        comment="periodic_x verification c2/d: converging-diverging cell, throat at Lx/2, lens left of it"),
    "px_cd120_W40_wt16_shift_lens_eps1.5um": dict(
        BASE, periodic=True, walls=dict(CD, x_throat=0.0), plugs=[{"x0": -2.5e-5, "x1": -1.0e-5}],
        comment="periodic_x verification c2: the cd cell and lens shifted by -Lx/2 (throat on the seam)"),
}

EXP_NAME = "px_T-20_iso_theta60_2d_ac1e-2"
T_FINAL = 2 * 86400.0
N_OUT = 9


def experiment_text() -> str:
    return f"""# experiment/periodic_x/{EXP_NAME}.opts   (studies/periodic_x/scripts/make_px_inputs.py)
#
# -periodic_x verification: isothermal at -20 C (T pinned on the two walls; the
# x faces are a periodic seam), no vapor reservoir (water conserved), wall angle
# 60 deg from Young (gamma_is = 0.3 - 0.109 cos theta), 2 days, a snapshot every 6 h.
# -dtmax is in the GEOMETRY file (2 tau_sub at its eps).
#
-temp -20.0
-humidity 1.0
-grad_temp0 0.0,0.0,0.0
-flag_BC_Tfix 1
-flag_BC_rhovfix 0                      # no vapor in or out: water conserved
-phi_open_bc 0
-t_final {T_FINAL:.6e}
-n_out {N_OUT}
-outp 0
-beta_sub0 {BETA0:.4e}                     # K&P beta0 at -20 C, alpha_c = 1e-2
-gamma_ia {G_IA}
-gamma_as {G_AS}
-gamma_is {gamma_is(60.0):.7f}                     # cos(theta) = +0.500000
-dtCFL 1
"""


def main() -> int:
    geo = ROOT / "inputs/geometry/periodic_x"
    exp = ROOT / "inputs/experiment/periodic_x"
    geo.mkdir(parents=True, exist_ok=True)
    exp.mkdir(parents=True, exist_ok=True)
    print(f"eps {EPS:.4e} m, dtmax = 2 tau_sub = {DTMAX:.4e} s")
    for name, spec in LAYOUTS.items():
        lay = ROOT / "preprocess/layouts" / f"{name}.json"
        lay.write_text(json.dumps(dict({"name": name}, **spec), indent=2) + "\n")
        subprocess.run([sys.executable, str(ROOT / "preprocess/build_ice_layout.py"), str(lay),
                        "--outdir", str(geo), "--no-pdf"], check=True, cwd=ROOT)
    (exp / f"{EXP_NAME}.opts").write_text(experiment_text())
    print(f"wrote {exp / (EXP_NAME + '.opts')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
