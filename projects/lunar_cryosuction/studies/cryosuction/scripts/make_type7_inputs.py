"""Write the experiments and tests file for the type-7 v5 set (user, 2026-11-13): 13 runs.

  7A-1..7  every configuration C1..C7, isothermal (-20 C on every face), collector at x = 0
  7C-1,4,7 C1, C4, C7 at G = 20 K/m, warm end at the collector (x = 0)
  7D-1,4,7 C1, C4, C7 at G = 20 K/m, warm end at the dead end (x = L)

Geometries: preprocess/layouts/type7v5_C*.json built with preprocess/build_ice_layout.py
(inputs/geometry/cryosuction/type7v5_C*_2D_L1800um_eps0.86um.opts).
Collector: rho_v = 0.990 rho_vs(T0) on x0 (decision D2; the solver pins rho_v to frac * rho_vs(temp0),
STATUS issue 5); dead end x = L sealed (no vapor flux); phi: 90-degree face on both ends. Walls y0, y1:
contact angle 60 deg. Isothermal set: T fixed on every face (as the earlier 7A). Gradient sets: T fixed on
the two ends only, walls insulating (as types 1, 2, 6). Output: 401 snapshots plus event snapshots.

    venv_cryosuction/bin/python studies/cryosuction/scripts/make_type7_inputs.py [--dry-run]
"""
from __future__ import annotations

import argparse
import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "inputs/experiment/cryosuction"
STUDY = ROOT / "studies/cryosuction"
TODAY = datetime.date.today().isoformat()
DAYS, N_OUT, G = 60, 401, 20.0
GEOM = "type7v5_C{}_2D_L1800um_eps0.86um"

COMMON = """-temp -20.0
-humidity 1.0
-flag_BC_rhovfix 1
-rhovfix_faces x0
-rhovfix_lo 0.990                       # collector: 0.990 rho_vs(T0) (decision D2)
-phi_open_bc 0                          # 90-degree phi face on the pore ends
-t_final {tf:.6e}                       # {days} days
-n_out {nout}                              # time-uniform snapshots
-beta_sub0 7.9408e5                     # K&P beta0 at -20 C, alpha_c = 1e-2
-gamma_ia 0.109
-gamma_as 0.3
-gamma_is 0.2455000                     # Young: cos(theta) = 0.5, theta = 60 deg
-dtCFL 1
-outp 0
-outp_dphi 0.5                          # extra snapshot when phi changed by 0.5 anywhere (e5a593e7)
-outp_dphi_min_dt 3600                  # at most one such snapshot an hour
"""


def experiment(name, head, thermal, dry):
    text = (f"# experiment/cryosuction/{name}.opts   (make_type7_inputs.py, {TODAY})\n#\n" + head +
            "# -dtmax is in the GEOMETRY file (2 tau_sub at its eps).\n#\n" + thermal +
            COMMON.format(tf=DAYS * 86400.0, days=DAYS, nout=N_OUT))
    path = EXP / f"{name}.opts"
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}")
    if not dry:
        path.write_text(text)
    return name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    iso = experiment(
        f"t7_T-20_h0.990_iso_theta60_{DAYS}d_n401_ac1e-2",
        "# Plan type 7A (v5): isothermal harvest; collector at x = 0, sealed dead end at x = L.\n",
        "-grad_temp0 0.0,0.0,0.0\n-flag_BC_Tfix 1                          # T fixed on every face\n", a.dry_run)
    warm0 = experiment(
        f"t7_T-20_h0.990_G{G:g}warmcollector_theta60_{DAYS}d_n401_ac1e-2",
        f"# Plan type 7C (v5): G = {G:g} K/m, warm end at the collector (x = 0); walls insulating.\n",
        f"-grad_temp0 {-G:.1f},0.0,0.0\n-flag_BC_Tfix 1\n-Tfix_faces x0,x1\n", a.dry_run)
    warmL = experiment(
        f"t7_T-20_h0.990_G{G:g}warmdeadend_theta60_{DAYS}d_n401_ac1e-2",
        f"# Plan type 7D (v5): G = {G:g} K/m, warm end at the dead end (x = L); walls insulating.\n",
        f"-grad_temp0 {G:.1f},0.0,0.0\n-flag_BC_Tfix 1\n-Tfix_faces x0,x1\n", a.dry_run)
    lines = [f"# Plan type 7 v5 (make_type7_inputs.py, {TODAY}): one smooth rough pore (L 1.8 mm, throats "
             f"50/65/45/60 um), seven ice configurations C1 (40 %) .. C7 (4 %). 13 runs, 2.39 M DoF each."]
    lines += [f"{GEOM.format(i)}:{iso}   # 7A-{i}" for i in range(1, 8)]
    lines += [f"{GEOM.format(i)}:{warm0}   # 7C-{i}" for i in (1, 4, 7)]
    lines += [f"{GEOM.format(i)}:{warmL}   # 7D-{i}" for i in (1, 4, 7)]
    path = STUDY / "type7v5_tests.txt"
    print(f"  {'(dry) ' if a.dry_run else ''}wrote {path.relative_to(ROOT)}")
    if not a.dry_run:
        path.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
