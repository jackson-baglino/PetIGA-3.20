"""Write the geometry, experiment and tests files for plan type 6 (ice moving
through a constriction).

Plan (SimulationPlan_v2, type 6):
  pore    dog-bone: wide sections W = 150 um, straight neck wc = 40, 80 or
          150 um (150 = no constriction), neck 300 um long, total 1.2 mm,
          symmetric about the pore axis (preprocess/build_geometry_dogbone.py)
  ice     lens 150 um long in the cold wide section
  vapor   local-saturation pair: each end at rho_vs of its own temperature
  thermal gradient, warm top. Depth is drawn downward (plan tool T7), so the
          top is x = 0: T(0) is the warm end, grad_temp0 = -G along x
  runs    6A-1..3  wc 40, 80, 150 um, G 100 K/m, theta 60, 60 days
          6A-4..6  wc 40, 80, 150 um, G 300 K/m, theta 60, 30 days
          6B-1     wc 40 um, G 300 K/m, theta 90, 30 days (no capillary force
                   on the flat neck walls)

Choices made here, not in the plan (printed, and in the file headers):
  shoulders  C-infinity, 100 um long each (steepest wall slope 1.10 at
             wc = 40 um); the straight neck is the full 300 um, so the wide
             sections are flat for 350 um at each end
  lens       150 um long, centered in the cold flat section (x = 950-1100 um),
             flat faces at t = 0 (as in types 3 and 4)
  eps        8.5840e-7 m, the production value (issue 1); comp_eps at the
             G = 300 speed (5e-10 m/s, Rave 20 um) allows up to 1.28e-6 m
  output     one snapshot a day at G = 100, two a day at G = 300

Local-saturation fractions come from rhov_bc.fracs (plan tool T4): the solver
pins rho_v to frac * rho_vs(temp0) (src/cryosuction_main.c, IGASetBoundaryValue
in the vapor BC block), so frac = rho_vs(T_end) / rho_vs(temp0).

    venv_cryosuction/bin/python studies/cryosuction/scripts/make_type6_inputs.py [--dry-run]
"""
from __future__ import annotations

import argparse
import datetime
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT / "studies/cryosuction/scripts"))
from make_type34_inputs import BETA0, DAY, EPS, EXP, G_AS, G_IA, GEO, MESH, STUDY, T0, gamma_is, tau_sub  # noqa: E402
from rhov_bc import fracs  # noqa: E402

LX, W = 1.2e-3, 1.5e-4
NECKS = (4.0e-5, 8.0e-5, 1.5e-4)
TODAY = datetime.date.today().isoformat()


def geom_name(wc):
    return f"dogbone_2D_L1200um_W150um_wc{wc*1e6:.0f}um_eps0.86um"


def build_geometry(wc, dtmax, dry):
    name = geom_name(wc)
    dat = MESH / f"dogbone_wc{wc*1e6:.0f}um.dat"
    opts = GEO / f"{name}.opts"
    png = ROOT / f"preprocess/dogbone_wc{wc*1e6:.0f}um.png"
    head = (f"Plan type 6 (make_type6_inputs.py, {TODAY}): dog-bone pore, neck wc = {wc*1e6:.0f} um"
            + (" (= W: no constriction)" if wc >= W else "") + ".")
    cmd = [sys.executable, str(ROOT / "preprocess/build_geometry_dogbone.py"), "--wc", f"{wc:.6e}",
           "--Lx", f"{LX:.6e}", "--W", f"{W:.6e}", "--eps", f"{EPS:.4e}",
           "--out", str(dat), "--opts", str(opts), "--png", str(png), "--header", head]
    print("  " + " ".join(c.replace(str(ROOT) + "/", "") for c in cmd[1:-2]) + " --header ...")
    if dry:
        return name
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    print("    " + r.stdout.replace("\n", "\n    ").rstrip())
    if r.returncode != 0:
        raise SystemExit(f"dog-bone builder failed for wc {wc}:\n{r.stderr[-1500:]}")
    opts.write_text(opts.read_text() + f"-dtmax {dtmax:.4e}\n")
    return name


def experiment(name, theta, G, days, per_day, purpose, dry):
    lo, hi, Tlo, Thi, _ = fracs(T0, -G, LX, "local-sat")
    gis = gamma_is(theta)
    text = f"""# experiment/cryosuction/{name}.opts   (make_type6_inputs.py, {TODAY})
#
# {purpose}
# Gradient G = {G:g} K/m along the pore, WARM at x = 0 (the top; depth is drawn
# downward): grad_temp0 = {-G:g} K/m, so T(0) = {Tlo:.6f} C and T(Lx) = {Thi:.6f} C
# for Lx = {LX*1e3:g} mm. T is fixed on both ends; the walls are insulating.
# Vapor: local-saturation pair, each end at rho_vs of its own temperature
# (rhov_bc.py --T0 {T0:g} --G {-G:g} --Lx {LX:g} --local-sat). The values depend on Lx:
# use this file only with the {LX*1e3:g} mm type-6 pores.
# theta = {theta:g} deg from Young (gamma_is = 0.3 - 0.109 cos theta), {days:g} days.
# -dtmax is in the GEOMETRY file (2 tau_sub depends on its eps).
#
-temp {T0:.1f}
-humidity 1.0
-grad_temp0 {-G:.1f},0.0,0.0
-flag_BC_Tfix 1
-Tfix_faces x0,x1
-flag_BC_rhovfix 1
-rhovfix_faces x0,x1
-rhovfix_lo {lo:.10f}                 # rho_vs({Tlo:.4f} C) / rho_vs({T0:g} C), warm end
-rhovfix_hi {hi:.10f}                 # rho_vs({Thi:.4f} C) / rho_vs({T0:g} C), cold end
-phi_open_bc 0                          # 90-degree phi face on both pore ends
-t_final {days*DAY:.6e}
-n_out {int(days * per_day) + 1}                              # {per_day:g} snapshot(s) a day
-beta_sub0 {BETA0:.4e}                     # K&P beta0 at -20 C, alpha_c = 1e-2
-gamma_ia {G_IA:g}
-gamma_as {G_AS:g}
-gamma_is {gis:.7f}                     # cos(theta) = {__import__('math').cos(__import__('math').radians(theta)):+.6f}
-dtCFL 1
-outp 0
"""
    print(f"  {'(dry) ' if dry else ''}wrote {(EXP / f'{name}.opts').relative_to(ROOT)}  "
          f"(rhovfix_lo {lo:.10f}, hi {hi:.10f})")
    if not dry:
        (EXP / f"{name}.opts").write_text(text)
    return name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print what would be written")
    a = ap.parse_args()
    dtmax = 2.0 * tau_sub(EPS)
    print(f"eps {EPS:.4e} m, dtmax = 2 tau_sub = {dtmax:.4e} s")
    print("\nType 6 geometry")
    geo = {wc: build_geometry(wc, dtmax, a.dry_run) for wc in NECKS}
    print("\nExperiments")
    e100 = experiment("t6ls_T-20_G100warmtop_theta60_60d_ac1e-2", 60, 100.0, 60, 1,
                      "Plan type 6A-1..3: lens crossing a constriction, G = 100 K/m.", a.dry_run)
    e300 = experiment("t6ls_T-20_G300warmtop_theta60_30d_ac1e-2", 60, 300.0, 30, 2,
                      "Plan type 6A-4..6: lens crossing a constriction, G = 300 K/m.", a.dry_run)
    e90 = experiment("t6ls_T-20_G300warmtop_theta90_30d_ac1e-2", 90, 300.0, 30, 2,
                     "Plan type 6B-1: theta = 90 control (flat faces on the straight neck walls).", a.dry_run)
    lines = [f"# Plan type 6 (make_type6_inputs.py, {TODAY}). RUNS.md IDs in the comments."]
    lines += [f"{geo[wc]}:{e100}   # 6A-{i}" for i, wc in enumerate(NECKS, 1)]
    lines += [f"{geo[wc]}:{e300}   # 6A-{i}" for i, wc in enumerate(NECKS, 4)]
    lines += [f"{geo[NECKS[0]]}:{e90}   # 6B-1"]
    path = STUDY / "type6_tests.txt"
    print(f"\n  {'(dry) ' if a.dry_run else ''}wrote {path.relative_to(ROOT)}")
    if not a.dry_run:
        path.write_text("\n".join(lines) + "\n")
    print("\n".join("    " + ln for ln in lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
