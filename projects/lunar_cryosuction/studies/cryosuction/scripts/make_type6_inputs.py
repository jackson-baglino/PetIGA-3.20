"""Write the geometry, experiment and tests files for plan type 6 (ice moving
through a constriction).

Plan (SimulationPlan_v2, type 6):
  pore    dog-bone: wide sections W = 150 um, straight neck wc = 40, 80 or
          150 um (150 = no constriction), neck 300 um long, total 1.2 mm in
          the plan (1.1 mm here, see below), symmetric about the pore axis (preprocess/build_geometry_dogbone.py)
  ice     lens 150 um long in the cold wide section
  vapor   local-saturation pair: each end at rho_vs of its own temperature
  thermal gradient, warm top. Depth is drawn downward (plan tool T7), so the
          top is x = 0: T(0) is the warm end, grad_temp0 = -G along x
  runs    6A-1..3  wc 40, 80, 150 um, G 100 K/m, theta 60, 60 days
          6A-4..6  wc 40, 80, 150 um, G 300 K/m, theta 60, 30 days
          6B-1     wc 40 um, G 300 K/m, theta 90, 30 days (no capillary force
                   on the flat neck walls)
          6B-2     wc 80 um, otherwise as 6B-1 (user, 2026-11-12)
          The first 6B-1 (2026-11-09) started from faces 0.375 um out of flat
          (two crossed 7.5 mm discs); the user voided it (2026-11-12). Theta-90
          geometries are now 'lens90flat' (exactly flat faces, -ice_flat_*);
          --only-6b-flat writes the 6B-1/6B-2 reruns with 401 snapshots
          (type6b_flat_tests.txt).

Choices made here, not in the plan (printed, and in the file headers):
  shoulders  C-infinity, 100 um long each (steepest wall slope 1.10 at
             wc = 40 um); the straight neck is the full 300 um
  length     1.1 mm: each wide section is flat for 300 um, the lens plus 75 um
             on either side of it; the warm side gets the same 300 um so the
             pore is symmetric (user, 2026-11-09)
  lens       150 um long on the axis, centered in the cold flat section
             (apexes x = 875-1025 um); faces start at the run's wall-BC angle
             (user, 2026-11-09), so the theta 60 runs start with concave faces
             (-bridge_concaveL/R, FREEZE 18270f6f) and 6B with flat ones.
             One geometry file per (neck, lens angle); they share the mesh.
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

LX, W = 1.1e-3, 1.5e-4
NECKS = (4.0e-5, 8.0e-5, 1.5e-4)
TODAY = datetime.date.today().isoformat()


def geom_name(wc, theta, w=W, seed=0.0):
    lens = "90flat" if theta == 90 else f"{theta:g}"
    sd = f"_seed{seed*1e6:g}um" if seed else ""
    return f"dogbone_2D_L{LX*1e6:.0f}um_W{w*1e6:.0f}um_wc{wc*1e6:.0f}um_lens{lens}{sd}_eps0.86um"


def build_geometry(wc, theta, dtmax, dry, keep_mesh=False, w=W, seed=0.0):
    """keep_mesh: the mesh file exists and is used by runs already made; build
    into a temporary file and require it to be byte-identical, so the shared
    mesh is never silently replaced."""
    name = geom_name(wc, theta, w, seed)
    final = MESH / (f"dogbone_wc{wc*1e6:.0f}um.dat" if w == W else f"dogbone_W{w*1e6:.0f}um_wc{wc*1e6:.0f}um.dat")
    keep_mesh = keep_mesh and final.exists()
    dat = final.with_suffix(".new.dat") if keep_mesh else final
    opts = GEO / f"{name}.opts"
    png = ROOT / f"preprocess/{name}.pdf"
    head = (f"Plan type 6 (make_type6_inputs.py, {TODAY}): dog-bone pore, neck wc = {wc*1e6:.0f} um"
            + (" (= W: no constriction)" if wc >= w else "") + (f", W = {w*1e6:.0f} um" if w != W else "") + f", lens faces at {theta:g} deg.")
    cmd = [sys.executable, str(ROOT / "preprocess/build_geometry_dogbone.py"), "--wc", f"{wc:.6e}",
           "--Lx", f"{LX:.6e}", "--W", f"{w:.6e}", "--eps", f"{EPS:.4e}", "--lens-theta", f"{theta:g}",
           "--out", str(dat), "--opts", str(opts), "--png", str(png), "--header", head]
    print("  " + " ".join(c.replace(str(ROOT) + "/", "") for c in cmd[1:-2]) + " --header ...")
    if dry:
        return name
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    print("    " + r.stdout.replace("\n", "\n    ").rstrip())
    if r.returncode != 0:
        raise SystemExit(f"dog-bone builder failed for wc {wc}:\n{r.stderr[-1500:]}")
    if keep_mesh:
        same = dat.read_bytes() == final.read_bytes()
        dat.unlink()
        if not same:
            raise SystemExit(f"rebuilt mesh differs from {final.relative_to(ROOT)}; not replacing it")
        print(f"    mesh identical to {final.relative_to(ROOT)} (kept)")
        opts.write_text(opts.read_text().replace(dat.name, final.name))
    extra = f"-dtmax {dtmax:.4e}\n"
    if seed:
        extra += (f"-ice_flat_seed {seed:.4e}                 # warm (low-x) face x_f - R_f + A cos(2 pi y/Ly): "
                  f"axis {2*seed*1e6:g} um ahead of the walls (solver fe8a7334)\n")
    opts.write_text(opts.read_text() + extra)
    return name


def experiment(name, theta, G, days, per_day, purpose, dry, n_out=None, event=False):
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
-n_out {n_out or int(days * per_day) + 1}                              # {f'{(n_out - 1) / days:.2f}' if n_out else f'{per_day:g}'} snapshot(s) a day
-beta_sub0 {BETA0:.4e}                     # K&P beta0 at -20 C, alpha_c = 1e-2
-gamma_ia {G_IA:g}
-gamma_as {G_AS:g}
-gamma_is {gis:.7f}                     # cos(theta) = {__import__('math').cos(__import__('math').radians(theta)):+.6f}
-dtCFL 1
-outp 0
""" + ("""-outp_dphi 0.5                          # extra snapshot when phi changed by 0.5 anywhere (e5a593e7)
-outp_dphi_min_dt 3600                  # at most one such snapshot an hour
""" if event else "")
    print(f"  {'(dry) ' if dry else ''}wrote {(EXP / f'{name}.opts').relative_to(ROOT)}  "
          f"(rhovfix_lo {lo:.10f}, hi {hi:.10f})")
    if not dry:
        (EXP / f"{name}.opts").write_text(text)
    return name


G_SET = ((10.0, 360), (20.0, 180), (50.0, 75))   # RunPlan_2026-11-12 6C: G [K/m], days


def plan_1112(dtmax, dry):
    """6C: wc 40/80/150, theta 60, G 10/20/50. 6D: theta-90 face stability in the straight pore."""
    print("\n6C geometries (existing lens60, shared meshes)")
    g60 = {wc: build_geometry(wc, 60, dtmax, dry, keep_mesh=True) for wc in NECKS}
    print("\n6D geometries")
    seed = 0.25e-6                                    # axis leads the walls by 0.5 um
    d150s = build_geometry(W, 90, dtmax, dry, keep_mesh=True, seed=seed)
    d150 = build_geometry(W, 90, dtmax, dry, keep_mesh=True)
    d75s = build_geometry(75e-6, 90, dtmax, dry, w=75e-6, seed=seed)
    d150_75 = build_geometry(W, 75, dtmax, dry, keep_mesh=True)
    print("\nExperiments")
    ex = {}
    for G, days in G_SET:
        ex[(60, G)] = experiment(f"t6ls_T-20_G{G:g}warmtop_theta60_{days}d_n401_ac1e-2", 60, G, days, None,
                                 f"Plan type 6C (RunPlan_2026-11-12): lens through a neck, G = {G:g} K/m.",
                                 dry, n_out=401, event=True)
    for G, days in ((10.0, 120), (20.0, 120), (50.0, 60)):
        ex[(90, G)] = experiment(f"t6ls_T-20_G{G:g}warmtop_theta90_{days}d_n401_ac1e-2", 90, G, days, None,
                                 f"Plan type 6D (RunPlan_2026-11-12): theta-90 face stability, G = {G:g} K/m.",
                                 dry, n_out=401, event=True)
    ex[(75, 50.0)] = experiment("t6ls_T-20_G50warmtop_theta75_60d_n401_ac1e-2", 75, 50.0, 60, None,
                                "Plan type 6D-6 (RunPlan_2026-11-12): theta 75, straight pore, G = 50 K/m.",
                                dry, n_out=401, event=True)
    c = [f"# Plan type 6C (make_type6_inputs.py --plan-1112, {TODAY}): theta-60 lens through necks of 40, 80 "
         f"and 150 um at G 10, 20, 50 K/m. Submit --into the type-6 batch."]
    k = 1
    for G, _ in G_SET:
        for wc in NECKS:
            c.append(f"{g60[wc]}:{ex[(60, G)]}   # 6C-{k}")
            k += 1
    d = [f"# Plan type 6D (make_type6_inputs.py --plan-1112, {TODAY}): theta-90 face stability in the straight "
         f"1100 um pore; flat faces, warm face seeded {seed*1e6:g} um (axis 0.5 um ahead). Submit --into the type-6 batch.",
         f"{d150s}:{ex[(90, 10.0)]}   # 6D-1",
         f"{d150s}:{ex[(90, 20.0)]}   # 6D-2",
         f"{d150s}:{ex[(90, 50.0)]}   # 6D-3",
         f"{d75s}:{ex[(90, 50.0)]}   # 6D-4 (W 75 um)",
         f"{d150}:{ex[(90, 50.0)]}   # 6D-5 (no seed)",
         f"{d150_75}:{ex[(75, 50.0)]}   # 6D-6 (theta 75)"]
    for name, lines in (("type6c_tests.txt", c), ("type6d_tests.txt", d)):
        path = STUDY / name
        print(f"\n  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}")
        print("\n".join("    " + ln for ln in lines))
        if not dry:
            path.write_text("\n".join(lines) + "\n")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print what would be written")
    ap.add_argument("--plan-1112", action="store_true",
                    help="write 6C (G 10/20/50, theta 60) and 6D (face stability) of RunPlan_2026-11-12, "
                         "type6c_tests.txt and type6d_tests.txt")
    ap.add_argument("--only-6b-flat", action="store_true",
                    help="write only the flat-face 6B runs (G 20, 50; 401 snapshots) and type6b_flat_tests.txt")
    a = ap.parse_args()
    dtmax = 2.0 * tau_sub(EPS)
    print(f"eps {EPS:.4e} m, dtmax = 2 tau_sub = {dtmax:.4e} s")
    if a.only_6b_flat:
        # User 2026-11-13: no new G 300 runs; 6B at the plan's gradients, durations as 6C.
        print("\n6B: theta 90, flat faces, G 20 and 50")
        g = {wc: build_geometry(wc, 90, dtmax, a.dry_run, keep_mesh=True) for wc in NECKS[:2]}
        e = {G: experiment(f"t6ls_T-20_G{G:g}warmtop_theta90_{days}d_n401_ac1e-2", 90, G, days, None,
                           f"Plan type 6B (2026-11-13): theta 90 lens, exactly flat starting faces, through a "
                           f"neck, G = {G:g} K/m.", a.dry_run, n_out=401, event=True)
             for G, days in ((20.0, 180), (50.0, 75))}
        path = STUDY / "type6b_flat_tests.txt"
        text = (f"# Plan type 6B (make_type6_inputs.py --only-6b-flat, {TODAY}): theta 90 lens with exactly flat "
                f"faces through the 40 and 80 um necks at G 20 and 50 K/m (pairs with 6C-4/5 and 6C-7/8).\n"
                f"# Replaces the voided 6B-1 (G 300, crossed-disc faces). Submit --into the type-6 batch.\n"
                f"{g[NECKS[0]]}:{e[20.0]}   # 6B-1\n{g[NECKS[1]]}:{e[20.0]}   # 6B-2\n"
                f"{g[NECKS[0]]}:{e[50.0]}   # 6B-3\n{g[NECKS[1]]}:{e[50.0]}   # 6B-4\n")
        print(f"\n  {'(dry) ' if a.dry_run else ''}wrote {path.relative_to(ROOT)}\n{text}")
        if not a.dry_run:
            path.write_text(text)
        return 0
    if a.plan_1112:
        return plan_1112(dtmax, a.dry_run)
    print("\nType 6 geometry")
    geo = {wc: build_geometry(wc, 60, dtmax, a.dry_run) for wc in NECKS}
    geo90 = build_geometry(NECKS[0], 90, dtmax, a.dry_run)
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
    lines += [f"{geo90}:{e90}   # 6B-1"]
    path = STUDY / "type6_tests.txt"
    print(f"\n  {'(dry) ' if a.dry_run else ''}wrote {path.relative_to(ROOT)}")
    if not a.dry_run:
        path.write_text("\n".join(lines) + "\n")
    print("\n".join("    " + ln for ln in lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
