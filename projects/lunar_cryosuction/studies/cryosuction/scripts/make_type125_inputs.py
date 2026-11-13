"""Write the geometry, experiment and tests files for plan types 1, 2 and 5
(SimulationPlan_v2, sections on types 1, 2 and 5).

Type 1  one ice cap on a flat regolith wall (y0), box 600 x 300 um; the far
        face y1 is not a wall (90-degree phi face, sealed, insulating).
Type 2  several caps on the same wall, box 600 or 900 x 300 um.
Type 5  a lens at a smooth throat in a straight 800 x 150 um pore (walls y0, y1).

Initial angles (user, 2026-11-09): droplets and caps start at liquid water's
angle on regolith, THETA_WATER = 40 deg for now; lenses start at the wall-BC
angle THETA_BC = 60 deg (decision D1), which every experiment here prescribes.

Caps (types 1, 2) use the solver's disc IC (-ic_type multi_grains): a disc of
arc radius R whose center sits R cos(theta) below the wall meets it at theta;
the solver lays the IC out on the rectangle, so no mesh or field file is
needed. Type 5 uses preprocess/build_ice_layout.py (mesh + -ic_type phi_file):
the lens follows the equal-curvature rule on the curved throat walls.

Orientation (as types 6, 7): depth is drawn downward, x = 0 is the top, a
"warm top" gradient G is grad_temp0 = -G along x; "warm third" is x < Lx/3.

Vapor sets (studies/cryosuction/scripts/rhov_bc.py, plan tool T4):
  closed     no vapor in or out (-flag_BC_rhovfix 0): water conserved
  ls         local-saturation pair: each end at rho_vs of its own temperature
  res<h>     reservoir pair: both ends at h rho_vs(T0), h = 1 + d0 chi (form ii)
  dh         vapor-gradient twin (isothermal): ends at 1 +/- dh/2, dh equal to
             the local-saturation contrast of the matching gradient run

eps: type 1 8.5840e-7 m (production value, plan: R_min/eps = 29); types 2 and 5
8.0e-7 m (plan: R_min/25 for the 20 um cap; half-throat 20 um / 25).
dtmax = 2 tau_sub(eps) (mesh_params.derived_pf_params via make_type34_inputs).

Gradients and days for types 1 and 2 follow the run plan of 2026-11-12
(docs/run_plan/build_runplan2.py): G_PROD, G_HIGH, G_ONE and the TYPE1 table
below; every type 1/2 experiment writes N_OUT = 401 time-uniform snapshots.

    venv_cryosuction/bin/python studies/cryosuction/scripts/make_type125_inputs.py [--dry-run] [--only 1,2,5]
"""
from __future__ import annotations

import argparse
import datetime
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_type34_inputs import BETA0, DAY, EXP, G_AS, G_IA, GEO, STUDY, T0, gamma_is, tau_sub  # noqa: E402
from rhov_bc import fracs  # noqa: E402

THETA_WATER, THETA_BC = 40.0, 60.0
EPS1, EPS25 = 8.5840e-7, 8.0e-7
LAYOUTS = ROOT / "preprocess/layouts"
TODAY = datetime.date.today().isoformat()
um = 1e-6


def write(path: Path, text: str, dry: bool):
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}")
    if not dry:
        path.write_text(text)


# ------------------------------------------------------------- caps (1, 2)
def cap_geometry(name, lx, ly, caps, eps, purpose, dry):
    """Box lx x ly, wall y0; caps = [(x_center, R_arc)] at THETA_WATER."""
    ct = math.cos(math.radians(THETA_WATER))
    nx, ny = math.ceil(lx * math.sqrt(2) / eps), math.ceil(ly * math.sqrt(2) / eps)
    half = [R * math.sin(math.radians(THETA_WATER)) for _, R in caps]
    for (x, _), b in zip(caps, half):
        assert b < x < lx - b, f"{name}: cap at {x} runs off the wall"
    for k in range(len(caps) - 1):
        gap = (caps[k + 1][0] - half[k + 1]) - (caps[k][0] + half[k])
        assert gap > 6 * eps, f"{name}: caps {k},{k+1} closer than 6 eps ({gap:.2e} m)"
    lines = "".join(f"#   cap at x = {x/um:.1f} um, arc radius {R/um:.1f} um, base {2*b/um:.1f} um, "
                    f"height {R*(1-ct)/um:.2f} um, area {R*R*(math.radians(THETA_WATER)-math.sin(math.radians(THETA_WATER))*ct):.4e} m^2\n"
                    for (x, R), b in zip(caps, half))
    dtmax = 2 * tau_sub(eps)
    text = f"""# =============================================================================
# geometry/cryosuction/{name}.opts   (make_type125_inputs.py, {TODAY})
#
# {purpose}
# Box {lx/um:.0f} x {ly/um:.0f} um. Regolith wall y0 (contact-angle BC from the experiment);
# y1 is not a wall (90-degree phi face, sealed, insulating). Caps start at theta =
# {THETA_WATER:g} deg (liquid water on regolith, user 2026-11-09, placeholder value): each is a
# disc of arc radius R centered R cos(theta) below the wall (-ice_grain_*):
{lines}# eps {eps:.4e} m, h = eps/sqrt(2): Nx {nx}, Ny {ny}. dtmax = 2 tau_sub = {dtmax:.4e} s.
# =============================================================================
-ic_type multi_grains
-dim 2
-Nx {nx}
-Ny {ny}
-Nz 1
-Lx {lx:.6e}
-Ly {ly:.6e}
-Lz 0
-ice_grain_cx {','.join(f'{x:.6e}' for x, _ in caps)}
-ice_grain_cy {','.join(f'{-R*ct:.6e}' for _, R in caps)}
-ice_grain_R  {','.join(f'{R:.6e}' for _, R in caps)}
-eps {eps:.4e}
-eps_valid_temp {T0:g}
-delt_t 1.0e-4
-periodic 0
-wall_faces y0
-dtmax {dtmax:.4e}
"""
    write(GEO / f"{name}.opts", text, dry)
    print(f"    Nx {nx} Ny {ny} ({(nx+2)*(ny+2)*3/1e6:.2f} M DOF)")
    return name


# --------------------------------------------------------------- throat (5)
def throat_geometry(name, wt, lens_c, eps, purpose, dry, caps=(), grains=()):
    lx, w, lens_len = 8.0e-4, 1.5e-4, 1.2e-4
    spec = {
        "name": name, "eps": eps, "Lx": lx, "mean_gap": w,
        "comment": purpose,
        "walls": {"flat": True},
        "throat": {"x": 4.0e-4, "R": 2.5e-4, "gap": wt},
        "theta_ice_init": THETA_WATER,
        "plugs": [{"x0": lens_c - 0.5 * lens_len, "x1": lens_c + 0.5 * lens_len, "theta": THETA_BC}],
        "caps": [{"x": x, "wall": wall, "R": R} for x, wall, R in caps],
        "grains": [{"x": x, "y": y, "R": R} for x, y, R in grains],
        "dtmax": 2 * tau_sub(eps),
    }
    lay = LAYOUTS / f"{name}.json"
    print(f"  {'(dry) ' if dry else ''}wrote {lay.relative_to(ROOT)}")
    if dry:
        return name
    lay.write_text(json.dumps(spec, indent=2) + "\n")
    r = subprocess.run([sys.executable, str(ROOT / "preprocess/build_ice_layout.py"), str(lay)],
                       capture_output=True, text=True, cwd=ROOT)
    print("    " + r.stdout.replace("\n", "\n    ").rstrip())
    if r.returncode != 0:
        raise SystemExit(f"build_ice_layout failed for {name} (rc {r.returncode}):\n{r.stderr[-1500:]}")
    return name


# -------------------------------------------------------------- experiments
def experiment(name, lx, G, vapor, days, purpose, dry, per_day=1, h=1.0, chi=0.0, dh=0.0, n_out=None):
    """G > 0: warm top (x = 0), T fixed on x0, x1; G = 0: isothermal (T fixed on every face).

    n_out: number of time-uniform snapshots; default days * per_day + 1."""
    if G:
        therm = f"-grad_temp0 {-G:.6g},0.0,0.0\n-flag_BC_Tfix 1\n-Tfix_faces x0,x1\n"
        tdesc = f"gradient {G:g} K/m, warm top (x = 0), T fixed on x0 and x1, other faces insulating"
    else:
        therm = "-grad_temp0 0.0,0.0,0.0\n-flag_BC_Tfix 1\n"
        tdesc = "isothermal at -20 C (T fixed on the boundary)"
    if vapor == "closed":
        vap = "-flag_BC_rhovfix 0                      # no vapor in or out: water conserved\n"
    elif vapor == "ls":
        lo, hi, Tlo, Thi, _ = fracs(T0, -G, lx, "local-sat")
        vap = (f"-flag_BC_rhovfix 1\n-rhovfix_faces x0,x1\n"
               f"-rhovfix_lo {lo:.10f}                 # rho_vs({Tlo:.4f} C) / rho_vs(-20 C)\n"
               f"-rhovfix_hi {hi:.10f}                 # rho_vs({Thi:.4f} C) / rho_vs(-20 C)\n")
    elif vapor == "res":
        lo, hi, *_ = fracs(T0, -G, lx, "h-open", h, chi)
        vap = (f"-flag_BC_rhovfix 1\n-rhovfix_faces x0,x1\n"
               f"-rhovfix_lo {lo:.10f}                 # h = {h:g}, 1 + d0 chi with chi = {chi:.4e} 1/m\n"
               f"-rhovfix_hi {hi:.10f}\n")
    elif vapor == "dh":
        vap = (f"-flag_BC_rhovfix 1\n-rhovfix_faces x0,x1\n"
               f"-rhovfix_lo {1+0.5*dh:.10f}                 # 1 + dh/2, dh = {dh:.6e}\n"
               f"-rhovfix_hi {1-0.5*dh:.10f}                 # 1 - dh/2\n")
    else:
        raise ValueError(vapor)
    gis = gamma_is(THETA_BC)
    nout = n_out or int(days * per_day) + 1
    text = f"""# experiment/cryosuction/{name}.opts   (make_type125_inputs.py, {TODAY})
#
# {purpose}
# Thermal: {tdesc}. Vapor set "{vapor}". Wall angle theta = {THETA_BC:g} deg (D1) from
# Young (gamma_is = 0.3 - 0.109 cos theta). {days:g} days, {nout} snapshots ({(nout - 1) / days:.4g} a day).
# -dtmax is in the GEOMETRY file (2 tau_sub depends on its eps).
#
-temp -20.0
-humidity 1.0
{therm}{vap}-phi_open_bc 0                          # 90-degree phi face on every non-wall face
-t_final {days*DAY:.6e}
-n_out {nout}
-beta_sub0 {BETA0:.4e}                     # K&P beta0 at -20 C, alpha_c = 1e-2
-gamma_ia {G_IA:g}
-gamma_as {G_AS:g}
-gamma_is {gis:.7f}                     # cos(theta) = {math.cos(math.radians(THETA_BC)):+.6f}
-dtCFL 1
-outp 0
""" + ("""-outp_dphi 0.5                          # extra snapshot when phi changed by 0.5 anywhere (solver e5a593e7)
-outp_dphi_min_dt 3600                  # at most one such snapshot an hour
""" if n_out else "")
    write(EXP / f"{name}.opts", text, dry)
    return name


def tests_file(path, title, rows, dry):
    lines = [f"# {title} (make_type125_inputs.py, {TODAY}). RUNS.md IDs in the comments."]
    lines += [f"{g}:{e}   # {i}" for i, g, e in rows]
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}  ({len(rows)} runs)")
    if not dry:
        path.write_text("\n".join(lines) + "\n")


def gname(t, s):
    return f"t{t}_{s}"


# ------------------------------------------------- gradients, types 1 and 2
# Mirrors docs/run_plan/build_runplan2.py (run plan of 2026-11-12): production
# gradients G_PROD (user, 2026-11-12), a few runs at G_HIGH, and the single-G
# series (position, size, vapor set, chains, spacing) at G_ONE. Every type 1
# and type 2 experiment writes N_OUT time-uniform snapshots (-outp 0) for movies.
G_PROD = (10, 20)
G_HIGH = 50
G_ONE = G_PROD[1]
N_OUT = 401
R0_RES = 5.0e-5          # 1B reservoir set: both ends at 1 + d0/R0 (the 50 um cap)


def thermal_tag(G):
    return f"G{G:g}warmtop" if G else "iso"


def exp12(t, lx, G, vapor, days, purpose, dry, l900=False, chi=0.0, dh_of=None):
    """Type 1/2 experiment named from its settings; vapor 'res' takes chi, 'dh' the G it twins."""
    if vapor == "dh":
        lo, hi, *_ = fracs(T0, -dh_of, lx, "local-sat")
        mid, kw = f"iso_dh-G{dh_of:g}", {"dh": lo - hi}
    else:
        tag = {"closed": "closed", "ls": "ls"}.get(vapor) or ("resR50" if t == 1 else "resavg")
        mid, kw = f"{thermal_tag(G)}_{tag}", ({"chi": chi} if vapor == "res" else {})
    name = f"t{t}_T-20_{mid}{'_L900' if l900 else ''}_theta60_{days:g}d_ac1e-2"
    # 10-day runs: 321 snapshots keep the interval (2700 s) above 2 dtmax (cryosuction_main.c output check)
    return experiment(name, lx, G, vapor, days, purpose, dry, n_out=N_OUT if days > 10 else 321, **kw)


# ------------------------------------------------------------------ type 1
# (ID, geometry tag, G [K/m], vapor set, days, dh twin of G): the plan's type 1 table.
TYPE1 = [
    ("1A-0", "R50_x300", 0, "closed", 10, None),
    ("1A-1", "R50_x300", G_PROD[0], "closed", 60, None),
    ("1A-2", "R50_x300", G_PROD[1], "closed", 60, None),
    ("1A-3", "R50_x300", G_HIGH, "closed", 30, None),
    ("1B-1", "R50_x500", G_ONE, "res", 60, None),
    ("1B-2", "R50_x300", G_ONE, "res", 60, None),
    ("1B-3", "R50_x100", G_ONE, "res", 60, None),
    ("1B-4", "R50_x500", G_ONE, "ls", 60, None),
    ("1B-5", "R50_x300", G_ONE, "ls", 60, None),
    ("1B-6", "R50_x100", G_ONE, "ls", 60, None),
    ("1B-7", "R50_x300", 0, "res", 10, None),
    ("1B-8", "R50_x300", G_HIGH, "ls", 30, None),
    ("1C-1", "R25_x300", G_ONE, "closed", 60, None),
    ("1C-2", "R100_x300", G_ONE, "closed", 60, None),
    ("1C-3", "R25_x300", G_ONE, "ls", 60, None),
    ("1C-4", "R100_x300", G_ONE, "ls", 60, None),
    ("1D-1", "R50_x300", 0, "dh", 60, G_ONE),     # twin of 1B-5
    ("1D-2", "R50_x300", 0, "dh", 30, G_HIGH),    # twin of 1B-8
]
PURPOSE1 = {
    "closed": "Plan 1A/1C: drift in a closed box; mass constant.",
    "res": "Plan 1B-1..3, 1B-7: both ends at 1 + d0/R0 (R0 = 50 um), draft 1's form ii.",
    "ls": "Plan 1B-4..6, 1B-8, 1C-3/4: local-saturation pair.",
    "dh": "Plan 1D: vapor-gradient twin (isothermal, ends at the local-saturation contrast of the gradient run).",
}


def type1(dry):
    lx, ly = 6.0e-4, 3.0e-4
    print("\nType 1 geometry (eps %.4e, dtmax %.4e)" % (EPS1, 2 * tau_sub(EPS1)))
    geo = {}
    for tag, x, R in (("R50_x300", 3.0e-4, 5.0e-5), ("R50_x100", 1.0e-4, 5.0e-5), ("R50_x500", 5.0e-4, 5.0e-5),
                      ("R25_x300", 3.0e-4, 2.5e-5), ("R100_x300", 3.0e-4, 1.0e-4)):
        geo[tag] = cap_geometry(f"t1_box600x300um_cap{tag}um_th40_eps0.86um", lx, ly, [(x, R)], EPS1,
                                f"Plan type 1: one cap, arc radius {R/um:.0f} um, centered at x = {x/um:.0f} um.", dry)
    print("\nType 1 experiments")
    made, rows = {}, []
    for rid, gtag, G, vapor, days, twin in TYPE1:
        key = (G, vapor, days, twin)
        if key not in made:
            purpose = PURPOSE1[vapor] if (G or vapor == "dh") else (
                "Plan 1A-0: control, the cap relaxes to theta and stays." if vapor == "closed"
                else "Plan 1B-7: isothermal, both ends at 1 + d0/R0: unstable balance.")
            made[key] = exp12(1, lx, G, vapor, days, purpose, dry, chi=1.0 / R0_RES, dh_of=twin)
        rows.append((rid, geo[gtag], made[key]))
    tests_file(STUDY / "type1_tests.txt", "Plan type 1", rows, dry)


# ------------------------------------------------------------------ type 2
def type2(dry):
    print("\nType 2 geometry (eps %.4e, dtmax %.4e)" % (EPS25, 2 * tau_sub(EPS25)))
    ly = 3.0e-4
    b40 = 4.0e-5 * math.sin(math.radians(THETA_WATER))
    rng = __import__("numpy").random.default_rng(7)
    r5 = [float(r) for r in rng.uniform(2.0e-5, 6.0e-5, 5)]
    geo = {
        "pairLW": (6e-4, [(2e-4, 5e-5), (4e-4, 2e-5)], "ripening pair, large cap (50 um) on the warm side"),
        "pairLC": (6e-4, [(2e-4, 2e-5), (4e-4, 5e-5)], "ripening pair, large cap on the cold side (2B-6)"),
        "chain3eq": (9e-4, [(2.5e-4, 4e-5), (4.5e-4, 4e-5), (6.5e-4, 4e-5)], "three equal caps (40 um)"),
        "chain3small": (9e-4, [(2.5e-4, 4e-5), (4.5e-4, 2e-5), (6.5e-4, 4e-5)], "three caps, small center"),
        "chain5rand": (9e-4, [(1.5e-4 * (k + 1), r) for k, r in enumerate(r5)],
                       "five caps, arc radii uniform 20-60 um (numpy default_rng(7))"),
        **{f"gap{g}": (6e-4, [(3e-4 - b40 - 0.5 * g * 1e-6, 4e-5), (3e-4 + b40 + 0.5 * g * 1e-6, 4e-5)],
                       f"two 40 um caps, {g} um apart at the wall") for g in (20, 80, 200)},
    }
    names = {}
    for tag, (lx, caps, what) in geo.items():
        names[tag] = cap_geometry(f"t2_box{lx/um:.0f}x300um_{tag}_th40_eps0.80um", lx, ly, caps, EPS25,
                                  f"Plan type 2: {what}.", dry)
    print("\nType 2 experiments")
    chi_avg = 0.5 * (1 / 2e-5 + 1 / 5e-5)
    ladder = (0.8, 1.6, 3.2, 8.0, 32.0)           # 2B: around G* = 1.6 K/m, kept from 2026-11-09
    e = {
        "iso_closed": exp12(2, 6e-4, 0, "closed", 60, "Plan 2A-1: ripening alone.", dry),
        "iso_res": exp12(2, 6e-4, 0, "res", 60, "Plan 2A-2: draft 1's open version, both ends at 1 + d0 chi_avg.",
                         dry, chi=chi_avg),
        **{f"G{g:g}": exp12(2, 6e-4, g, "closed", 60, f"Plan 2B: closed box at G = {g:g} K/m (G* = 1.6 K/m).", dry)
           for g in ladder},
        "chain": exp12(2, 9e-4, G_ONE, "closed", 60, f"Plan 2C: chains at G = {G_ONE} K/m (900 um box).", dry,
                       l900=True),
        "gap": exp12(2, 6e-4, G_ONE, "closed", 60, f"Plan 2D: spacing at G = {G_ONE} K/m.", dry),
    }
    rows = [("2A-1", names["pairLW"], e["iso_closed"]), ("2A-2", names["pairLW"], e["iso_res"])]
    rows += [(f"2B-{i}", names["pairLW"], e[f"G{g:g}"]) for i, g in enumerate(ladder, 1)]
    rows += [("2B-6", names["pairLC"], e["G3.2"])]
    rows += [(f"2C-{i}", names[t], e["chain"]) for i, t in enumerate(("chain3eq", "chain3small", "chain5rand"), 1)]
    rows += [(f"2D-{i}", names[f"gap{g}"], e["gap"]) for i, g in enumerate((20, 80, 200), 1)]
    tests_file(STUDY / "type2_tests.txt", "Plan type 2", rows, dry)


# ------------------------------------------------------------------ type 5
def type5(dry):
    print("\nType 5 geometry (eps %.4e, dtmax %.4e)" % (EPS25, 2 * tau_sub(EPS25)))
    lx = 8.0e-4
    g = {
        "wt40": throat_geometry("t5_pore800x150um_wt40um_lens120c_eps0.80um", 4e-5, 4e-4, EPS25,
                                "Plan 5A-1..5: 40 um throat, lens (120 um on the axis) centered in it.", dry),
        "wt80": throat_geometry("t5_pore800x150um_wt80um_lens120c_eps0.80um", 8e-5, 4e-4, EPS25,
                                "Plan 5A-6..8: 80 um throat, lens centered.", dry),
        "above": throat_geometry("t5_pore800x150um_wt40um_lens120up_eps0.80um", 4e-5, 3e-4, EPS25,
                                 "Plan 5C-1: lens starts 100 um above the throat (toward the warm top).", dry),
        "below": throat_geometry("t5_pore800x150um_wt40um_lens120dn_eps0.80um", 4e-5, 5e-4, EPS25,
                                 "Plan 5C-2: lens starts 100 um below the throat.", dry),
        "5B": throat_geometry("t5_pore800x150um_wt40um_lens120c_caps_grain_eps0.80um", 4e-5, 4e-4, EPS25,
                              "Plan 5B: lens + 2 caps + 1 grain.", dry,
                              caps=[(1.8e-4, "bottom", 4.0e-5), (6.4e-4, "top", 3.0e-5)],
                              grains=[(6.6e-4, 5.5e-5, 2.0e-5)]),
    }
    print("\nType 5 experiments")
    e = {"res0": experiment("t5_T-20_iso_res1_theta60_120d_ac1e-2", lx, 0, "res", 120,
                            "Plan 5A/5C at G = 0: both ends at h = 1.", dry)}
    for G in (1, 3, 10, 30):
        e[f"ls{G}"] = experiment(f"t5_T-20_G{G}warmtop_ls_theta60_120d_ac1e-2", lx, G, "ls", 120,
                                 f"Plan 5A/5C: local-saturation pair at G = {G} K/m.", dry)
    e["closed0"] = experiment("t5_T-20_iso_closed_theta60_120d_ac1e-2", lx, 0, "closed", 120,
                              "Plan 5B-1: closed, isothermal.", dry)
    e["closed3"] = experiment("t5_T-20_G3warmtop_closed_theta60_120d_ac1e-2", lx, 3, "closed", 120,
                              "Plan 5B-2: closed, G = 3 K/m.", dry)
    rows = [(f"5A-{i}", g["wt40"], e[k]) for i, k in enumerate(("res0", "ls1", "ls3", "ls10", "ls30"), 1)]
    rows += [(f"5A-{i}", g["wt80"], e[k]) for i, k in enumerate(("res0", "ls3", "ls10"), 6)]
    rows += [("5B-1", g["5B"], e["closed0"]), ("5B-2", g["5B"], e["closed3"])]
    rows += [("5C-1", g["above"], e["ls3"]), ("5C-2", g["below"], e["ls3"])]
    tests_file(STUDY / "type5_tests.txt", "Plan type 5", rows, dry)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", default="1,2,5")
    a = ap.parse_args()
    want = set(a.only.split(","))
    print(f"theta: caps start at {THETA_WATER:g} deg (water placeholder), lenses at {THETA_BC:g} deg (BC, D1)")
    if "1" in want:
        type1(a.dry_run)
    if "2" in want:
        type2(a.dry_run)
    if "5" in want:
        type5(a.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
