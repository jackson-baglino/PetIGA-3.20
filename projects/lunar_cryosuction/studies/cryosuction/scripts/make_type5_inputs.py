#!/usr/bin/env python3
"""make_type5_inputs.py: layouts, meshes, experiment files and the tests file for plan type 5
as redesigned on 2026-11-13 (SI draft 3.2, "vapor boundary conditions for the throat
simulations") and approved by the user on 2026-11-17.

Question: how does ice aggregate in the pore space before harvesting? A lens in a pore that
narrows to a throat moves by curvature toward the narrower gap; a temperature gradient pulls it
toward the warm end. Does the throat hold the lens, and where does it come to rest?

Every pore is 400 um long and 150 um at its widest, with regolith walls (contact-angle condition,
theta 60 deg) on y = 0 and y = 150 um. The lens is 100 um long on the axis and starts with its
center 70 um on the COLD side of the throat (x = 270 um; x = 0 is the warm end). eps 0.80 um
(half the 40 um throat / 25, the type 5 rule of make_type125_inputs.py), dtmax = 2 tau_sub(eps).

  5P   periodic cell (-periodic_x 1), sinhyp walls, G = 0, no vapor in or out (water conserved):
       5P-1 throat 40 um, 5P-2 throat 80 um, 5P-0 straight periodic channel (control: the lens
       should not move).
  5G   hyperbolic pore g = sqrt(wt^2 + c^2 xi^2) ("hyp", g = W at both ends), T fixed at both
       ends with a gradient G (warm at x = 0), each end held at h_net rho_vs(T_end): h_net is the
       constant end humidity that returns the lens to its initial area after 300 d
       (type5_model.EndsPore.h_net, SI Eq. S6). G = 0.5 G_crit (the throat holds the lens) and
       2 G_crit (it pushes the lens through), G_crit from type5_model for each throat, rounded to
       0.01 K/m (user, 2026-11-17).
  5S   straight pore 150 um wide, same lens, the same four G values, each with its own h_net:
       how far the ends and the gradient alone move the lens.

All runs 300 days (user: extend by 300 d with resume_run.sh if needed), 401 time-uniform
snapshots plus event snapshots (-outp_dphi 0.5).

    venv_cryosuction/bin/python studies/cryosuction/scripts/make_type5_inputs.py [--dry-run]
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
from make_type34_inputs import BETA0, DAY, EXP, G_AS, G_IA, STUDY, T0, gamma_is, tau_sub  # noqa: E402
from rhov_bc import fracs  # noqa: E402
from type5_model import EndsPore, consts, side_words  # noqa: E402

LX, W, LENS, START, THETA = 4.0e-4, 1.5e-4, 1.0e-4, 7.0e-5, 60.0
THROATS = (4.0e-5, 8.0e-5)
GFRAC = (0.5, 2.0)
EPS = 8.0e-7
DAYS, N_OUT = 300.0, 401
LAYOUTS = ROOT / "preprocess/layouts"
TODAY = datetime.date.today().isoformat()
um = 1e-6


def write(path: Path, text: str, dry: bool):
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}")
    if not dry:
        path.write_text(text)


def geometry(name, walls, periodic, purpose, dry):
    xc = 0.5 * LX + START
    spec = {"name": name, "eps": EPS, "Lx": LX, "comment": purpose, "walls": walls,
            "periodic": periodic, "theta_ice_init": THETA, "caps": [], "grains": [],
            "plugs": [{"x0": xc - 0.5 * LENS, "x1": xc + 0.5 * LENS, "end_faces": True, "faces": "local"}],
            "dtmax": 2 * tau_sub(EPS)}
    if walls.get("flat"):
        spec["mean_gap"] = W
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


def experiment(name, G, h, purpose, dry):
    """G = 0 and h None: periodic, closed. Otherwise: T fixed on x0/x1 with grad_temp0 = -G, ends at
    h rho_vs(T_end) (rhov_bc.py --h-face)."""
    if h is None:
        body = ("-grad_temp0 0.0,0.0,0.0\n-flag_BC_Tfix 1                        # T fixed on the walls"
                " (the x faces are a periodic seam): isothermal\n"
                "-flag_BC_rhovfix 0                      # no vapor in or out: water conserved\n")
        desc = "isothermal at -20 C, periodic in x, no vapor in or out"
    else:
        lo, hi, Tlo, Thi, _ = fracs(T0, -G, LX, "h-face", h)
        body = (f"-grad_temp0 {-G:.6g},0.0,0.0\n-flag_BC_Tfix 1\n-Tfix_faces x0,x1\n"
                f"-flag_BC_rhovfix 1\n-rhovfix_faces x0,x1\n"
                f"-rhovfix_lo {lo:.10f}                 # h_net rho_vs({Tlo:.6f} C) / rho_vs(-20 C), warm end\n"
                f"-rhovfix_hi {hi:.10f}                 # h_net rho_vs({Thi:.6f} C) / rho_vs(-20 C), cold end\n")
        desc = (f"gradient {G:g} K/m, warm at x = 0, T fixed on x0 and x1, walls insulating; each end at "
                f"h_net = {h:.10f} times rho_vs of its own temperature (rhov_bc.py --T0 -20 --G {-G:g} "
                f"--Lx {LX:g} --h-face {h:.10f})")
    text = f"""# experiment/cryosuction/{name}.opts   (make_type5_inputs.py, {TODAY})
#
# {purpose}
# {desc}.
# Wall angle theta = {THETA:g} deg from Young (gamma_is = 0.3 - 0.109 cos theta). {DAYS:g} days,
# {N_OUT} snapshots. -dtmax is in the GEOMETRY file (2 tau_sub depends on its eps).
# Use only with the 400 um type-5 pores (the end values depend on Lx).
#
-temp -20.0
-humidity 1.0
{body}-phi_open_bc 0                          # 90-degree phi face on every non-wall face
-t_final {DAYS * DAY:.6e}
-n_out {N_OUT}
-beta_sub0 {BETA0:.4e}                     # K&P beta0 at -20 C, alpha_c = 1e-2
-gamma_ia {G_IA:g}
-gamma_as {G_AS:g}
-gamma_is {gamma_is(THETA):.7f}                     # cos(theta) = {math.cos(math.radians(THETA)):+.6f}
-dtCFL 1
-outp 0
-outp_dphi 0.5                          # extra snapshot when phi changed by 0.5 anywhere (solver e5a593e7)
-outp_dphi_min_dt 3600                  # at most one such snapshot an hour
"""
    write(EXP / f"{name}.opts", text, dry)
    return name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    dry = a.dry_run
    k = consts(T0, 1e-2)
    print(f"Type 5: Lx {LX/um:g} um, W {W/um:g} um, lens {LENS/um:g} um at +{START/um:g} um (cold side), "
          f"theta {THETA:g}, eps {EPS:.4e}, dtmax {2 * tau_sub(EPS):.4e} s, {DAYS:g} d")

    print("\nGeometry")
    tag = f"lens100s70_eps{EPS / um:.2f}um"
    gP = {wt: geometry(f"t5_cell400x150um_wt{wt / um:.0f}um_sinhyp_{tag}", {"profile": "sinhyp", "W": W, "wt": wt},
                       True, f"Type 5P: periodic sinhyp cell, throat {wt / um:.0f} um at Lx/2.", dry)
          for wt in THROATS}
    gP0 = geometry(f"t5_cell400x150um_flat_{tag}", {"flat": True}, True,
                   "Type 5P-0: periodic straight channel, 150 um (control).", dry)
    gG = {wt: geometry(f"t5_pore400x150um_wt{wt / um:.0f}um_hyp_{tag}", {"profile": "hyp", "W": W, "wt": wt},
                       False, f"Type 5G: hyperbolic pore, throat {wt / um:.0f} um at Lx/2, ends 150 um.", dry)
          for wt in THROATS}
    gS = geometry(f"t5_pore400x150um_flat_{tag}", {"flat": True}, False,
                  "Type 5S: straight pore, 150 um (control for 5G).", dry)

    print("\nG_crit and h_net (type5_model.EndsPore, 300 d)")
    Gs, rows = {}, []
    for wt in THROATS:
        pore = EndsPore(wt, k, W=W, L=LX, theta=THETA, lens=LENS, start=START)
        gc = pore.g_crit()
        Gs[wt] = [round(f * gc, 2) for f in GFRAC]
        print(f"  throat {wt / um:.0f} um: G_crit {gc:.4f} K/m -> G {Gs[wt]}")
    straight = EndsPore(W, k, W=W, L=LX, theta=THETA, lens=LENS, start=START)

    print("\nExperiments")
    e0 = experiment(f"t5_T-20_iso_periodic_closed_theta60_{DAYS:g}d_ac1e-2", 0.0, None,
                    "Type 5P: periodic cell at G = 0, water conserved.", dry)
    rows += [("5P-1", gP[THROATS[0]], e0), ("5P-2", gP[THROATS[1]], e0), ("5P-0", gP0, e0)]
    count = {"G": 0, "S": 0}
    for wt in THROATS:
        pore = EndsPore(wt, k, W=W, L=LX, theta=THETA, lens=LENS, start=START)
        for f, G in zip(GFRAC, Gs[wt]):
            for kind, model, geo in (("G", pore, gG[wt]), ("S", straight, gS)):
                h, passes = model.h_net(G, DAYS)
                r = model.run(h, G, DAYS)
                ptag = f"wt{wt / um:.0f}" if kind == "G" else "flat"
                name = f"t5_T-20_G{G:g}warmtop_hnet-{ptag}_theta60_{DAYS:g}d_ac1e-2"
                where = ("" if r["t"] >= DAYS * DAY - 1 else f", a face reaches an end at {r['t'] / DAY:.0f} d")
                purpose = (f"Type 5{kind}: {'throat ' + str(round(wt / um)) + ' um' if kind == 'G' else 'straight pore'}"
                           f", G = {G:g} K/m ({f:g} G_crit of the {wt / um:.0f} um throat). Model (type5_model): "
                           f"h_net {h:.10f} ({passes} passes), lens center after {DAYS:g} d "
                           f"{side_words(r['xm'])}{where}")
                print(f"  {purpose}")
                experiment(name, G, h, purpose, dry)
                count[kind] += 1
                rows.append((f"5{kind}-{count[kind]}", geo, name))
    lines = [f"# Plan type 5, redesign of 2026-11-13, approved 2026-11-17 (make_type5_inputs.py, {TODAY}).",
             "# RUNS.md IDs in the comments."]
    lines += [f"{g}:{e}   # {i}" for i, g, e in rows]
    tf = STUDY / "type5_tests.txt"
    print(f"  {'(dry) ' if dry else ''}wrote {tf.relative_to(ROOT)}  ({len(rows)} runs)")
    if not dry:
        tf.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
