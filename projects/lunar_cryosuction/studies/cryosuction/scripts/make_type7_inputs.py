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

Saturated set (--set sat, user 2026-11-17): 21 runs on C1, C4, C7, 120 days. Vapor starts saturated
everywhere (-humidity 1.0 gives rho_vs(T(x)) in the air) and the collector on x0 holds the saturation
density of its OWN temperature (rhov_bc.fracs "local-sat"), so the temperature gradient is the only
driver of harvest. Temperature: both ends fixed at the IC values and the walls y0, y1 held at the same
linear profile (-Twall_faces y0,y1 -Twall_profile linear, symmetric Nitsche), as the regolith around
the pore would set it.

  7S0-1,4,7          G = 0 control (same thermal conditions with dT = 0)
  7SD-G-1,4,7        G = 10, 20, 50 K/m, warm end at the dead end (x = L): ice should move to the collector
  7SC-G-1,4,7        G = 10, 20, 50 K/m, warm end at the collector (x = 0)

Pumped warm collector (--set pump, user 2026-11-18): replaces 7SC. A warm collector saturated at its own
T is the wettest point of the pore and feeds it (7SC ice grew 13-359 %). A heated device that collects must
also remove vapor: the collector holds H_PUMP = 0.99 times the saturation density of its own temperature
(rhov_bc.fracs "h-face"; the 1 % deficit of batch 1, now at the face's own T). 9 runs: C1, C4, C7 x G 10,
20, 50 K/m, warm end at the collector, walls and ends as 7S. Controls: the batch-1 isothermal runs 7A-1/4/7
(0.99 rho_vs(-20 C), which is the face's own T at G = 0).

  7SP-G-1,4,7        G = 10, 20, 50 K/m, warm pumped collector (x = 0)

    venv_cryosuction/bin/python studies/cryosuction/scripts/make_type7_inputs.py [--set v5|sat|pump] [--dry-run]
"""
from __future__ import annotations

import argparse
import datetime
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "inputs/experiment/cryosuction"
STUDY = ROOT / "studies/cryosuction"
TODAY = datetime.date.today().isoformat()
DAYS, N_OUT, G = 60, 401, 20.0
GEOM = "type7v5_C{}_2D_L1800um_eps0.86um"
LX, T0 = 1.8e-3, -20.0                  # pore length (preprocess/layouts/type7v5_C*.json) and mean T
SAT_DAYS, SAT_G, SAT_CONFIGS = 120, (10.0, 20.0, 50.0), (1, 4, 7)
H_PUMP = 0.99                           # pumped collector: humidity at the face's own temperature

COMMON = """-temp -20.0
-humidity 1.0
-flag_BC_rhovfix 1
-rhovfix_faces x0
{collector}-phi_open_bc 0                          # 90-degree phi face on the pore ends
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


D2 = "-rhovfix_lo 0.990                       # collector: 0.990 rho_vs(T0) (decision D2)\n"


def experiment(name, head, thermal, dry, days=DAYS, collector=D2):
    text = (f"# experiment/cryosuction/{name}.opts   (make_type7_inputs.py, {TODAY})\n#\n" + head +
            "# -dtmax is in the GEOMETRY file (2 tau_sub at its eps).\n#\n" + thermal +
            COMMON.format(tf=days * 86400.0, days=days, nout=N_OUT, collector=collector))
    path = EXP / f"{name}.opts"
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}")
    if not dry:
        path.write_text(text)
    return name


def sat_experiment(G, warm, dry, h=None):
    """Saturated set: grad_temp0 = +G (warm dead end, x = L) or -G (warm collector, x = 0); G = 0 control.
    h: pumped collector at h times the saturation density of its own temperature (None: saturated)."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from rhov_bc import fracs
    g = G if warm == "deadend" else -G
    lo, _, Tlo, Thi, _ = fracs(T0, g, LX, "local-sat" if h is None else "h-face", 1.0 if h is None else h)
    tcold, twarm = min(Tlo, Thi), max(Tlo, Thi)
    s0 = LX if g > 0 else 0.0
    if G == 0:
        name, what = f"t7s_T-20_hsat_iso_wallT_theta60_{SAT_DAYS}d_n401_ac1e-2", "G = 0 control"
    else:
        htag = "hsat" if h is None else f"hface{h:.3f}"
        name = f"t7s_T-20_{htag}_G{G:g}warm{warm}_wallT_theta60_{SAT_DAYS}d_n401_ac1e-2"
        what = f"G = {G:g} K/m, warm end at the {'dead end (x = L)' if warm == 'deadend' else 'collector (x = 0)'}"
    holds = ("the saturation density of its own temperature" if h is None
             else f"{h:g} times the saturation density of its own temperature (pumped)")
    head = (("# Plan type 7S (saturated, user 2026-11-17): " if h is None else "# Plan type 7SP (pumped, user 2026-11-18): ") + f"{what}. "
            f"Vapor saturated at start; the collector\n"
            f"# holds {holds} ({Tlo:.6f} C); T fixed on both ends and\n"
            f"# on the walls (linear wall profile from {twarm:.6f} C at x = {s0:g} m to {tcold:.6f} C).\n")
    thermal = (f"-grad_temp0 {g:.1f},0.0,0.0\n-flag_BC_Tfix 1\n-Tfix_faces x0,x1                        # ends at the IC values\n"
               f"-Twall_faces y0,y1                       # walls: T_w = T_init + dT (1 - |x - s0| / L), Nitsche\n"
               f"-Twall_profile linear\n-Twall_axis 0\n-Twall_s0 {s0:.6e}\n-Twall_length {LX:.6e}\n"
               f"-Twall_dT {abs(g) * LX:.6e}                    # G Lx\n-Twall_Tinit {tcold:.10f}             # cold end\n")
    if h is None:
        coll = (f"-rhovfix_lo {lo:.10f}                 # rho_vs({Tlo:.6f} C) / rho_vs(-20 C): saturated at the "
                f"collector's own T (rhov_bc.py --T0 -20 --G {g:g} --Lx {LX:g} --local-sat)\n")
    else:
        coll = (f"-rhovfix_lo {lo:.10f}                 # {h:g} rho_vs({Tlo:.6f} C) / rho_vs(-20 C): pumped collector "
                f"(rhov_bc.py --T0 -20 --G {g:g} --Lx {LX:g} --h-face {h:g})\n")
    return experiment(name, head, thermal, dry, days=SAT_DAYS, collector=coll)


def main_sat(dry) -> int:
    ctrl = sat_experiment(0.0, "deadend", dry)
    rows = [(f"7S0-{i}", ctrl, i) for i in SAT_CONFIGS]
    for warm, tag in (("deadend", "SD"), ("collector", "SC")):
        for G in SAT_G:
            name = sat_experiment(G, warm, dry)
            rows += [(f"7{tag}-{G:g}-{i}", name, i) for i in SAT_CONFIGS]
    lines = [f"# Plan type 7S (make_type7_inputs.py --set sat, {TODAY}): saturated vapor, collector saturated at its own T,",
             "# linear wall temperature; C1, C4, C7 x (G 0; G 10, 20, 50 K/m warm dead end / warm collector). 21 runs, 120 d."]
    lines += [f"{GEOM.format(i)}:{e}   # {rid}" for rid, e, i in rows]
    path = STUDY / "type7sat_tests.txt"
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}  ({len(rows)} runs)")
    if not dry:
        path.write_text("\n".join(lines) + "\n")
    return 0


def main_pump(dry) -> int:
    rows = []
    for G in SAT_G:
        name = sat_experiment(G, "collector", dry, h=H_PUMP)
        rows += [(f"7SP-{G:g}-{i}", name, i) for i in SAT_CONFIGS]
    lines = [f"# Plan type 7SP (make_type7_inputs.py --set pump, {TODAY}): warm collector pumped to {H_PUMP:g} of the",
             "# saturation density of its own T (replaces 7SC); C1, C4, C7 x G 10, 20, 50 K/m. 9 runs, 120 d."]
    lines += [f"{GEOM.format(i)}:{e}   # {rid}" for rid, e, i in rows]
    path = STUDY / "type7pump_tests.txt"
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}  ({len(rows)} runs)")
    if not dry:
        path.write_text("\n".join(lines) + "\n")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--set", choices=("v5", "sat", "pump"), default="v5")
    a = ap.parse_args()
    if a.set == "sat":
        return main_sat(a.dry_run)
    if a.set == "pump":
        return main_pump(a.dry_run)
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
