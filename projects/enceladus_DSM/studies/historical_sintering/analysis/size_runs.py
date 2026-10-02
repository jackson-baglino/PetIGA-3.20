#!/usr/bin/env python3
"""size_runs.py -- mesh, timestep, run length and cost for replicating the
Kingery (1960) and Thomas et al. (1994) two-sphere experiments.

Run from enceladus_DSM/:
    python studies/historical_sintering/analysis/size_runs.py

METHOD
------
eps   : set by the neck-resolution floor  r/R >= sqrt(12*eps/R)  (see
        studies/sinter_exponent/PLAN.md), placed at 0.9x the experiment's first
        measured relative neck so that point is resolved with margin:
            eps = R * (0.9*u_first)^2 / 12
        The K&P bounds from comp_eps.py are 1-2 decades looser at alpha_c = 1e-3,
        so the floor is what binds.
mesh  : comp_eps.py rule h = eps/sqrt(2); domain Lx = 4.4R (two tangent grains
        + 0.2R padding each end, as in the Molaro tangent geometry), Ly = 1.2R.
t     : model neck trajectory t(u) calibrated on the mesh_pair FINE arm
        (alpha_c = 1e-3, -20 C, tangent start, R_eff = 84.4 um), which went
        u = 0.194 -> 0.303 in 15.05 -> 78.6 h, i.e. t ~ u^3.73. Scaled to each
        case by the kinetic-limit law t ~ beta_sub * R^2 / d0.
steps : t_final / dtmax + 100, with dtmax = 2*tau_sub (project rule since
        2026-10-02) and t_final = 1.25 x the model time to reach the last point.
cost  : 3 dof/node, 200k dof/rank (scripts/lib/alloc.sh), 8/15/23 s per step,
        $0.012 per core-hour (Tier 1, scripts/HPC/hpc_cost.sh).
"""

import csv
import math
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                      # enceladus_DSM/

ALPHA = 1e-3
R_REF, U_REF, T_REF_H = 84.4e-6, 0.3027, 78.6   # mesh_pair fine arm, end point
M_EXP = 3.73                                     # t ~ u^m over u = 0.19-0.30
S_PER_STEP = (8.0, 15.0, 23.0)
DOF_PER_RANK = 2.0e5
USD_PER_CORE_H = 0.012

CASES = [
    # name, T [C], R [m], first u, last u, first t [h], last t [h], floor anchor
    ("Kingery 1960",        -17.8, 110e-6, 0.2485, 0.3329, 0.0269, 0.1384, "pt1"),
    ("Thomas 1994",         -20.0, 120e-6, 0.0992, 0.4433, 3.3163, 190.05, "pt1"),
    ("Thomas 1994 (pt3)",   -20.0, 120e-6, 0.1811, 0.4433, 20.153, 190.05, "pt3"),
]


def comp_eps(T, R, eps, Lx, Ly):
    """Run comp_eps.py with Eq. 45 pinned so that it returns exactly `eps`."""
    out = subprocess.run(
        [sys.executable, str(ROOT / "preprocess/comp_eps.py"),
         "--Lx", str(Lx), "--Ly", str(Ly), "--Rave", str(R), "--T0", str(T),
         "--alpha", str(ALPHA), "--Dchannel", "ice", "--vn_feature", str(2 * eps)],
        capture_output=True, text=True, check=True).stdout

    def grab(pat):
        return float(re.search(pat + r"\s*([0-9.eE+-]+)", out).group(1))

    return dict(tau=grab(r"τ_sub  \(relaxation\)  ="),
                Nx=int(grab(r"Nx = ceil\(Lx·√2/ε\)\s*=")),
                Ny=int(grab(r"Ny = ceil\(Ly·√2/ε\)\s*=")),
                beta=grab(r"-beta_sub0"), d0=grab(r"-d0_sub0"))


def main():
    ref = comp_eps(-20.0, R_REF, 0.24e-6, 4.4 * R_REF, 1.2 * R_REF)
    rows = []
    for name, T, R, uf, ul, tf, tl, anchor in CASES:
        u_floor = 0.9 * uf
        eps = R * u_floor**2 / 12
        Lx, Ly = 4.4 * R, 1.2 * R
        p = comp_eps(T, R, eps, Lx, Ly)
        scale = (R / R_REF)**2 * (p["beta"] / p["d0"]) / (ref["beta"] / ref["d0"])
        t_of = lambda u: T_REF_H * (u / U_REF)**M_EXP * scale
        t_final = 1.25 * t_of(ul)
        dtmax = 2 * p["tau"]
        steps = t_final * 3600 / dtmax + 100
        nodes = p["Nx"] * p["Ny"]
        ranks = math.ceil(3 * nodes / DOF_PER_RANK)
        row = dict(case=name, T_C=T, R_um=R * 1e6, anchor=anchor,
                   u_floor=round(u_floor, 4), eps_um=round(eps * 1e6, 4),
                   h_nm=round(eps / math.sqrt(2) * 1e9, 1),
                   Lx_um=round(Lx * 1e6), Ly_um=round(Ly * 1e6),
                   Nx=p["Nx"], Ny=p["Ny"], nodes_M=round(nodes / 1e6, 2),
                   ranks=ranks, tau_sub_s=round(p["tau"], 1), dtmax_s=round(dtmax, 1),
                   t_model_first_h=round(t_of(uf), 2), t_model_last_h=round(t_of(ul), 1),
                   exp_window_h=round(tl - tf, 3),
                   model_window_h=round(t_of(ul) - t_of(uf), 1),
                   slowdown=round((t_of(ul) - t_of(uf)) / (tl - tf), 1),
                   t_final_h=round(t_final, 1), steps=round(steps))
        for s in S_PER_STEP:
            wall = steps * s / 3600
            row[f"wall_h@{s:g}s"] = round(wall, 1)
            row[f"usd@{s:g}s"] = round(wall * ranks * USD_PER_CORE_H, 1)
        rows.append(row)

    out = HERE.parent / "data" / "run_sizing.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
