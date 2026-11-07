"""Plan type 3 (lens in a flat pore): measured face speed against the series law.

For each run folder of a type 3 batch:

  A(t)    ice area per unit depth, tot_ice in SSA_evo.dat (every step),
          resampled on a 0.5-day grid
  com(t)  ice center of mass on the pore axis (ice_bodies.csv, daily snapshots)
  U       area-based speed of an exchanging face, positive = growth:
            both ends held (3A, 3B):  U = (dA/dt) / (2H)
            x0 end only (3D):          U = (dA/dt) / H     (the x1 pocket is sealed)
  l       distance from the x0-side face to the x0 end, as the face of the
          rectangle with the same area and centroid: l = com_x - A/(2H)
          (= (Lx - A/H)/2 for a centered lens)

Theory (plan Section 3, velocity-study README), no fitted parameter:

    v = d0 chi / (beta0 + K l),  chi = 2 cos(theta) / H,  K = rho_ice / (D_v rho_vs)

with theta from Young (cos theta = (gamma_as - gamma_is)/gamma_ia, staged opts),
beta0 = -beta_sub0, d0, D_v(T) and rho_vs(T) as the solver computes them
(preprocess/mesh_params, pplib). A second prediction uses the contact angle
MEASURED at that time (contact_angle.csv, mean of the four contact lines). The
effective kinetic resistance the run implies is

    beta_eff = d0 chi / U - K l      (compare with beta0)

Reported in the window --t0 .. end (default 40 d, after the 90 deg -> theta
shape relaxation, which 3C-1 shows lasts ~30 d) and at fixed l (--at-l).
Prints every input it used. Writes <out>/type3_speeds.csv and
<out>/type3_speeds.png (checked with figcheck).

    venv_cryosuction/bin/python studies/cryosuction/scripts/analyze_type3.py \
        --batch <batch folder> --out <folder> [--t0 40] [--at-l 1.6e-4]
"""
from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "postprocess"))
sys.path.insert(0, str(ROOT))
import pplib  # noqa: E402
from preprocess.mesh_params.physics import (  # noqa: E402
    capillary_length_m, dv_m2_per_s, rho_vs_sat_kg_per_m3)

DAY = 86400.0


def run_meta(run):
    o = pplib.read_opts(run)
    f = lambda k, d=float("nan"): pplib.opt_float(o, k, d)  # noqa: E731
    cos_th = (f("-gamma_as") - f("-gamma_is")) / f("-gamma_ia")
    faces = o.get("-rhovfix_faces", "")
    if not pplib.opt_bool(o, "-flag_BC_rhovfix"):
        mode = "closed"
    elif faces.strip() == "x0":
        mode = "x0"
    else:
        mode = "both"
    return dict(theta=math.degrees(math.acos(cos_th)), cos=cos_th, H=f("-Ly"), Lx=f("-Lx"),
                beta0=f("-beta_sub0"), T=f("-temp"), mode=mode)


def series(run, m, t0_d):
    ssa = pplib.load_ssa(run)
    order = np.argsort(ssa[:, pplib.TIME])
    t, A = ssa[order, pplib.TIME], ssa[order, pplib.ICE]
    tg = np.arange(0.0, t[-1], 0.5 * DAY)
    Ag = np.interp(tg, t, A)
    dA = np.gradient(Ag, tg)
    b = np.genfromtxt(os.path.join(run, "ice_bodies.csv"), delimiter=",", names=True)
    com = np.interp(tg, b["t_s"], b["com_x"])
    H = m["H"]
    U = dA / (2 * H) if m["mode"] == "both" else dA / H
    ell = com - Ag / (2 * H)
    ca = os.path.join(run, "contact_angle.csv")
    th_meas = np.full_like(tg, np.nan)
    if os.path.isfile(ca):
        import pandas as pd
        c = pd.read_csv(ca, comment="#")
        g = c.groupby("time_s")["theta_deg"].mean()
        th_meas = np.interp(tg, g.index.to_numpy(), g.to_numpy())
    return tg, Ag, U, ell, th_meas


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--batch", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--t0", type=float, default=40.0, help="start of the averaging window [d]")
    p.add_argument("--at-l", type=float, default=1.6e-4, help="fixed face distance l [m]")
    a = p.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    runs = sorted(d for d in glob.glob(os.path.join(a.batch, "*__*"))
                  if os.path.isfile(os.path.join(d, "ice_bodies.csv")))
    rows, curves = [], []
    for run in runs:
        m = run_meta(run)
        if m["mode"] == "closed":
            continue
        d0 = capillary_length_m(m["T"])
        K = pplib.RHO_ICE / (dv_m2_per_s(m["T"]) * rho_vs_sat_kg_per_m3(m["T"]))
        beta0 = abs(m["beta0"])
        tg, Ag, U, ell, thm = series(run, m, a.t0)
        chi = 2 * m["cos"] / m["H"]
        v_th = d0 * chi / (beta0 + K * ell)
        v_thm = d0 * 2 * np.cos(np.radians(thm)) / m["H"] / (beta0 + K * ell)
        w = (tg >= a.t0 * DAY) & (tg <= tg[-1] - DAY)
        name = os.path.basename(run).split("__")[1].replace("_T-20", "").replace("_ac1e-2", "")
        name = os.path.basename(run).split("_L500")[0].replace("channel_2D_", "") + " " + name
        r = dict(run=name, mode=m["mode"], theta_young=m["theta"], H_um=m["H"] * 1e6,
                 A0_um2=Ag[0] * 1e12, A_end_um2=Ag[-1] * 1e12,
                 dA_frac=(Ag[-1] - Ag[0]) / Ag[0])
        if abs(m["cos"]) < 1e-9:
            r.update(U_mean_nm_d=np.mean(U[w]) * 1e9 * DAY, v_young_nm_d=0.0)
        else:
            r.update(U_mean_nm_d=np.mean(U[w]) * 1e9 * DAY,
                     v_young_nm_d=np.mean(v_th[w]) * 1e9 * DAY,
                     ratio_young=np.mean(U[w]) / np.mean(v_th[w]),
                     ratio_measured_theta=np.mean(U[w]) / np.mean(v_thm[w]),
                     theta_meas_end=thm[-1],
                     beta_eff=float(np.median(d0 * chi / U[w] - K * ell[w])),
                     beta0=beta0, K_l_mid=float(K * np.median(ell[w])))
            # fixed l: first crossing of --at-l after day 1
            k = np.flatnonzero((tg > DAY) & (np.diff(np.r_[ell[0], ell] > a.at_l) != 0))
            if k.size:
                i = k[0]
                r.update(U_at_l_nm_d=U[i] * 1e9 * DAY, v_at_l_nm_d=v_th[i] * 1e9 * DAY,
                         t_at_l_d=tg[i] / DAY)
        r.update(l_start_um=ell[0] * 1e6, l_end_um=ell[-1] * 1e6)
        rows.append(r)
        curves.append((name, tg, U, v_th, m))
        print(f"{name}: theta {m['theta']:.1f} H {m['H']*1e6:.0f} um mode {m['mode']} "
              f"d0 {d0:.5e} K {K:.4e} beta0 {beta0:.4e} | A {Ag[0]*1e12:.0f} -> {Ag[-1]*1e12:.0f} um2 | "
              f"U {r['U_mean_nm_d']:+.2f} nm/d vs theory {r['v_young_nm_d']:+.2f} "
              f"(ratio {r.get('ratio_young', float('nan')):.4f}; with measured theta "
              f"{r.get('ratio_measured_theta', float('nan')):.4f}); beta_eff {r.get('beta_eff', float('nan')):.3e}")

    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(a.out, "type3_speeds.csv"), "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=keys)
        wtr.writeheader()
        wtr.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.0, 6.0), sharex=True, constrained_layout=True)
    styles = {"both": "-", "x0": ":"}
    for name, tg, U, v_th, m in curves:
        ls = styles[m["mode"]]
        ln, = ax.plot(tg / DAY, U * 1e9 * DAY, lw=1.6 if ls == "-" else 2.4, ls=ls, label=name)
        ax.plot(tg / DAY, v_th * 1e9 * DAY, lw=0.9, ls="--", color=ln.get_color())
        if abs(m["cos"]) > 1e-9:
            ax2.plot(tg[1:] / DAY, U[1:] / v_th[1:], lw=1.6 if ls == "-" else 2.4, ls=ls, color=ln.get_color())
    for a_ in (ax, ax2):
        a_.axvline(a.t0, color="0.6", lw=0.8)
    ax2.axhline(1.0, color="k", lw=0.6)
    ax2.set_ylim(0.9, 1.2)
    ax2.set_xlabel("time [d]")
    ax.set_ylabel("face speed U [nm/day]")
    ax2.set_ylabel("measured / series law")
    ax.set_title("thin: series law with Young θ; dotted: x0 supply only (3D), on top of 3A", fontsize=8)
    fig.legend(loc="outside right center", fontsize=7, frameon=False)
    bad = check_figure(fig)
    fig.savefig(os.path.join(a.out, "type3_speeds.png"), dpi=200)
    print("figcheck:", bad or "OK")
    print("wrote", os.path.join(a.out, "type3_speeds.csv"), "and type3_speeds.png")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
