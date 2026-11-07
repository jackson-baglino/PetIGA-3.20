"""Plan type 4 (lens in a wedge): speed of each face against the series law.

Geometry. x0 is the narrow end, x1 the wide end; both are vapor faces at
saturation (h = 1). The walls are straight lines y_bot(x), y_top(x) from the
geometry opts (-wall_bot_y0/-slope, -wall_top_y0/-slope; a flat channel has
none and the width is -Ly), so the width is w(x) = w0 + 2 m x with
m = tan(omega).

Per face and snapshot (ice_bodies.csv and ice_bodies_contacts.csv):

  x_c     contact-line x (mean of the y0 and y1 contacts of that face)
  w_c     chord width w(x_c)
  s       sagitta: axis crossing (face_lo / face_hi) minus x_c, signed so that
          s > 0 means the face bulges into the ice (concave, it grows)
  x_eq    area-equivalent flat position: x_c moved by A_seg / w_c into the ice
          (lo face: +x, hi face: -x), A_seg the circular segment of chord w_c
          and sagitta |s| (sign of s). Nothing about theta is assumed: the
          shape is measured. For a flat pore this is the same face as the
          type 3 area-based one (checked here: the two faces' rates add up to
          dA/dt from SSA_evo.dat, printed as "area check").
  U       growth speed of the face: -dx_eq/dt (lo), +dx_eq/dt (hi)

Theory (the type 3 series law with the wedge put in; no fitted parameter):

  v = d0 chi / (beta0 + K l_eff),   K = rho_ice / (D_v rho_vs)
  chi    = 2 cos(theta_eff) / w_c,  theta_eff = theta - omega (narrow face),
           theta + omega (wide face): the wall is tilted by omega relative
           to the face normal. At theta = 90 deg this is 1/r, r the distance
           from the apex: the exact apex-centred arc (4B-1).
  l_eff  = w_eq * integral(dx / w(x)) from the face to its end
         = (w_eq / 2m) ln(w_eq / w_end) for a taper, l for a flat pore.
           Steady 1D diffusion through a channel whose width changes, per
           unit face width. Kinetics on the arc length instead of the chord
           are neglected (arc/chord is 1.05 at 60 deg; beta0 is ~7 % of the
           resistance, so ~0.3 % in v).

A second prediction uses the MEASURED curvature 1/R of the segment
(R = (w_c^2/4 + s^2) / (2 s)). Means over --t0 .. the window end. The
window ends 2 days before the face's contact line comes within --end-gap of
its end face (ice then touches a vapor face and the run is no longer a lens
in a pore; the solver's resolution.dat shows the reservoir-face phi rising).

    venv_cryosuction/bin/python studies/cryosuction/scripts/analyze_type4.py \
        --batch <batch folder> --out <folder> [--t0 40] [--end-gap 5e-6]

Prints every input. Writes <out>/type4_faces.csv, <out>/type4_series.csv and
<out>/type4_faces.png (checked with figcheck).
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


def geom_opts(run):
    """Wall lines, theta, T and beta0 from the run's staged opts."""
    o = pplib.read_opts(run)
    f = lambda k, d=float("nan"): pplib.opt_float(o, k, d)  # noqa: E731
    Lx, Ly = f("-Lx"), f("-Ly")
    if math.isfinite(f("-wall_bot_slope")):
        b0, bs, t0, ts = f("-wall_bot_y0"), f("-wall_bot_slope"), f("-wall_top_y0"), f("-wall_top_slope")
    else:
        b0, bs, t0, ts = 0.0, 0.0, Ly, 0.0
    w0, dw = t0 - b0, ts - bs                    # w(x) = w0 + dw x, dw = 2 m
    cos_th = (f("-gamma_as") - f("-gamma_is")) / f("-gamma_ia")
    return dict(Lx=Lx, w0=w0, dw=dw, omega=math.degrees(math.atan(dw / 2)),
                theta=math.degrees(math.acos(cos_th)), T=f("-temp"), beta0=abs(f("-beta_sub0")))


def run_id(g):
    """RUNS.md / type4_tests.txt ID from omega and theta."""
    om, th = round(g["omega"], 1), round(g["theta"])
    if th == 90:
        return "4B-1"
    if th == 85:
        return "4B-2"
    return {0.0: "4A-0", 2.5: "4A-1", 5.0: "4A-2", 10.0: "4A-3", 20.0: "4A-4"}.get(om, "4?")


def width(g, x):
    return g["w0"] + g["dw"] * x


def l_eff(g, x, side):
    w_f = width(g, x)
    x_end = 0.0 if side == "lo" else g["Lx"]
    if abs(g["dw"]) < 1e-12:
        return np.abs(x - x_end)
    w_end = width(g, x_end)
    return w_f / g["dw"] * np.abs(np.log(w_f / w_end))


def segment_area(w, s):
    """Signed area of the circular segment with chord w and sagitta s."""
    s_abs = np.abs(s)
    with np.errstate(divide="ignore", invalid="ignore"):
        R = (w ** 2 / 4 + s_abs ** 2) / (2 * s_abs)
        half = np.arcsin(np.clip(w / (2 * R), -1, 1))
        a = R ** 2 * (half - np.sin(half) * np.cos(half))
        # a sagitta larger than the half chord is more than a half disc
        a = np.where(s_abs > w / 2, np.pi * R ** 2 - a, a)
    return np.where(s_abs < 1e-12, 0.0, np.sign(s) * a), R


def faces(run, g):
    import pandas as pd
    b = pd.read_csv(os.path.join(run, "ice_bodies.csv"))
    c = pd.read_csv(os.path.join(run, "ice_bodies_contacts.csv"))
    b = b.sort_values("area_m2").drop_duplicates("t_s", keep="last").set_index("t_s").sort_index()  # largest body
    out = {}
    for side in ("lo", "hi"):
        rows = []
        for t, grp in c.groupby("t_s"):
            if t not in b.index:
                continue
            mid = b.loc[t, "com_x"]
            pts = grp[grp.x < mid] if side == "lo" else grp[grp.x > mid]
            ax = b.loc[t, f"face_{side}"]
            if len(pts) < 2 or not np.isfinite(ax):
                rows.append((t, np.nan, np.nan, np.nan, np.nan, np.nan))
                continue
            xc = pts.x.mean()
            wc = width(g, xc)
            s = (ax - xc) if side == "lo" else (xc - ax)          # > 0: bulges into the ice
            A, R = segment_area(wc, s)
            xeq = xc + A / wc if side == "lo" else xc - A / wc
            rows.append((t, xc, wc, s, xeq, R))
        out[side] = np.array(rows, dtype=float)
    return b, out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--batch", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--t0", type=float, default=40.0, help="start of the averaging window [d]")
    p.add_argument("--end-gap", type=float, default=5e-6, help="window ends 2 d before the contact line is this close to its end [m]")
    a = p.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    runs = sorted(d for d in glob.glob(os.path.join(a.batch, "*__*"))
                  if os.path.isfile(os.path.join(d, "ice_bodies_contacts.csv")))
    rows, series, curves = [], [], []
    for run in runs:
        g = geom_opts(run)
        d0 = capillary_length_m(g["T"])
        K = pplib.RHO_ICE / (dv_m2_per_s(g["T"]) * rho_vs_sat_kg_per_m3(g["T"]))
        beta0 = g["beta0"]
        name = os.path.basename(run).split("_gap")[0].split("_L270")[0].replace("_2D", "") \
            + f" θ{g['theta']:.0f}"
        name = f"{run_id(g)} {name}"
        b, F = faces(run, g)
        print(f"\n{name}: omega {g['omega']:.3f} deg, w(x) = {g['w0']*1e6:.2f} um + {g['dw']:.5f} x, "
              f"Lx {g['Lx']*1e6:.1f} um, theta {g['theta']:.2f}, d0 {d0:.5e}, K {K:.4e}, beta0 {beta0:.4e}")
        # area check: the two faces' rates against dA/dt
        ssa = pplib.load_ssa(run)
        o = np.argsort(ssa[:, pplib.TIME])
        tA, A = ssa[o, pplib.TIME], ssa[o, pplib.ICE]
        for side in ("lo", "hi"):
            arr = F[side]
            t, xc, wc, s, xeq, R = arr.T
            sign = -1.0 if side == "lo" else 1.0
            th_eff = g["theta"] - g["omega"] if side == "lo" else g["theta"] + g["omega"]
            chi = 2 * math.cos(math.radians(th_eff)) / wc
            chi_m = np.where(np.isfinite(R) & (R > 0), np.sign(s) / R, 0.0)
            le = l_eff(g, xeq, side)
            ok = np.isfinite(xeq)
            # window end: 2 d before the contact line is within end-gap of its end
            dist = xc if side == "lo" else g["Lx"] - xc
            hit = np.flatnonzero(~ok | (dist < a.end_gap))
            t_hit = t[hit[0]] if hit.size else np.inf
            t_end = min(t_hit - 2 * DAY, t[ok][-1] if ok.any() else 0) if np.isfinite(t_hit) else t[ok][-1]
            # 5-day centered speed (daily snapshots) for the figure and series file
            U = np.full_like(t, np.nan)
            for i in range(len(t)):
                j0, j1 = i - 2, i + 2
                if j0 >= 0 and j1 < len(t) and ok[j0] and ok[j1]:
                    U[i] = sign * (xeq[j1] - xeq[j0]) / (t[j1] - t[j0])
            v_th = d0 * chi / (beta0 + K * le)
            v_thm = d0 * chi_m / (beta0 + K * le)
            for i in range(len(t)):
                series.append(dict(run=name, face=side, t_d=t[i] / DAY, x_contact_um=xc[i] * 1e6,
                                   chord_um=wc[i] * 1e6, sagitta_um=s[i] * 1e6, x_eq_um=xeq[i] * 1e6,
                                   l_eff_um=le[i] * 1e6, U_nm_d=U[i] * 1e9 * DAY,
                                   law_young_nm_d=v_th[i] * 1e9 * DAY, law_measured_chi_nm_d=v_thm[i] * 1e9 * DAY))
            w = ok & (t >= a.t0 * DAY) & (t <= t_end)
            r = dict(run=name, face=side, omega=g["omega"], theta=g["theta"], theta_eff=th_eff,
                     window_d=f"{a.t0:g}-{t_end/DAY:.0f}" if w.sum() > 1 else "none",
                     t_contact_d=t_hit / DAY if np.isfinite(t_hit) else np.nan,
                     x_eq_start_um=xeq[0] * 1e6)
            if w.sum() > 1:
                i0, i1 = np.flatnonzero(w)[[0, -1]]
                Um = sign * (xeq[i1] - xeq[i0]) / (t[i1] - t[i0])
                vm = np.trapezoid(v_th[w], t[w]) / (t[i1] - t[i0])
                vmm = np.trapezoid(v_thm[w], t[w]) / (t[i1] - t[i0])
                # acceleration: slope of a straight-line fit to U in the window
                uu = np.isfinite(U) & w
                acc = np.polyfit(t[uu] / DAY, U[uu] * 1e9 * DAY, 1)[0] if uu.sum() > 3 else np.nan
                # kinetic resistance the run implies (type 3 found ~8.8e5 s/m), and the
                # share of the law's resistance that is diffusion
                with np.errstate(divide="ignore", invalid="ignore"):
                    be = d0 * chi[uu] / U[uu] - K * le[uu]
                r.update(beta_eff=float(np.median(be)) if uu.any() else np.nan,
                         diffusion_share=float(np.median(K * le[w] / (beta0 + K * le[w]))))
                # the same ratio over the last 30 days of the window (slow relaxers)
                wl = w & (t >= t[i1] - 30 * DAY)
                j0, j1 = np.flatnonzero(wl)[[0, -1]]
                vl = np.trapezoid(v_th[wl], t[wl]) / (t[j1] - t[j0])
                r["ratio_young_last30d"] = sign * (xeq[j1] - xeq[j0]) / (t[j1] - t[j0]) / vl if vl else np.nan
                r.update(U_mean_nm_d=Um * 1e9 * DAY, law_young_nm_d=vm * 1e9 * DAY,
                         ratio_young=Um / vm if vm else np.nan,
                         ratio_measured_chi=Um / vmm if vmm else np.nan,
                         accel_nm_d2=acc, chord_um_window=f"{wc[i0]*1e6:.1f}-{wc[i1]*1e6:.1f}",
                         l_eff_um_window=f"{le[i0]*1e6:.1f}-{le[i1]*1e6:.1f}",
                         sagitta_um_end=s[i1] * 1e6, R_meas_um_end=R[i1] * 1e6,
                         R_young_um=wc[i1] / (2 * math.cos(math.radians(th_eff))) * 1e6
                         if abs(math.cos(math.radians(th_eff))) > 1e-9 else np.inf)
            rows.append(r)
            curves.append((name, side, t, U, v_th, g))
            print(f"  {side}: theta_eff {th_eff:.1f}, x_eq {xeq[0]*1e6:.2f} -> {np.nanmax(xeq[ok]) * 1e6 if side == 'hi' else np.nanmin(xeq[ok]) * 1e6:.2f} um, "
                  f"contact within {a.end_gap*1e6:g} um of its end at {r['t_contact_d']:.1f} d, window {r['window_d']} d | "
                  f"U {r.get('U_mean_nm_d', np.nan):+.2f} nm/d, law {r.get('law_young_nm_d', np.nan):+.2f} "
                  f"(ratio {r.get('ratio_young', np.nan):.4f}, last 30 d {r.get('ratio_young_last30d', np.nan):.4f}; measured chi {r.get('ratio_measured_chi', np.nan):.4f}); "
                  f"dU/dt {r.get('accel_nm_d2', np.nan):+.4f} nm/d^2; beta_eff {r.get('beta_eff', np.nan):.3e}, "
                  f"diffusion share {r.get('diffusion_share', np.nan):.3f}")
        # area check (both faces exchanging; uses the time range where both exist)
        lo, hi = F["lo"], F["hi"]
        okb = np.isfinite(lo[:, 4]) & np.isfinite(hi[:, 4]) & (lo[:, 0] >= a.t0 * DAY)
        if okb.sum() > 2:
            i0, i1 = np.flatnonzero(okb)[[0, -1]]
            t0_, t1_ = lo[i0, 0], lo[i1, 0]
            dA = np.interp(t1_, tA, A) - np.interp(t0_, tA, A)
            # area between the two x_eq positions in the wedge, start and end
            def band(x1, x2, g=g):
                return g["w0"] * (x2 - x1) + g["dw"] / 2 * (x2 ** 2 - x1 ** 2)
            dA_faces = band(lo[i1, 4], hi[i1, 4]) - band(lo[i0, 4], hi[i0, 4])
            print(f"  area check {t0_/DAY:.0f}-{t1_/DAY:.0f} d: dA (SSA_evo) {dA*1e12:.1f} um^2, "
                  f"from x_eq faces {dA_faces*1e12:.1f} um^2 (ratio {dA_faces/dA if dA else np.nan:.4f})")
            for r in rows[-2:]:
                r["area_check"] = dA_faces / dA if dA else np.nan

    def dump(path, recs):
        keys = list(dict.fromkeys(k for r in recs for k in r))
        with open(path, "w", newline="") as fh:
            wtr = csv.DictWriter(fh, fieldnames=keys)
            wtr.writeheader()
            wtr.writerows(recs)
    dump(os.path.join(a.out, "type4_faces.csv"), rows)
    dump(os.path.join(a.out, "type4_series.csv"), series)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure
    fig, axs = plt.subplots(2, 2, figsize=(7.5, 5.6), sharex=True, constrained_layout=True)
    for name, side, t, U, v_th, _g in curves:
        col = 0 if side == "lo" else 1
        ln, = axs[0, col].plot(t / DAY, U * 1e9 * DAY, lw=1.4, label=name if side == "lo" else None)
        axs[0, col].plot(t / DAY, v_th * 1e9 * DAY, lw=0.8, ls="--", color=ln.get_color())
        with np.errstate(divide="ignore", invalid="ignore"):
            rat = np.where(np.abs(v_th) > 1e-14, U / v_th, np.nan)
        axs[1, col].plot(t / DAY, rat, lw=1.4, color=ln.get_color())
    for col, ttl in ((0, "narrow-side face (toward x = 0)"), (1, "wide-side face (toward x = Lx)")):
        axs[0, col].set_title(ttl, fontsize=8)
        axs[1, col].axhline(1.0, color="k", lw=0.6)
        axs[1, col].set_ylim(0.6, 1.4)
        axs[1, col].set_xlabel("time [d]")
        for r_ in (0, 1):
            axs[r_, col].axvline(a.t0, color="0.6", lw=0.8)
    axs[0, 0].set_ylabel("growth speed U [nm/day]")
    axs[1, 0].set_ylabel("measured / series law")
    fig.legend(loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    bad = check_figure(fig)
    fig.savefig(os.path.join(a.out, "type4_faces.png"), dpi=200)
    print("\nfigcheck:", bad or "OK")
    print("wrote type4_faces.csv, type4_series.csv, type4_faces.png in", a.out)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
