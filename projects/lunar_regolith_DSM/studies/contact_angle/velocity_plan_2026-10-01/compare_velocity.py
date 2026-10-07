#!/usr/bin/env python3
"""Measured interface velocities of the velocity study against the theory.

Reads the per-run CSVs that postprocess/meniscus_velocity.py (channel) and
postprocess/wedge_gt_velocity.py (wedge) leave in each run folder, for every
batch{A,B,C,D} run under --root, and writes into --out:

    velocity_summary.csv          one row per run (and per meniscus on the wedge)
    fig_channel_summary.png       v against theta, v against sigma_inf
    fig_channel_time_maps.png     v(t) maps: time across, theta or sigma_inf up
    fig_channel_phase_diagram.png theta x sigma_inf, coloured by v
    fig_wedge_summary.png         as the channel summary, both menisci
    fig_wedge_time_maps_ac*.png   as the channel maps, both menisci

WHAT "VELOCITY" MEANS HERE

  channel  U = (dA/dt) / (2H): the rate at which the mean meniscus position
           advances. It is the quantity the series-resistance relation
           predicts (vapour flux through the channel cross-section H), and it
           is insensitive to the meniscus still relaxing from the clipped-disc
           IC to its equilibrium arc, which the mid-plane and contact-line
           velocities are not. (dA/dt)/arc-length, the mean NORMAL velocity,
           is smaller by arc/2H = 1.03 to 1.24.
  wedge    the centreline velocity of each meniscus, from wedge_gt_velocity.py.
           The band IC is the theta = 90 shape, so at any other angle the
           centreline also carries the shape relaxation; only its late-time
           value and its slope against sigma_inf are comparable to theory.

THEORY (make_theory_figures.py), with beta = beta_sub0, no fitted offset, and
evaluated at each run's MEASURED meniscus position:

    channel   v = (sigma_inf + 2 d0 cos(theta)/H) / (beta + K l),  l = (Lx - A/H)/2
    wedge     v_in/out = (sigma_inf + d0 (cos(theta) +/- sin(alpha)) / (r sin(alpha)))
                         / (beta + K r ln(r/r_L or r_R/r))

v > 0 is growth throughout.

Usage:
    python3 compare_velocity.py --root <velocity_study folder> [--out <dir>]
"""
import argparse
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import make_theory_figures as T    # constants, theory curvatures and the house style

LX = 375e-6                        # channel_2D_H125um_eps0.86um
ALPHAS = (1e-3, 1e-2)
NM_DAY = T.NM_DAY
DAY = T.DAY

# growth = blue (ice), sublimation = orange, neutral grey at zero
DIVERGING = LinearSegmentedColormap.from_list(
    "sublime_grow", ["#8a3210", "#eb6834", "#f6c9b3", "#e6e9ed",
                     "#b5d2f2", "#2a78d6", "#0f3462"])
DIVERGING.set_bad("#ffffff")
C_A = {1e-3: T.BLUE, 1e-2: T.ORANGE}           # the two alpha_c series
C_M = {"inner": T.BLUE, "outer": T.ORANGE}     # the two wedge menisci
AC_LABEL = {1e-3: r"$\alpha_c = 10^{-3}$", 1e-2: r"$\alpha_c = 10^{-2}$"}


def beta0(ac):
    return T.BETA_HK1 / ac


def v_channel(th_deg, sig, ell, ac):
    return (sig + 2.0 * T.D0 * np.cos(np.radians(th_deg)) / T.H) / (beta0(ac) + T.K * ell)


def v_wedge(which, th_deg, sig, r, ac):
    c = np.cos(np.radians(th_deg))
    if which == "inner":
        return (sig + T.D0 * (c + T.SA) / (r * T.SA)) / (beta0(ac) + T.K * r * np.log(r / T.R_L))
    return (sig + T.D0 * (c - T.SA) / (r * T.SA)) / (beta0(ac) + T.K * r * np.log(T.R_R / r))


# --- loading -----------------------------------------------------------------
NAME = re.compile(r"(channel|wedge).*theta(\d+)(?:_sig([mp])(\d)e-5)?_ac(1e-[23])$")


def load(root):
    """Every run under root/batch_*/ as a dict with its time series."""
    runs = []
    for d in sorted(glob.glob(os.path.join(root, "batch_*", "*__*"))):
        m = NAME.search(d)
        if not m:
            continue
        geom, th, s, n, ac = m.groups()
        sig = 0.0 if s is None else (1 if s == "p" else -1) * int(n) * 1e-5
        r = dict(geom=geom, theta=int(th), sigma=sig, ac=float(ac), dir=d)
        if geom == "channel":
            f = os.path.join(d, "meniscus_velocity.csv")
            if not os.path.exists(f):
                print("  missing", f)
                continue
            c = np.genfromtxt(f, delimiter=",", skip_header=4, names=True)
            r["t"] = c["time_s"]
            r["ell"] = 0.5 * (LX - c["area_m2"] / T.H)
            r["v"] = {"mean": c["dA_dt_m2_s"] / (2.0 * T.H)}
            r["v_th"] = {"mean": v_channel(r["theta"], sig, r["ell"], r["ac"])}
        else:
            f = os.path.join(d, "wedge_gt_velocity.csv")
            if not os.path.exists(f):
                print("  missing", f)
                continue
            c = np.genfromtxt(f, delimiter=",", names=True)
            ok = np.isfinite(c["vn_left_meas"]) & np.isfinite(c["vn_right_meas"])
            c = c[ok]
            r["t"] = c["time"]
            r["r"] = {"inner": c["r_left"], "outer": c["r_right"]}
            r["v"] = {"inner": c["vn_left_meas"], "outer": c["vn_right_meas"]}
            r["v_th"] = {k: v_wedge(k, r["theta"], sig, r["r"][k], r["ac"])
                         for k in ("inner", "outer")}
        runs.append(r)
    return runs


def window(r, lo, hi):
    """Mask for lo..hi as fractions of the run's own span."""
    t = r["t"]
    return (t >= lo * t.max()) & (t <= hi * t.max()) & (t > DAY)


def avg(r, key, lo=0.5, hi=1.0, theory=False):
    w = window(r, lo, hi)
    return float(np.mean((r["v_th"] if theory else r["v"])[key][w]))


def pick(runs, geom, ac, sweep):
    """The theta sweep (sigma = 0) or the sigma sweep (theta = 60), in order."""
    sel = [r for r in runs if r["geom"] == geom and r["ac"] == ac]
    if sweep == "theta":
        return sorted((r for r in sel if r["sigma"] == 0.0), key=lambda r: r["theta"])
    return sorted((r for r in sel if r["theta"] == 60), key=lambda r: r["sigma"])


def fit_beta(runs, ac):
    """Least-squares beta in v = F/(beta + K l) over the channel runs."""
    sel = [r for r in runs if r["geom"] == "channel" and r["ac"] == ac]
    F = np.array([r["sigma"] + 2 * T.D0 * np.cos(np.radians(r["theta"])) / T.H for r in sel])
    V = np.array([avg(r, "mean") for r in sel])
    L = np.array([float(np.mean(r["ell"][window(r, 0.5, 1.0)])) for r in sel])
    bs = np.linspace(0.2, 5.0, 48001) * beta0(ac)
    err = [np.sum((V - F / (b + T.K * L)) ** 2) for b in bs]
    return bs[int(np.argmin(err))]


# --- figures -----------------------------------------------------------------
def save(fig, out, name):
    fig.savefig(os.path.join(out, name), dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", name)


def sweep_axis(ax, sweep):
    if sweep == "theta":
        ax.set_xlabel(r"contact angle $\theta$ [deg]")
        ax.set_xticks([30, 60, 90, 120, 150])
    else:
        ax.set_xlabel(r"reservoir supersaturation $\sigma_\infty$ [$10^{-5}$]")
    ax.axhline(0, color=T.MUTED, lw=0.8, zorder=1)


def xval(r, sweep):
    return r["theta"] if sweep == "theta" else r["sigma"] * 1e5


def fig_channel_summary(runs, out):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), constrained_layout=True)
    for ax, sweep, letter in zip(axes, ("theta", "sigma"), "ab"):
        for ac in ALPHAS:
            rr = pick(runs, "channel", ac, sweep)
            if not rr:
                continue
            x = np.array([xval(r, sweep) for r in rr])
            ax.plot(x, [avg(r, "mean", theory=True) * NM_DAY for r in rr],
                    color=C_A[ac], lw=2.0, zorder=2,
                    label="theory, " + AC_LABEL[ac])
            ax.plot(x, [avg(r, "mean") * NM_DAY for r in rr], "o", ms=9,
                    mfc=C_A[ac], mec="white", mew=2.0, zorder=3,
                    label="simulation, " + AC_LABEL[ac])
        sweep_axis(ax, sweep)
        ax.set_ylabel("meniscus velocity [nm/day]")
        T.panel(ax, letter)
    axes[0].legend(fontsize=12, loc="upper right")
    axes[0].text(0.03, 0.05, r"$\sigma_\infty = 0$", transform=axes[0].transAxes, color=T.MUTED)
    axes[1].text(0.03, 0.92, r"$\theta = 60^\circ$", transform=axes[1].transAxes, color=T.MUTED)
    save(fig, out, "fig_channel_summary.png")


def time_map(ax, rr, key, sweep, t_end):
    """Rows = runs, columns = time, colour = velocity. Returns the mappable."""
    tg = np.linspace(1.0, t_end, 240)                  # days
    Z = np.full((len(rr), tg.size), np.nan)
    for i, r in enumerate(rr):
        td = r["t"] / DAY
        m = td >= 1.0
        Z[i] = np.interp(tg, td[m], r["v"][key][m] * NM_DAY, left=np.nan, right=np.nan)
    vmax = np.nanpercentile(np.abs(Z), 99)
    im = None
    for i in range(len(rr)):                           # one band per run, gap between
        im = ax.imshow(Z[i:i + 1], aspect="auto", cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                       extent=(tg[0], tg[-1], i - 0.44, i + 0.44), interpolation="nearest")
    ax.set_ylim(-0.6, len(rr) - 0.4)
    ax.set_yticks(range(len(rr)))
    if sweep == "theta":
        ax.set_yticklabels(["%d°" % r["theta"] for r in rr])
        ax.set_ylabel(r"contact angle $\theta$")
    else:
        ax.set_yticklabels(["%+d" % round(r["sigma"] * 1e5) if r["sigma"] else "0" for r in rr])
        ax.set_ylabel(r"$\sigma_\infty$ [$10^{-5}$]")
    ax.set_xlim(tg[0], tg[-1])
    ax.grid(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    return im


def fig_channel_time_maps(runs, out):
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 8.6), constrained_layout=True)
    for i, ac in enumerate(ALPHAS):
        for j, sweep in enumerate(("theta", "sigma")):
            ax = axes[i, j]
            rr = pick(runs, "channel", ac, sweep)
            if not rr:
                ax.set_axis_off()
                continue
            im = time_map(ax, rr, "mean", sweep, 90.0)
            cb = fig.colorbar(im, ax=ax, pad=0.02)
            cb.set_label("velocity [nm/day]", fontsize=13)
            cb.outline.set_visible(False)
            ax.set_title(AC_LABEL[ac] + (r",  $\sigma_\infty = 0$" if sweep == "theta"
                                         else r",  $\theta = 60^\circ$"),
                         fontsize=14, loc="left", color=T.MUTED)
            if i == 1:
                ax.set_xlabel("time [days]")
            T.panel(ax, "abcd"[2 * i + j])
    save(fig, out, "fig_channel_time_maps.png")


def fig_channel_phase(runs, out):
    """theta x sigma_inf coloured by velocity: theory field, simulations on top."""
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.0), constrained_layout=True)
    th = np.linspace(20, 160, 281)
    sg = np.linspace(-3.6e-5, 3.6e-5, 241)
    TH, SG = np.meshgrid(th, sg)
    for ax, ac, letter in zip(axes, ALPHAS, "ab"):
        Z = v_channel(TH, SG, T.L_CH, ac) * NM_DAY
        vmax = np.abs(Z).max()
        im = ax.pcolormesh(TH, SG * 1e5, Z, cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                           shading="auto", rasterized=True)
        ax.contour(TH, SG * 1e5, Z, levels=[0.0], colors=[T.INK], linewidths=1.6)
        sel = [r for r in runs if r["geom"] == "channel" and r["ac"] == ac]
        if sel:
            ax.scatter([r["theta"] for r in sel], [r["sigma"] * 1e5 for r in sel],
                       c=[avg(r, "mean", 0.05, 0.25) * NM_DAY for r in sel],
                       cmap=DIVERGING, vmin=-vmax, vmax=vmax, s=190,
                       edgecolors="white", linewidths=2.2, zorder=3)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("velocity [nm/day]", fontsize=13)
        cb.outline.set_visible(False)
        ax.set_xlabel(r"contact angle $\theta$ [deg]")
        ax.set_ylabel(r"reservoir supersaturation $\sigma_\infty$ [$10^{-5}$]")
        ax.set_xticks([30, 60, 90, 120, 150])
        ax.grid(False)
        ax.set_title(AC_LABEL[ac], fontsize=14, loc="left", color=T.MUTED)
        ax.text(0.03, 0.93, "growth", transform=ax.transAxes, color="white")
        ax.text(0.97, 0.04, "sublimation", transform=ax.transAxes, ha="right", color="white")
        T.panel(ax, letter)
    save(fig, out, "fig_channel_phase_diagram.png")


def fig_wedge_summary(runs, out):
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0), constrained_layout=True)
    for i, ac in enumerate(ALPHAS):
        for j, sweep in enumerate(("theta", "sigma")):
            ax = axes[i, j]
            rr = pick(runs, "wedge", ac, sweep)
            x = np.array([xval(r, sweep) for r in rr])
            for key in ("inner", "outer"):
                ax.plot(x, [avg(r, key, theory=True) * NM_DAY for r in rr],
                        color=C_M[key], lw=2.0, zorder=2, label="theory, %s meniscus" % key)
                ax.plot(x, [avg(r, key) * NM_DAY for r in rr], "o", ms=9,
                        mfc=C_M[key], mec="white", mew=2.0, zorder=3,
                        label="simulation, %s meniscus" % key)
            sweep_axis(ax, sweep)
            ax.set_ylabel("centreline velocity [nm/day]")
            ax.set_title(AC_LABEL[ac] + (r",  $\sigma_\infty = 0$" if sweep == "theta"
                                         else r",  $\theta = 60^\circ$"),
                         fontsize=14, loc="left", color=T.MUTED)
            T.panel(ax, "abcd"[2 * i + j])
    axes[0, 0].legend(fontsize=12)
    save(fig, out, "fig_wedge_summary.png")


def fig_wedge_time_maps(runs, out, ac):
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 8.6), constrained_layout=True)
    for i, key in enumerate(("inner", "outer")):
        for j, sweep in enumerate(("theta", "sigma")):
            ax = axes[i, j]
            rr = pick(runs, "wedge", ac, sweep)
            if not rr:
                ax.set_axis_off()
                continue
            im = time_map(ax, rr, key, sweep, 150.0)
            cb = fig.colorbar(im, ax=ax, pad=0.02)
            cb.set_label("velocity [nm/day]", fontsize=13)
            cb.outline.set_visible(False)
            ax.set_title("%s meniscus, %s" % (key, r"$\sigma_\infty = 0$" if sweep == "theta"
                                              else r"$\theta = 60^\circ$"),
                         fontsize=14, loc="left", color=T.MUTED)
            if i == 1:
                ax.set_xlabel("time [days]")
            T.panel(ax, "abcd"[2 * i + j])
    save(fig, out, "fig_wedge_time_maps_ac%s.png" % ("1e-3" if ac == 1e-3 else "1e-2"))


def write_summary(runs, out):
    path = os.path.join(out, "velocity_summary.csv")
    with open(path, "w") as fh:
        fh.write("# second half of each run; v > 0 = growth; theory uses beta_sub0 at the "
                 "measured meniscus position\n")
        fh.write("geometry,alpha_c,theta_deg,sigma_inf,meniscus,t_end_d,"
                 "v_meas_m_s,v_theory_m_s,v_meas_nm_day,v_theory_nm_day,meas_over_theory\n")
        for r in sorted(runs, key=lambda r: (r["geom"], r["ac"], r["sigma"] != 0, r["theta"], r["sigma"])):
            for key in r["v"]:
                vm, vt = avg(r, key), avg(r, key, theory=True)
                ratio = vm / vt if abs(vt) > 1e-15 else float("nan")
                fh.write("%s,%g,%d,%g,%s,%.1f,%.4e,%.4e,%.2f,%.2f,%.3f\n"
                         % (r["geom"], r["ac"], r["theta"], r["sigma"], key,
                            r["t"].max() / DAY, vm, vt, vm * NM_DAY, vt * NM_DAY, ratio))
    print("wrote velocity_summary.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="folder holding the batch_* folders")
    ap.add_argument("--out", default=None, help="output folder (default: <root>/figures)")
    args = ap.parse_args()
    out = args.out or os.path.join(args.root, "figures")
    os.makedirs(out, exist_ok=True)

    runs = load(args.root)
    print("%d runs: %d channel, %d wedge" % (
        len(runs), sum(r["geom"] == "channel" for r in runs),
        sum(r["geom"] == "wedge" for r in runs)))
    for ac in ALPHAS:
        b = fit_beta(runs, ac)
        print("  alpha_c = %g: fitted beta = %.3e s/m = %.3f x beta_sub0" % (ac, b, b / beta0(ac)))

    write_summary(runs, out)
    fig_channel_summary(runs, out)
    fig_channel_time_maps(runs, out)
    fig_channel_phase(runs, out)
    fig_wedge_summary(runs, out)
    for ac in ALPHAS:
        fig_wedge_time_maps(runs, out, ac)


if __name__ == "__main__":
    main()
