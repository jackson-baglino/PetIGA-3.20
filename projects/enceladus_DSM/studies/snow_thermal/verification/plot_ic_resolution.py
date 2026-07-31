#!/usr/bin/env python3
"""Plot the layered-benchmark IC resolution study.

Reads ic_resolution_study.csv (written by run_ic_resolution_study.sh) and
produces a four-panel figure characterising how much of the k_eff benchmark's
error belongs to the initial condition rather than to the cell-problem solver.

Usage:
    python3 plot_ic_resolution.py [--theme light|dark|both] [--csv PATH]
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Categorical slots 1-3 of the validated reference palette. Fixed order, never
# cycled. Checked with postprocess/validate_palette.py: worst pair dE = 10.0
# under protanopia, above the 8.0 threshold.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]

K_ICE, K_AIR = 2.29, 0.02          # W/m/K, as hardcoded in material_properties.c
PHI_EXACT = 0.5                    # the continuous profile is exactly 50% ice

THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e",
                  muted="#8a8985", grid="#e3e2dd"),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7",
                  muted="#8a8985", grid="#333330"),
}


def load(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            rows.append({k: (v if k == "sweep" else float(v)) for k, v in r.items()})
    return rows


def k_parallel(phi):
    """Arithmetic mean — conduction along the layers."""
    return phi * K_ICE + (1.0 - phi) * K_AIR


def k_perp(phi):
    """Harmonic mean — conduction across the layers."""
    return 1.0 / (phi / K_ICE + (1.0 - phi) / K_AIR)


def thin_log_xticks(ax):
    """Decade majors only. Matplotlib's default log minor labels collide badly
    when the data spans well under a decade."""
    from matplotlib.ticker import LogFormatterSciNotation, LogLocator, NullFormatter
    ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=(3.0, 4.0, 6.0, 7.0, 8.0, 9.0)))
    ax.xaxis.set_major_formatter(LogFormatterSciNotation())
    ax.xaxis.set_minor_formatter(NullFormatter())


def style_axes(ax, t):
    ax.set_facecolor(t["surface"])
    ax.grid(True, which="both", color=t["grid"], linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(t["grid"])
    ax.tick_params(colors=t["ink2"], labelsize=9)
    ax.xaxis.label.set_color(t["ink2"])
    ax.yaxis.label.set_color(t["ink2"])
    ax.title.set_color(t["ink"])


def make_figure(rows, theme):
    t = THEMES[theme]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
    fig.patch.set_facecolor(t["surface"])

    err = lambda r: abs(r["phi_bar_solver"] - PHI_EXACT)

    # ---------------------------------------------------------------- panel 1
    ax = axes[0, 0]
    style_axes(ax, t)
    s = sorted([r for r in rows if r["sweep"] == "C_fixed_eps"], key=lambda r: r["dx"])
    x = np.array([r["dx"] for r in s])
    y = np.array([err(r) for r in s])
    ax.loglog(x, y, "o-", color=SERIES[0], lw=2, ms=7, zorder=3,
              label=r"measured, $\epsilon=2\times10^{-5}$ m")
    ref = y[-1] * (x / x[-1]) ** 2
    ax.loglog(x, ref, "--", color=t["muted"], lw=1.4, zorder=2, label=r"slope 2 ($\Delta x^2$)")
    ax.set_xlabel(r"element size $\Delta x$  [m]")
    ax.set_ylabel(r"$|\bar\phi - 0.5|$")
    ax.set_title("1.  Refining the mesh at fixed $\\epsilon$\nsecond-order convergence", fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"])
    thin_log_xticks(ax)

    # ---------------------------------------------------------------- panel 2
    ax = axes[0, 1]
    style_axes(ax, t)
    s = sorted([r for r in rows if r["sweep"] == "B_fixed_mesh"], key=lambda r: r["eps"])
    x = np.array([r["eps"] for r in s])
    y = np.array([err(r) for r in s])
    ax.loglog(x, y, "s-", color=SERIES[1], lw=2, ms=7, zorder=3, label="measured, $N=256$")
    ref = y[-1] * (x[-1] / x)
    ax.loglog(x, ref, "--", color=t["muted"], lw=1.4, zorder=2, label=r"slope $-1$ ($1/\epsilon$)")
    ax.set_xlabel(r"interface width $\epsilon$  [m]")
    ax.set_ylabel(r"$|\bar\phi - 0.5|$")
    ax.set_title("2.  Sharpening $\\epsilon$ at fixed mesh\nerror grows as $1/\\epsilon$", fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"])

    # ---------------------------------------------------------------- panel 3
    ax = axes[1, 0]
    style_axes(ax, t)
    labels = {"A_follow": r"A  mesh follows $\epsilon$",
              "B_fixed_mesh": "B  fixed mesh",
              "C_fixed_eps": r"C  fixed $\epsilon$"}
    marks = {"A_follow": "^", "B_fixed_mesh": "s", "C_fixed_eps": "o"}
    allx, ally = [], []
    for i, (name, lab) in enumerate(labels.items()):
        s = [r for r in rows if r["sweep"] == name]
        gx = np.array([r["dx"] ** 2 / r["eps"] for r in s])
        gy = np.array([err(r) for r in s])
        allx.extend(gx); ally.extend(gy)
        ax.loglog(gx, gy, marks[name], color=SERIES[i], ms=8, lw=0, zorder=3,
                  label=lab, markeredgecolor=t["surface"], markeredgewidth=1.2)
    allx, ally = np.array(allx), np.array(ally)
    c = np.exp(np.mean(np.log(ally) - np.log(allx)))
    xs = np.array([allx.min(), allx.max()])
    ax.loglog(xs, c * xs, "--", color=t["muted"], lw=1.4, zorder=2,
              label=rf"$|\bar\phi-0.5| = {c:.0f}\,\Delta x^2/\epsilon$")
    ax.set_xlabel(r"$\Delta x^2/\epsilon$  [m]")
    ax.set_ylabel(r"$|\bar\phi - 0.5|$")
    ax.set_title("3.  All 12 runs collapse onto one law", fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"])

    # ---------------------------------------------------------------- panel 4
    ax = axes[1, 1]
    style_axes(ax, t)
    s = sorted(rows, key=lambda r: r["elem_per_band_5_95"])
    n = np.array([r["elem_per_band_5_95"] for r in s])
    phi = np.array([r["phi_bar_solver"] for r in s])
    bias_par = np.abs(k_parallel(phi) - k_parallel(PHI_EXACT)) / k_parallel(PHI_EXACT)
    bias_per = np.abs(k_perp(phi) - k_perp(PHI_EXACT)) / k_perp(PHI_EXACT)

    # Both series, and they very nearly coincide. At phi = 0.5 the relative
    # sensitivity is identical for the two means,
    #     d ln k/d phi = 2(k_i - k_a)/(k_i + k_a) = 1.9654,
    # for ANY contrast ratio -- the arithmetic mean is (k_i+k_a)/2 with slope
    # (k_i-k_a), and the harmonic mean is 2 k_i k_a/(k_i+k_a) with slope
    # (k_i-k_a) k_perp^2/(k_i k_a), and the two ratios reduce to the same thing.
    # Second-order terms separate them by at most 1.6% (at the coarsest mesh,
    # |d phi| = 8e-3) and by 0.01% at the finest, so the curves overlap; hence
    # open squares over filled circles rather than two filled markers.
    ax.loglog(n, bias_per, "s", markerfacecolor="none", markeredgecolor=SERIES[1],
              markeredgewidth=1.8, ms=11, lw=0, zorder=4, label=r"$k_\perp$ (harmonic)")
    ax.loglog(n, bias_par, "o", color=SERIES[0], ms=7, lw=0, zorder=3,
              label=r"$k_\parallel$ (arithmetic)", markeredgecolor=t["surface"], markeredgewidth=1.0)
    # comp_eps.py's rule h = eps/sqrt(2) puts ~8.5 elements across the 5-95%
    # band (it reports this as "~7.5"). The same mesh is 13.0 across the 1-99%
    # band -- do not mix the conventions.
    ax.axvline(8.5, color=t["muted"], ls=":", lw=1.4, zorder=1)
    ax.annotate("project mesh rule\n$h=\\epsilon/\\sqrt{2}$  ($\\approx$8.5 elem)",
                xy=(8.0, bias_par.min() * 1.6), fontsize=8.5, color=t["ink2"], ha="right")
    ax.set_xlabel(r"elements across the $\phi=0.05\ldots0.95$ band")
    ax.set_ylabel(r"$|\Delta k| / k$ implied by the IC alone")
    ax.set_title("4.  Error the benchmark would misread as solver error\n"
                 "if compared against $\\phi=0.5$ instead of measured $\\bar\\phi$",
                 fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"], loc="upper right")

    fig.suptitle("Layered k$_{eff}$ benchmark — how much error belongs to the initial condition",
                 fontsize=13, color=t["ink"], x=0.012, ha="left", y=0.985)
    fig.text(0.012, 0.945,
             "1 mm periodic cell, 50/50 ice/air slab (-ic_type ice_slab). "
             "The continuous profile is exactly 50% ice for any $\\epsilon$; the discrete field is not.",
             fontsize=9, color=t["ink2"], ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=os.path.join(HERE, "ic_resolution_study.csv"))
    ap.add_argument("--theme", default="both", choices=["light", "dark", "both"])
    args = ap.parse_args()

    rows = load(args.csv)
    themes = ["light", "dark"] if args.theme == "both" else [args.theme]
    for th in themes:
        fig = make_figure(rows, th)
        out = os.path.join(HERE, f"ic_resolution_{th}.png")
        fig.savefig(out, dpi=160, facecolor=fig.get_facecolor())
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
