#!/usr/bin/env python3
"""Compare k_eff(t) between two eps on an identical initial condition.

Reads the k_eff.csv written by each arm of run_eps_convergence_dsm.sh and plots
the ratio R(t) = k_iso(t; safety 1.0) / k_iso(t; safety 0.5).

  R -> 1 as necks form   the eps dependence is a start-up transient introduced
                         by the initial condition, not by the physics
  R stays flat           the dependence is structural; an eps ladder is needed

Usage:
    python3 plot_eps_convergence_dsm.py --s05 <run_dir> --s10 <run_dir> [--out DIR]
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Validated categorical slots 1-3; see postprocess/validate_palette.py.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e",
                  muted="#8a8985", grid="#e3e2dd"),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7",
                  muted="#8a8985", grid="#333330"),
}


def load(run_dir):
    path = os.path.join(run_dir, "k_eff.csv")
    if not os.path.isfile(path):
        raise SystemExit(f"no k_eff.csv in {run_dir}")
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit(f"{path} is empty")
    return {
        "t":     np.array([float(r["time"]) for r in rows]),
        "step":  np.array([int(r["step"]) for r in rows]),
        "k_iso": np.array([float(r["k_iso"]) for r in rows]),
        "k00":   np.array([float(r["k_00"]) for r in rows]),
        "k11":   np.array([float(r["k_11"]) for r in rows]),
        "phi":   np.array([float(r["phi_bar"]) for r in rows]),
    }


def style(ax, t):
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


def make_figure(a, b, theme):
    """a = safety 0.5, b = safety 1.0"""
    t = THEMES[theme]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    fig.patch.set_facecolor(t["surface"])
    hr = 1.0 / 3600.0

    # --- 1. k_iso(t), both arms -------------------------------------------
    ax = axes[0]; style(ax, t)
    ax.plot(a["t"] * hr, a["k_iso"], "-", color=SERIES[0], lw=2,
            label=r"safety 0.5  ($\epsilon=1\,\mu$m)")
    ax.plot(b["t"] * hr, b["k_iso"], "-", color=SERIES[1], lw=2,
            label=r"safety 1.0  ($\epsilon=2\,\mu$m)")
    ax.set_xlabel("time [h]")
    ax.set_ylabel(r"$k_{iso}$  [W m$^{-1}$ K$^{-1}$]")
    ax.set_title("1.  Effective conductivity", fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"])

    # --- 2. the ratio, on a common time grid -------------------------------
    ax = axes[1]; style(ax, t)
    lo, hi = max(a["t"][0], b["t"][0]), min(a["t"][-1], b["t"][-1])
    tg = np.linspace(lo, hi, 400)
    R = np.interp(tg, b["t"], b["k_iso"]) / np.interp(tg, a["t"], a["k_iso"])
    ax.plot(tg * hr, R, "-", color=SERIES[2], lw=2.2, zorder=3)
    ax.axhline(1.0, color=t["muted"], ls="--", lw=1.4, zorder=2)
    ax.annotate("no eps dependence", xy=(tg[-1] * hr, 1.0),
                xytext=(-4, 6), textcoords="offset points",
                fontsize=8.5, color=t["ink2"], ha="right")
    ax.set_xlabel("time [h]")
    ax.set_ylabel(r"$k_{iso}(\epsilon=2\mu m)\,/\,k_{iso}(\epsilon=1\mu m)$")
    ax.set_title("2.  Does the eps dependence decay?\n"
                 "toward 1 = start-up transient; flat = structural",
                 fontsize=10.5, loc="left")
    if len(R):
        ax.text(0.03, 0.06,
                f"R(start) = {R[0]:.3f}\nR(end)   = {R[-1]:.3f}",
                transform=ax.transAxes, fontsize=9, color=t["ink2"],
                va="bottom", ha="left", family="monospace")

    # --- 3. anisotropy and ice fraction ------------------------------------
    ax = axes[2]; style(ax, t)
    ax.plot(a["t"] * hr, a["k00"] / a["k11"], "-", color=SERIES[0], lw=1.8,
            label=r"$k_{00}/k_{11}$, safety 0.5")
    ax.plot(b["t"] * hr, b["k00"] / b["k11"], "-", color=SERIES[1], lw=1.8,
            label=r"$k_{00}/k_{11}$, safety 1.0")
    ax.axhline(1.0, color=t["muted"], ls="--", lw=1.2, zorder=2)
    ax2 = ax.twiny()          # a second X axis is fine; never a second Y
    ax2.set_visible(False)
    ax.set_xlabel("time [h]")
    ax.set_ylabel(r"$k_{00}/k_{11}$   (1 = isotropic)")
    ax.set_title("3.  Anisotropy — a 28-grain cell is not isotropic,\n"
                 "and that is a REV question, not an eps one",
                 fontsize=10.5, loc="left")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=t["ink2"])

    fig.suptitle("Does the eps dependence of $k_{eff}$ survive DSM?  "
                 "0.5 mm packing, identical initial condition, 6 h at $-20\\,^\\circ$C",
                 fontsize=13, color=t["ink"], x=0.008, ha="left", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s05", required=True, help="run dir for safety 0.5")
    ap.add_argument("--s10", required=True, help="run dir for safety 1.0")
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--theme", default="both", choices=["light", "dark", "both"])
    args = ap.parse_args()

    a, b = load(args.s05), load(args.s10)

    lo, hi = max(a["t"][0], b["t"][0]), min(a["t"][-1], b["t"][-1])
    tg = np.linspace(lo, hi, 400)
    R = np.interp(tg, b["t"], b["k_iso"]) / np.interp(tg, a["t"], a["k_iso"])
    print(f"samples: safety 0.5 -> {len(a['t'])},  safety 1.0 -> {len(b['t'])}")
    print(f"overlap: {lo:.4g} .. {hi:.4g} s")
    print(f"R(start) = {R[0]:.4f}   R(end) = {R[-1]:.4f}   "
          f"min = {R.min():.4f}   max = {R.max():.4f}")
    drop = (R[0] - R[-1]) / max(abs(R[0] - 1.0), 1e-30)
    print(f"fraction of the initial excess removed by t_end: {100*drop:.1f}%")

    # Machine-readable, so the conclusion is not only in a picture.
    out_csv = os.path.join(args.out, "eps_convergence_dsm.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_s", "k_iso_s05", "k_iso_s10", "ratio"])
        for x in tg:
            k0 = float(np.interp(x, a["t"], a["k_iso"]))
            k1 = float(np.interp(x, b["t"], b["k_iso"]))
            w.writerow([f"{x:.6e}", f"{k0:.9e}", f"{k1:.9e}", f"{k1/k0:.9e}"])
    print("wrote", out_csv)

    for th in (["light", "dark"] if args.theme == "both" else [args.theme]):
        fig = make_figure(a, b, th)
        p = os.path.join(args.out, f"eps_convergence_dsm_{th}.png")
        fig.savefig(p, dpi=160, facecolor=fig.get_facecolor())
        plt.close(fig)
        print("wrote", p)


if __name__ == "__main__":
    main()
