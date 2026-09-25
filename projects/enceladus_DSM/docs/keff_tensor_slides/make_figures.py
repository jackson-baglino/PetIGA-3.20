#!/usr/bin/env python3
"""Build every figure for the k_eff tensor-law advisor slides (slides.md).

    venv_enceladus/bin/python docs/keff_tensor_slides/make_figures.py

Each panel is drawn at the size it is shown on the slide: half the width of a
16:9 slide, 6.4 x 4.6 in, with 18 pt axis labels. No equations are drawn in
the figures; the equations are in slides.md as LaTeX for IguanaTeX.

Data sources, none of them fitted here:
  analytic          disk:     studies/keff_sharp_limit/disk/keff_disk_analytic.py
                    laminate: studies/keff_sharp_limit/verification/keff_laminate_analytic.py
  disk eps ladder   studies/keff_sharp_limit/disk/keff_disk.csv
  disk f-sweep      studies/keff_sharp_limit/disk_sweep/keff_disk_sweep.csv   (HPC)
  laminate runs     ~/SimulationResults/enceladus_DSM/scratch/iceslab_2D_L1mm_eps20um_keff/
  pilot packing     pilot_rewiden.csv (pilot_rewiden.py)

Figures whose data is missing are skipped with a message, except the f-sweep
figures, which draw the first-order prediction alone and say so in the file
name (_prediction_only) until the sweep CSV exists.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
OUT = HERE / "figures"
sys.path.insert(0, str(PROJ / "studies/keff_sharp_limit/disk"))
sys.path.insert(0, str(PROJ / "studies/keff_sharp_limit/verification"))
import keff_disk_analytic as disk                              # noqa: E402
import keff_laminate_analytic as lam                            # noqa: E402

K_I, K_A = 2.29, 0.02
L_DISK = 1.0e-3
EPS_R_PROD = 0.02          # production: eps = 1 um on R_ave = 50 um

# ---- style ------------------------------------------------------------------
W, H = 6.4, 4.6            # one half-slide panel, inches
FS_LABEL, FS_TICK, FS_LEG, FS_NOTE = 18, 15, 13.5, 14
LW, LW_THIN, MS = 2.5, 1.6, 8
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#dddddd"
# Okabe-Ito (colour-blind safe); the repo standard, preprocess/figstyle.py.
C_ARI = "#D55E00"          # arithmetic law
C_TEN = "#0072B2"          # tensor law
C_SHP = "#009E73"          # thresholded (sharp) geometry
C_EXACT = INK              # exact / sharp-interface reference
TINT = {"arithmetic": {0.04: "#E8956A", 0.02: C_ARI, 0.01: "#F2C4A8"},
        "tensor": {0.04: "#5A9BD0", 0.02: C_TEN, 0.01: "#A9CBE8"}}
MARK = {0.04: "^", 0.02: "o", 0.01: "s"}

plt.rcParams.update({
    "font.size": FS_TICK, "axes.labelsize": FS_LABEL, "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK, "legend.fontsize": FS_LEG, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.linewidth": 1.0, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "savefig.facecolor": "white", "figure.facecolor": "white",
    "mathtext.fontset": "dejavusans",
})


def panel(w=W, h=H):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    return fig, ax


def plain_log_x(ax, ticks):
    from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def legend_below(fig, ax, ncol=2, extra=1.1, **kw):
    """Put the legend under the axes, growing the figure so the axes keep their size."""
    w, h = fig.get_size_inches()
    fig.set_size_inches(w, h + extra)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=ncol,
              fontsize=FS_LEG - 1, columnspacing=1.2, handlelength=2.2, **kw)


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.tight_layout()
    for ext, kw in (("png", {"dpi": 300}), ("pdf", {"metadata": {"CreationDate": None}})):
        fig.savefig(OUT / f"{name}.{ext}", **kw)
    plt.close(fig)
    print(f"  wrote figures/{name}.png/.pdf")


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def logistic(u):
    return 1.0 / (1.0 + np.exp(-u))


def k_arith(p):
    return K_A + (K_I - K_A) * p


def k_harm(p):
    return 1.0 / (p / K_I + (1.0 - p) / K_A)


# ---- slide 1 ----------------------------------------------------------------
def fig1a_slab():
    """Laminate: k across the layers vs eps, arithmetic law, exact at every eps."""
    L = 1.0e-3
    root = Path.home() / "SimulationResults/enceladus_DSM/scratch/iceslab_2D_L1mm_eps20um_keff"
    meas = []
    for d in sorted(root.glob("*_sharplimit_arith_L*")):
        f = d / "k_eff.csv"
        if f.is_file():
            r = read_csv(f)[-1]
            meas.append((L / int(d.name.rsplit("_L", 1)[1]), float(r["k_11"])))
    # keep the latest run per rung
    meas = sorted({e: k for e, k in meas}.items())

    fig, ax = panel()
    e = np.geomspace(L / 600, L / 45, 200)
    ax.plot(e * 1e6, [lam.k_perp(x, 0.5, L) for x in e], color=C_ARI, lw=LW,
            label="arithmetic law, closed form")
    if meas:
        ax.plot([m[0] * 1e6 for m in meas], [m[1] for m in meas], "o", color=C_ARI,
                ms=MS + 2, mec="white", mew=2, label="arithmetic law, solver")
    ks = lam.k_perp_sharp(0.5)
    ax.axhline(ks, color=C_EXACT, lw=LW_THIN + 0.4, ls="--", label="exact (sharp interface)")
    ax.set_xscale("log")
    plain_log_x(ax, [2, 5, 10, 20])
    ax.set_xlim(1.6, 24)
    ax.set_xlabel(r"interface width $\varepsilon$  [$\mu$m]")
    ax.set_ylabel(r"$k_\perp$  [W m$^{-1}$ K$^{-1}$]")
    ax.set_ylim(0.035, 0.068)
    for denom, dx, dy in ((50, -14, 4), (128, -6, 12), (512, 0, 12)):
        x = L / denom
        k = lam.k_perp(x, 0.5, L)
        ax.annotate(f"+{k / ks - 1:.0%}", (x * 1e6, k), xytext=(dx, dy),
                    textcoords="offset points", ha="right" if dx else "center",
                    fontsize=FS_NOTE, color=INK)
    ax.legend(loc="upper left")
    save(fig, "fig1a_slab_bias")
    if not meas:
        print("  (no laminate arithmetic runs found; drew the closed form only)")


def fig1b_pilot():
    path = HERE / "pilot_rewiden.csv"
    if not path.is_file():
        print("  skip fig1b: run pilot_rewiden.py first")
        return
    rows = read_csv(path)
    fig, ax = panel()
    shades = {}
    times = sorted({round(float(r["t_days"]), 2) for r in rows if float(r["t_days"]) > 0.5})
    for t, col in zip(times, ("#E8956A", C_ARI)):
        rs = [r for r in rows if abs(float(r["t_days"]) - t) < 0.01]
        e = np.array([float(r["eps"]) for r in rs]) * 1e6
        k = np.array([float(r["k_arithmetic"]) for r in rs])
        sl, ic = np.polyfit(e, k, 1)
        xx = np.linspace(0, e.max() * 1.05, 50)
        ax.plot(xx, ic + sl * xx, color=col, lw=LW_THIN, ls="--")
        ax.plot(e, k, "o", color=col, ms=MS + 2, mec="white", mew=2,
                label=f"day {t:.0f}: bias +{k[0] / ic - 1:.0%}")
        ax.plot([0], [ic], "D", color=col, ms=MS, mfc="white", mew=2)
        shades[t] = (k[0], ic)
    ax.axvline(1.0, color=MUTED, lw=1.2, ls=":")
    ax.text(1.05, ax.get_ylim()[0] + 0.02 * np.ptp(ax.get_ylim()), "production",
            color=MUTED, fontsize=FS_NOTE, va="bottom")
    ax.set_xlim(-0.1, 3.3)
    ax.set_xlabel(r"interface width $\varepsilon$  [$\mu$m]")
    ax.set_ylabel(r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]")
    ax.legend(loc="upper left", title="arithmetic law, pilot packing",
              title_fontsize=FS_LEG)
    save(fig, "fig1b_pilot_bias")
    (t0, (m0, c0)), (t1, (m1, c1)) = sorted(shades.items())
    print(f"  pilot: day {t0:.2f} measured {m0:.4f} extrapolated {c0:.4f} (+{m0/c0-1:.1%}); "
          f"day {t1:.2f} measured {m1:.4f} extrapolated {c1:.4f} (+{m1/c1-1:.1%})")
    print(f"  pilot rise: measured {m1/m0-1:+.1%}, extrapolated {c1/c0-1:+.1%}")


# ---- slide 2 ----------------------------------------------------------------
def _disk_field(eps, n=1200):
    x = (np.arange(n) + 0.5) / n * L_DISK
    X, Y = np.meshgrid(x, x)
    r = np.hypot(X - L_DISK / 2, Y - L_DISK / 2)
    R = 250e-6
    if eps == 0:
        return (r <= R).astype(float), x
    return 0.5 - 0.5 * np.tanh(0.5 / eps * (r - R)), x


def _disk_panel(eps, name, inset):
    R = 250e-6
    phi, x = _disk_field(eps)
    fig, ax = plt.subplots(figsize=(W, H))
    ext = [0, L_DISK * 1e3, 0, L_DISK * 1e3]
    im = ax.imshow(phi, origin="lower", extent=ext, cmap="Blues", vmin=0, vmax=1.15,
                   interpolation="bilinear")
    ax.set_xlabel("x  [mm]")
    ax.set_ylabel("y  [mm]")
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(True)
    ax.text(0.5, 0.5, "ice", color="white", fontsize=FS_LABEL, ha="center", va="center",
            transform=ax.transAxes, fontweight="bold")
    ax.text(0.1, 0.9, "air", color=INK, fontsize=FS_LABEL, ha="center", va="center",
            transform=ax.transAxes)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=[0, 0.5, 1])
    cb.set_label(r"ice indicator $\chi$" if eps == 0 else r"phase field $\varphi$",
                 fontsize=FS_LABEL)
    cb.ax.set_ylim(0, 1)
    cb.ax.tick_params(labelsize=FS_TICK)
    if inset:
        half = 30e-6
        cx, cy = (0.5 * L_DISK + R), 0.5 * L_DISK
        axin = ax.inset_axes([0.64, 0.05, 0.33, 0.33])
        n = 400
        xs = np.linspace(cx - half, cx + half, n)
        XX, YY = np.meshgrid(xs, np.linspace(cy - half, cy + half, n))
        rr = np.hypot(XX - L_DISK / 2, YY - L_DISK / 2)
        pz = ((rr <= R).astype(float) if eps == 0
              else 0.5 - 0.5 * np.tanh(0.5 / eps * (rr - R)))
        axin.imshow(pz, origin="lower", cmap="Blues", vmin=0, vmax=1.15,
                    extent=[(cx - half) * 1e3, (cx + half) * 1e3,
                            (cy - half) * 1e3, (cy + half) * 1e3])
        axin.set_xticks([]); axin.set_yticks([])
        for s in axin.spines.values():
            s.set_visible(True); s.set_color(INK); s.set_linewidth(1.5)
        ax.indicate_inset_zoom(axin, edgecolor=INK, lw=1.5, alpha=1)
        axin.text(0.5, -0.04, r"60 $\mu$m zoom", transform=axin.transAxes, ha="center",
                  va="top", fontsize=FS_NOTE - 1, color=INK)
    fig.tight_layout()
    save(fig, name)


def fig2_disks():
    _disk_panel(0.0, "fig2a_disk_sharp", inset=True)
    _disk_panel(EPS_R_PROD * 250e-6, "fig2b_disk_diffuse", inset=True)


def fig2c_profile():
    """1-D profile across the interface: phi vs sharp step, in units of eps."""
    u = np.linspace(-8, 8, 600)
    fig, ax = panel(W, H * 0.8)
    ax.plot(u, np.where(u < 0, 0, 1), color=C_EXACT, lw=LW, ls="--", label=r"sharp  $\chi$")
    ax.plot(u, logistic(u), color=C_TEN, lw=LW, label=r"diffuse  $\varphi$")
    ax.axvspan(-4.6, 4.6, color=C_TEN, alpha=0.08, lw=0)
    ax.text(0, 0.5, r"visible band $\approx 9\,\varepsilon$", ha="center", va="center",
            fontsize=FS_NOTE, color=INK, bbox=dict(fc="white", ec="none", pad=2))
    ax.set_xlabel(r"distance across interface  $s/\varepsilon$")
    ax.set_ylabel("ice fraction")
    ax.set_yticks([0, 0.5, 1])
    ax.set_xlim(-8, 8)
    ax.text(-7.6, 0.08, "air", fontsize=FS_NOTE, color=MUTED)
    ax.text(7.6, 0.85, "ice", fontsize=FS_NOTE, color=MUTED, ha="right")
    ax.legend(loc="upper left")
    save(fig, "fig2c_profile")


# ---- slides 3 and 6: the f-sweep --------------------------------------------
def _sweep():
    p = PROJ / "studies/keff_sharp_limit/disk_sweep/keff_disk_sweep.csv"
    return read_csv(p) if p.is_file() else []


def _pred_err(f, ratio):
    """First-order relative bias of the arithmetic law, dipole field."""
    R = L_DISK * math.sqrt(f / math.pi)
    return disk.arith_slope(R, L_DISK) * ratio * R / disk.k_sharp(f)


def _ratio_label(r):
    return (r"$\varepsilon/R$ = " + f"{r:g}" + ("  (production)" if r == EPS_R_PROD else ""))


def fig_sweep(law, fade_arith=False):
    rows = [r for r in _sweep() if r["law"] == law]
    tag = "fig3" if law == "arithmetic" else "fig6"
    suffix = "" if rows else "_prediction_only"
    if not rows:
        if law == "tensor":
            print(f"  skip {tag}: no sweep CSV yet (there is no tensor prediction to draw)")
            return
        print(f"  {tag}: no sweep CSV yet, drawing the first-order prediction only")
    ff = np.linspace(0.02, 0.55, 200)

    # (a) k_eff vs f
    fig, ax = panel()
    ax.plot(ff, [disk.k_sharp(f) for f in ff], color=C_EXACT, lw=LW,
            label="exact (Rayleigh)", zorder=5)
    if fade_arith:
        ar = [r for r in _sweep() if r["law"] == "arithmetic"
              and float(r["eps_over_R"]) == EPS_R_PROD]
        if ar:
            ax.plot([float(r["f"]) for r in ar], [float(r["k"]) for r in ar], "-o",
                    color=C_ARI, alpha=0.35, lw=LW, ms=MS - 2,
                    label=r"arithmetic law, $\varepsilon/R$ = 0.02")
    for ratio in (0.04, 0.02, 0.01):
        rs = sorted((r for r in rows if float(r["eps_over_R"]) == ratio),
                    key=lambda r: float(r["f"]))
        prod = ratio == EPS_R_PROD
        if rs:
            ax.plot([float(r["f"]) for r in rs], [float(r["k"]) for r in rs],
                    "-" + MARK[ratio], color=TINT[law][ratio], lw=LW + (1 if prod else -0.5),
                    ms=MS + (2 if prod else 0), mec="white", mew=1.5,
                    zorder=6 if prod else 4, label=f"{law}, " + _ratio_label(ratio))
        elif law == "arithmetic":
            ax.plot(ff, [disk.k_sharp(f) * (1 + _pred_err(f, ratio)) for f in ff], ":",
                    color=TINT[law][ratio], lw=LW + (1 if prod else -0.5),
                    label="first-order prediction, " + _ratio_label(ratio))
    ax.set_xlabel(r"ice area fraction $f$")
    ax.set_ylabel(r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]")
    ax.set_xlim(0, 0.56)
    ax.legend(loc="upper left", fontsize=FS_LEG - 1.5)
    save(fig, f"{tag}a_keff_vs_f_{law}{suffix}")

    # (b) relative error vs f
    fig, ax = panel()
    ax.axhline(0, color=C_EXACT, lw=LW_THIN)
    for ratio in (0.04, 0.02, 0.01):
        prod = ratio == EPS_R_PROD
        col = TINT[law][ratio]
        if law == "arithmetic":
            ax.plot(ff, [100 * _pred_err(f, ratio) for f in ff], "--", color=col,
                    lw=LW_THIN + (0.8 if prod else 0))
        rs = sorted((r for r in rows if float(r["eps_over_R"]) == ratio),
                    key=lambda r: float(r["f"]))
        if rs:
            ax.plot([float(r["f"]) for r in rs], [100 * float(r["rel_err"]) for r in rs],
                    "-" + MARK[ratio], color=col, lw=LW + (1 if prod else -0.5),
                    ms=MS + (2 if prod else 0), mec="white", mew=1.5,
                    label=_ratio_label(ratio))
        elif law == "arithmetic":
            ax.plot([], [], "--", color=col, lw=LW, label=_ratio_label(ratio))
    if law == "arithmetic":
        ax.plot([], [], "--", color=MUTED, lw=LW_THIN, label="first-order prediction")
    ax.set_xlabel(r"ice area fraction $f$")
    ax.set_ylabel(r"error vs exact  [%]")
    ax.set_xlim(0, 0.56)
    if law == "tensor":
        ax.set_ylim(-5, 50)   # same scale as slide 3 would hide nothing; keep it comparable
    ax.legend(loc="upper left", title=f"{law} law", title_fontsize=FS_LEG)
    save(fig, f"{tag}{'b' if law == 'arithmetic' else 'c'}_error_vs_f_{law}{suffix}")


def fig6b_ladder():
    rows = read_csv(PROJ / "studies/keff_sharp_limit/disk/keff_disk.csv")
    R = 250e-6
    ks = disk.k_sharp(disk.area_fraction(R, L_DISK))
    s_pred = disk.arith_slope(R, L_DISK)
    fig, ax = panel()
    for law, col, mk, lab in (("arith", C_ARI, "o", "arithmetic law"),
                              ("tensor", C_TEN, "s", "tensor law")):
        rs = sorted((r for r in rows if r["interp"] == law), key=lambda r: float(r["eps"]))
        e = np.array([float(r["eps"]) for r in rs]) / R
        err = np.array([abs(0.5 * (float(r["k_00"]) + float(r["k_11"])) / ks - 1) for r in rs])
        ax.plot(e, 100 * err, mk + "-", color=col, lw=LW, ms=MS + 2, mec="white", mew=1.5,
                label=lab, zorder=5)
        if law == "tensor":
            ee = np.geomspace(e.min() * 0.8, e.max() * 1.25, 50)
            ax.plot(ee, 100 * err[0] * (ee / e[0]) ** 2, ":", color=col, lw=LW_THIN,
                    label=r"slope 2 guide ($\propto\varepsilon^2$)")
    ee = np.geomspace(0.004, 0.05, 50)
    ax.plot(ee, 100 * s_pred * ee * R / ks, "--", color=C_ARI, lw=LW_THIN,
            label=r"first-order prediction ($\propto\varepsilon$)")
    ax.axvline(EPS_R_PROD, color=MUTED, lw=1.2, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    plain_log_x(ax, [0.005, 0.01, 0.02, 0.04])
    ax.text(EPS_R_PROD * 1.06, 0.0045, "production", color=MUTED, fontsize=FS_NOTE)
    ax.set_ylim(0.003, 40)
    ax.set_xlabel(r"interface width $\varepsilon/R$")
    ax.set_ylabel(r"|error| vs exact  [%]")
    legend_below(fig, ax, ncol=2, extra=1.3)
    save(fig, "fig6b_error_vs_eps_ladder")


# ---- slide 4: the excesses of the arithmetic law ------------------------------
U = np.linspace(-7, 7, 1401)
H_U = (U >= 0).astype(float)
P_U = logistic(U)


def _step(ax, y_air, y_ice, label):
    ax.plot([U[0], 0, 0, U[-1]], [y_air, y_air, y_ice, y_ice], color=C_EXACT, lw=LW,
            ls="--", label=label, zorder=5)


def _iface_axes(ax, ylab):
    ax.set_xlabel(r"distance across interface  $s/\varepsilon$")
    ax.set_ylabel(ylab)
    ax.set_xlim(U[0], U[-1])
    ax.text(0.03, 0.97, "air", transform=ax.transAxes, fontsize=FS_NOTE, color=MUTED,
            va="top")
    ax.text(0.97, 0.97, "ice", transform=ax.transAxes, fontsize=FS_NOTE, color=MUTED,
            va="top", ha="right")


def fig4a_resistivity_arith():
    fig, ax = panel()
    r_sharp = 1.0 / np.where(H_U > 0, K_I, K_A)
    r_ar = 1.0 / k_arith(P_U)
    ax.fill_between(U, r_ar, r_sharp, color=C_ARI, alpha=0.18, lw=0,
                    label="missing resistance")
    _step(ax, 1 / K_A, 1 / K_I, "sharp interface")
    ax.plot(U, r_ar, color=C_ARI, lw=LW, label="arithmetic law")
    _iface_axes(ax, r"resistivity $1/K$  [m K W$^{-1}$]")
    ax.set_ylim(-2, 58)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62))
    save(fig, "fig4a_resistivity_across_arithmetic")


def fig4b_conductivity_arith():
    fig, ax = panel()
    k_sharp = np.where(H_U > 0, K_I, K_A)
    k_ar = k_arith(P_U)
    ax.fill_between(U, k_ar, k_sharp, where=U < 0, color=C_ARI, alpha=0.18, lw=0)
    ax.fill_between(U, k_ar, k_sharp, where=U >= 0, color=C_ARI, alpha=0.18, lw=0)
    _step(ax, K_A, K_I, "sharp interface")
    ax.plot(U, k_ar, color=C_ARI, lw=LW, label="arithmetic law")
    ax.text(-2.0, 0.35, "+", fontsize=26, ha="center", va="center", color=INK)
    ax.text(2.0, 1.95, "−", fontsize=26, ha="center", va="center", color=INK)
    ax.text(0.73, 0.40, "equal areas:\nthey cancel", transform=ax.transAxes, ha="center",
            fontsize=FS_NOTE, color=INK)
    _iface_axes(ax, r"conductivity $K$  [W m$^{-1}$ K$^{-1}$]")
    ax.set_ylim(-0.25, 2.6)
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62))
    save(fig, "fig4b_conductivity_along_arithmetic")


# ---- slide 5: the tensor law --------------------------------------------------
def fig5a_tensor_profile():
    fig, ax = panel()
    _step(ax, K_A, K_I, "sharp interface")
    ax.plot(U, k_arith(P_U), color=C_ARI, lw=LW,
            label="tensor law, along (arithmetic)\n= old law in every direction")
    ax.plot(U, k_harm(P_U), color=C_TEN, lw=LW + 0.5,
            label="tensor law, across (harmonic)")
    ax.set_yscale("log")
    _iface_axes(ax, r"conductivity  [W m$^{-1}$ K$^{-1}$]")
    ax.set_ylim(0.012, 6)
    legend_below(fig, ax, ncol=1, extra=1.25)
    save(fig, "fig5a_tensor_K_across_interface")


def fig5b_resistivity_tensor():
    fig, ax = panel()
    r_sharp = 1.0 / np.where(H_U > 0, K_I, K_A)
    r_h = 1.0 / k_harm(P_U)
    ax.fill_between(U, r_h, r_sharp, color=C_TEN, alpha=0.18, lw=0)
    _step(ax, 1 / K_A, 1 / K_I, "sharp interface")
    ax.plot(U, 1.0 / k_arith(P_U), color=C_ARI, lw=LW_THIN, alpha=0.6,
            label="arithmetic law")
    ax.plot(U, r_h, color=C_TEN, lw=LW, label="tensor law (across)")
    ax.text(-2.2, 38, "−", fontsize=26, ha="center", va="center", color=INK)
    ax.text(2.2, 12, "+", fontsize=26, ha="center", va="center", color=INK)
    ax.text(0.74, 0.40, "equal areas:\nthey cancel", transform=ax.transAxes, ha="center",
            fontsize=FS_NOTE, color=INK)
    _iface_axes(ax, r"resistivity across  [m K W$^{-1}$]")
    ax.set_ylim(-2, 58)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.93))
    save(fig, "fig5b_resistivity_across_tensor")


def fig5c_layers():
    """Schematic: the band as a stack of thin layers parallel to the interface."""
    fig, ax = plt.subplots(figsize=(W, H))
    n = 14
    cmap = plt.get_cmap("Blues")
    for i in range(n):
        p = logistic(-6 + 12 * (i + 0.5) / n)
        ax.add_patch(Rectangle((0.08, 0.1 + 0.8 * i / n), 0.56, 0.8 / n,
                               fc=cmap(0.05 + 0.85 * p), ec="white", lw=1.0))
    ax.text(0.36, 0.05, "air side", ha="center", va="center", fontsize=FS_NOTE, color=INK)
    ax.text(0.36, 0.95, "ice side", ha="center", va="center", fontsize=FS_NOTE, color=INK)
    ax.add_patch(FancyArrowPatch((0.12, 0.5), (0.60, 0.5), arrowstyle="-|>",
                                 mutation_scale=28, lw=3, color=C_ARI))
    ax.text(0.36, 0.555, "along: layers in parallel", ha="center", fontsize=FS_NOTE,
            color=INK, bbox=dict(fc="white", ec="none", pad=2))
    ax.add_patch(FancyArrowPatch((0.72, 0.12), (0.72, 0.88), arrowstyle="-|>",
                                 mutation_scale=28, lw=3, color=C_TEN))
    ax.text(0.76, 0.62, "across:\nlayers in\nseries", ha="left", va="center",
            fontsize=FS_NOTE, color=INK)
    ax.text(0.76, 0.30, "→ harmonic\n    mean", ha="left", va="center",
            fontsize=FS_NOTE, color=C_TEN, fontweight="bold")
    ax.text(0.36, 0.44, "→ arithmetic mean", ha="center", fontsize=FS_NOTE,
            color=C_ARI, fontweight="bold", bbox=dict(fc="white", ec="none", pad=2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    save(fig, "fig5c_layer_schematic")


# ---- slide 7: the pilot packing, three laws -----------------------------------
def fig7_pilot():
    path = HERE / "pilot_rewiden.csv"
    if not path.is_file():
        print("  skip fig7: run pilot_rewiden.py first")
        return
    rows = read_csv(path)
    times = sorted({round(float(r["t_days"]), 2) for r in rows if float(r["t_days"]) > 0.5})
    allk = [float(r[k]) for r in rows for k in ("k_arithmetic", "k_tensor", "k_sharp")
            if float(r["t_days"]) > 0.5]
    ylim = (min(allk) * 0.95, max(allk) * 1.03)
    for t, tag in zip(times, "ab"):
        rs = [r for r in rows if abs(float(r["t_days"]) - t) < 0.01]
        e = np.array([float(r["eps"]) for r in rs]) * 1e6
        fig, ax = panel()
        for key, col, mk, lab in (("k_arithmetic", C_ARI, "o", "arithmetic law"),
                                  ("k_tensor", C_TEN, "s", "tensor law"),
                                  ("k_sharp", C_SHP, "D", "thresholded at 1/2 (sharp)")):
            k = np.array([float(r[key]) for r in rs])
            ax.plot(e, k, mk + "-", color=col, lw=LW, ms=MS + 2, mec="white", mew=1.5,
                    label=lab)
            print(f"  pilot day {t:.2f} {key}: {k[0]:.4f} -> {k[-1]:.4f} "
                  f"({k[-1]/k[0]-1:+.1%} over 3x eps)")
        ax.axvline(1.0, color=MUTED, lw=1.2, ls=":")
        ax.set_xlim(0.8, 3.2)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r"interface width $\varepsilon$  [$\mu$m]")
        ax.set_ylabel(r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]")
        ax.legend(loc="upper left", title=f"pilot packing, day {t:.0f}",
                  title_fontsize=FS_LEG, fontsize=FS_LEG - 1)
        save(fig, f"fig7{'cd'['ab'.index(tag)]}_pilot_rewiden_day{t:.0f}")


# ---- slide 7: the PetIGA replay of the pilot packings -------------------------
HPC = Path.home() / "SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering"
REPLAY_DIRS = (HPC / "batch_2026-09-16__13.16.11_pilot_keff",     # arithmetic, in-line
               HPC / "batch_2026-09-24__12.18.18_keff_replay")   # tensor, replayed
# Seed 2's leg-1 tensor replay was lost to a node failure: no 1-day baseline.
ALL_SEEDS = ("1", "2", "3", "4")
MOVIE_SEED = "3"                 # the seed shown as a movie in the talk
SEED_LS = {"1": "-", "2": "-", "3": "--", "4": ":"}
DAY = 86400.0


def _replay():
    sys.path.insert(0, str(PROJ / "studies/keff_sintering/coefficient_fix"))
    import compare_laws as cl
    runs: dict = {}
    for b in REPLAY_DIRS:
        if not b.is_dir():
            return None, None
        for seed, laws in cl.collect(b).items():
            for law, parts in laws.items():
                runs.setdefault(seed, {}).setdefault(law, []).extend(parts)
    merged = {sd: {law: cl.merge_legs(parts) for law, parts in laws.items()}
              for sd, laws in runs.items() if sd in ALL_SEEDS}
    return merged, cl


def _ensemble(runs, cl):
    """Seeds whose BOTH laws start by the baseline, as compare_laws.py requires.

    Seed 2's leg-1 tensor replay was lost to a node failure, so until it is
    rerun its tensor curve starts at day 16 and it is left out; once a full
    k_eff_tensor.csv exists it joins automatically."""
    return tuple(sd for sd in ALL_SEEDS if sd in runs and
                 all(law in runs[sd] and
                     float(runs[sd][law]["time"][0]) <= cl.BASELINE_TOL * DAY
                     for law in ("arith", "tensor")))


def fig7_replay():
    runs, cl = _replay()
    if not runs:
        print("  skip fig7a/b: replay batch directories not found")
        return
    from matplotlib.lines import Line2D
    SEEDS = tuple(sd for sd in _ensemble(runs, cl) if sd != MOVIE_SEED)
    movie_in = MOVIE_SEED in _ensemble(runs, cl)
    ens_label = ", ".join(sorted(SEEDS + ((MOVIE_SEED,) if movie_in else ())))
    styles = (("arith", C_ARI, "arithmetic law"), ("tensor", C_TEN, "tensor law"))

    def seed_handles(extra_movie_label):
        h = [Line2D([], [], color=C_ARI, lw=LW, label="arithmetic law"),
             Line2D([], [], color=C_TEN, lw=LW, label="tensor law"),
             Line2D([], [], color=INK, lw=LW + 1.5, label=extra_movie_label)]
        h += [Line2D([], [], color=MUTED, lw=LW - 1, ls=SEED_LS[sd], label=f"seed {sd}")
              for sd in SEEDS]
        return h

    # (a) absolute k_eff(t): seed 2 (the movie) bold, the others thin and styled
    fig, ax = panel()
    for law, col, _ in styles:
        for sd in SEEDS:
            a = runs[sd][law]
            ax.plot(np.asarray(a["time"]) / DAY, a["k_iso"], color=col, lw=LW - 1.2,
                    ls=SEED_LS[sd], alpha=0.55)
        a = runs[MOVIE_SEED][law]
        ax.plot(np.asarray(a["time"]) / DAY, a["k_iso"], color=col, lw=LW + 1.5,
                solid_capstyle="round", zorder=6)
    t2 = np.asarray(runs[MOVIE_SEED]["tensor"]["time"]) / DAY
    if not movie_in:
        ax.annotate(f"seed {MOVIE_SEED} tensor replay\nlost to a node failure;\nrerun pending",
                    xy=(t2[0], runs[MOVIE_SEED]["tensor"]["k_iso"][0]), xytext=(4.5, 0.47),
                    fontsize=FS_NOTE - 2, color=MUTED, va="center",
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.2))
    ax.axvspan(0, 1, color=MUTED, alpha=0.10, lw=0)
    ax.set_xlim(0, 30)
    ax.set_xlabel("time  [days]")
    ax.set_ylabel(r"$k_\mathrm{eff}$  [W m$^{-1}$ K$^{-1}$]")
    legend_below(fig, ax, ncol=3, extra=1.2, handles=seed_handles(f"seed {MOVIE_SEED} (movie)"))
    save(fig, "fig7a_pilot_keff_vs_time")

    # (b) rise from the 1-day baseline: the ensemble, with seed 2 bold
    fig, ax = panel()
    for law, col, lab in styles:
        rises = []
        if movie_in:
            a = runs[MOVIE_SEED][law]
            kb, ke, r, tb = cl.rise(a, 1.0)
            t = np.asarray(a["time"]) / DAY
            m = t >= tb / DAY
            ax.plot(t[m], 100 * (np.asarray(a["k_iso"])[m] / kb - 1), color=col,
                    lw=LW + 1.5, zorder=6)
            rises.append(r)
            print(f"  seed {MOVIE_SEED} (movie) {law}: rise {r:+.1f}%")
        for sd in SEEDS:
            a = runs[sd][law]
            kb, ke, r, tb = cl.rise(a, 1.0)
            t = np.asarray(a["time"]) / DAY
            m = t >= tb / DAY
            ax.plot(t[m], 100 * (np.asarray(a["k_iso"])[m] / kb - 1), color=col,
                    lw=LW - 1.2, ls=SEED_LS[sd], alpha=0.7)
            rises.append(r)
        mean, sd_ = np.mean(rises), np.std(rises, ddof=1)
        ax.plot([], [], color=col, lw=LW, label=f"{lab}: +{mean:.1f}%  (sd {sd_:.1f})")
        print(f"  replay {law}: rises {', '.join(f'{x:+.1f}%' for x in rises)}; "
              f"mean {mean:+.1f}% sd {sd_:.1f}")
    if movie_in:
        ax.plot([], [], color=INK, lw=LW + 1.5, label=f"seed {MOVIE_SEED} (movie)")
    else:
        a = runs[MOVIE_SEED]["arith"]
        kb, ke, r2, tb = cl.rise(a, 1.0)
        t = np.asarray(a["time"]) / DAY
        m = t >= tb / DAY
        ax.plot(t[m], 100 * (np.asarray(a["k_iso"])[m] / kb - 1), color=C_ARI, lw=LW + 1.5,
                zorder=6, label=f"seed {MOVIE_SEED} (movie), arithmetic: +{r2:.1f}%")
    ax.axhline(0, color=C_EXACT, lw=LW_THIN)
    ax.set_xlim(0, 30)
    ax.set_xlabel("time  [days]")
    ax.set_ylabel(r"rise of $k_\mathrm{eff}$ since day 1  [%]")
    h, _ = ax.get_legend_handles_labels()
    h += [Line2D([], [], color=MUTED, lw=LW - 1, ls=SEED_LS[sd], label=f"seed {sd}")
          for sd in SEEDS]
    legend_below(fig, ax, ncol=2, extra=1.5, handles=h,
                 title=f"rise by day 30 (mean and sd of seeds {ens_label})",
                 title_fontsize=FS_LEG - 1)
    save(fig, "fig7b_pilot_rise")

    for sd in sorted(SEEDS + (MOVIE_SEED,)):
        a, t_ = runs[sd]["arith"], runs[sd]["tensor"]
        print(f"  seed {sd}: k(0) arithmetic {a['k_iso'][0]:.4f}, tensor {t_['k_iso'][0]:.4f}"
              f"  (t0 = {a['time'][0]/DAY:.2f} / {t_['time'][0]/DAY:.2f} d)")
    if movie_in:
        return
    # A partial movie seed's only like-for-like comparison: both laws over the tensor window.
    at, tt = runs[MOVIE_SEED]["arith"], runs[MOVIE_SEED]["tensor"]
    t0, t1 = float(tt["time"][0]), float(tt["time"][-1])
    ta = np.asarray(at["time"])
    ka0 = np.interp(t0, ta, at["k_iso"]); ka1 = np.interp(t1, ta, at["k_iso"])
    print(f"  seed {MOVIE_SEED}, day {t0/DAY:.1f} -> {t1/DAY:.1f}: arithmetic {ka0:.4f} -> {ka1:.4f} "
          f"({ka1/ka0-1:+.1%}), tensor {tt['k_iso'][0]:.4f} -> {tt['k_iso'][-1]:.4f} "
          f"({tt['k_iso'][-1]/tt['k_iso'][0]-1:+.1%})")


def main():
    print("slide 1"); fig1a_slab(); fig1b_pilot()
    print("slide 2"); fig2_disks(); fig2c_profile()
    print("slide 3"); fig_sweep("arithmetic")
    print("slide 4"); fig4a_resistivity_arith(); fig4b_conductivity_arith()
    print("slide 5"); fig5a_tensor_profile(); fig5b_resistivity_tensor(); fig5c_layers()
    print("slide 6"); fig_sweep("tensor", fade_arith=True); fig6b_ladder()
    print("slide 7"); fig7_replay(); fig7_pilot()


if __name__ == "__main__":
    main()
