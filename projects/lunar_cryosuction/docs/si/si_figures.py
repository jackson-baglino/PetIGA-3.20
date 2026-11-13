"""Schematics for the supplementary-methods documents (docs/si/).

    venv_cryosuction/bin/python docs/si/si_figures.py <out_dir>

fig_throat: the pore with a throat, a lens spanning it, its warm and cold faces,
the gas columns to the two ends and the end vapor condition; checked with
postprocess/figcheck.py at the placed width.
"""
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "postprocess"))
from figcheck import check_figure  # noqa: E402

INK, GRAY, ROCK, ICE, ICE_E = "#222222", "#555555", "#d9cbb3", "#9ecae1", "#1f4e8c"
W_IN = 6.3


def fig_throat(out: Path) -> int:
    plt.rcParams.update({"font.family": "Arial", "font.size": 8})
    L, W, wt = 400.0, 150.0, 40.0
    c = math.sqrt(W * W - wt * wt) / (L / 2)
    x = np.linspace(0, L, 801)
    xi = x - L / 2
    g = np.sqrt(wt * wt + (c * xi) ** 2)
    yb, yt = (W - g) / 2, (W + g) / 2
    fig, ax = plt.subplots(figsize=(W_IN, 3.4))
    ax.fill_between(x, -12, yb, color=ROCK, lw=0)
    ax.fill_between(x, yt, W + 12, color=ROCK, lw=0)
    ax.plot(x, yb, color=INK, lw=0.9)
    ax.plot(x, yt, color=INK, lw=0.9)
    # lens: faces as arcs meeting the walls; warm face at xw, cold at xc (drawn flat-ish arcs)
    xw, xc = 165.0, 255.0
    def gap(xx):
        return math.sqrt(wt * wt + (c * (xx - L / 2)) ** 2)
    ys = np.linspace(0, 1, 60)
    def face(x0, sign):
        gg = gap(x0)
        y = (W - gg) / 2 + ys * gg
        bulge = sign * 0.18 * gg * (1 - (2 * ys - 1) ** 2)
        return x0 + bulge, y
    fwx, fwy = face(xw, +1)
    fcx, fcy = face(xc, -1)
    mask = (x >= xw - 10) & (x <= xc + 10)
    px = np.r_[fwx, x[(x > xw) & (x < xc)], fcx[::-1], x[(x > xw) & (x < xc)][::-1]]
    py = np.r_[fwy, yt[(x > xw) & (x < xc)], fcy[::-1], yb[(x > xw) & (x < xc)][::-1]]
    ax.fill(px, py, color=ICE, ec="none", zorder=2)
    ax.plot(fwx, fwy, color=ICE_E, lw=1.2, zorder=3)
    ax.plot(fcx, fcy, color=ICE_E, lw=1.2, zorder=3)
    del mask
    # ends
    for xe in (0, L):
        ax.plot([xe, xe], [0 - 2, W + 2], color=INK, lw=1.0, ls=(0, (3, 2)))
    ax.text(4, W + 20, "end at $x = 0$ (warm):\n$\\rho_v = h\\,\\rho_{vs}(T(0))$", ha="left", va="bottom",
            color=INK)
    ax.text(L - 4, W + 20, "end at $x = L$ (cold):\n$\\rho_v = h\\,\\rho_{vs}(T(L))$", ha="right",
            va="bottom", color=INK)
    # labels (face labels below the pore, with leaders; nothing written across a wall line)
    ax.text((xw + xc) / 2, W / 2, "ice lens\n(area $A$)", ha="center", va="center", color=ICE_E,
            gid="overlap_ok")
    ax.annotate("throat: narrowest gap $w_t$", (L / 2, (W + wt) / 2 + 1), (L / 2 + 70, W + 6), color=INK,
                ha="left", va="center", fontsize=7.5,
                arrowprops=dict(arrowstyle="->", color=INK, lw=0.6))
    ax.annotate("warm face $w$", (xw - 2, (W - gap(xw)) / 2 + 6), (xw - 40, -28), color=ICE_E,
                arrowprops=dict(arrowstyle="->", color=ICE_E, lw=0.7), ha="center", va="top")
    ax.annotate("cold face $c$", (xc + 2, (W - gap(xc)) / 2 + 6), (xc + 40, -28), color=ICE_E,
                arrowprops=dict(arrowstyle="->", color=ICE_E, lw=0.7), ha="center", va="top")
    ya = W / 2
    ax.annotate("", (3, ya), (xw - 22, ya), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.6))
    ax.text((xw - 19) / 2, ya + 4, "gas column of $w$", ha="center", va="bottom", color=GRAY)
    ax.annotate("", (xc + 22, ya), (L - 3, ya), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.6))
    ax.text((xc + 19 + L) / 2, ya + 4, "gas column of $c$", ha="center", va="bottom", color=GRAY)
    ax.annotate("", (14, yb[28] + 1), (14, yt[28] - 1), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.6))
    ax.text(18, W / 2 - 14, "gap $g(x)$", ha="left", va="top", color=INK)
    ax.text(60, W + 6, "regolith", ha="center", va="center", color=GRAY, fontsize=7.5)
    ax.annotate("", (L, -48), (0, -48), arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.6))
    ax.text(L / 2, -52, "$x$; with a gradient, $T$ decreases toward $x = L$", ha="center", va="top", color=GRAY)
    ax.set_xlim(-10, L + 10)
    ax.set_ylim(-68, W + 55)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(0.01, 0.01, 0.99, 0.99)
    bad = check_figure(fig, min_pt=7.0, placed_width_in=W_IN)
    for b in bad:
        print("  FAIL:", b)
    fig.savefig(out, dpi=300)
    print(f"{out}: {len(bad)} problems")
    return len(bad)



def fig_thermal(out: Path) -> int:
    """Pore with heated walls: the walls are faces of the mesh, the regolith is outside it; s runs along
    the wall from the heated end; inset: T_w(s) at three times (erfc profile)."""
    from math import erfc
    plt.rcParams.update({"font.family": "Arial", "font.size": 8})
    fig = plt.figure(figsize=(W_IN, 2.6))
    ax = fig.add_axes((0.01, 0.06, 0.56, 0.9))
    Lx, W = 300.0, 100.0
    ax.fill_between([0, Lx], W, W + 22, color=ROCK, lw=0)
    ax.fill_between([0, Lx], -22, 0, color=ROCK, lw=0)
    ax.plot([0, Lx], [W, W], color="#b2182b", lw=2.0)
    ax.plot([0, Lx], [0, 0], color="#b2182b", lw=2.0)
    ax.plot([0, 0], [0, W], color=INK, lw=1.0, ls=(0, (3, 2)))
    ax.plot([Lx, Lx], [0, W], color=INK, lw=1.0, ls=(0, (3, 2)))
    ax.text(Lx / 2 + 20, W + 11, "regolith (outside the mesh)", ha="center", va="center", color=GRAY)
    ax.text(Lx / 2, W / 2 + 12, "pore: ice and gas (the mesh)", ha="center", va="center", color=INK)
    ax.text(Lx / 2, -11, "wall face $\Gamma_w$: $T = T_w(s,t)$", ha="center", va="center", color="#b2182b")
    ax.annotate("", (8, W / 2 - 18), (Lx - 8, W / 2 - 18), arrowprops=dict(arrowstyle="<-", color=GRAY, lw=0.6))
    ax.text(Lx / 2, W / 2 - 23, "$s = x + z_0$: depth below the heated surface", ha="center", va="top",
            color=GRAY)
    ax.annotate("", (40, W + 19), (40, W + 1), arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))
    ax.text(45, W + 15, "$n$", ha="left", va="center", color=INK)
    ax.annotate("", (40, -19), (40, -1), arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))
    ax.text(45, -15, "$n$", ha="left", va="center", color=INK)
    ax.text(4, W / 2 + 32, "toward the\nsurface", ha="left", va="center", color="#b2182b")
    ax.set_xlim(-6, Lx + 6)
    ax.set_ylim(-26, W + 26)
    ax.axis("off")
    a2 = fig.add_axes((0.67, 0.22, 0.31, 0.68))
    s = np.linspace(0, 1, 200)
    for tau, ls in ((0.01, "-"), (0.05, "--"), (0.2, ":")):
        a2.plot(s, [erfc(v / (2 * math.sqrt(tau))) for v in s], color=INK, ls=ls, lw=1.0,
                label=f"$\\kappa_r (t - t_0) = {tau:g}\\,\\ell_s^2$")
    a2.set_xlabel("$s/\\ell_s$")
    a2.set_ylabel("$(T_w - T_{init})/\\Delta T$")
    a2.legend(fontsize=7, frameon=False, loc="upper right")
    a2.tick_params(labelsize=7, length=2)
    bad = check_figure(fig, min_pt=7.0, placed_width_in=W_IN)
    for b in bad:
        print("  FAIL:", b)
    fig.savefig(out, dpi=300)
    print(f"{out}: {len(bad)} problems")
    return len(bad)


if __name__ == "__main__":
    d = Path(sys.argv[1])
    d.mkdir(parents=True, exist_ok=True)
    sys.exit(fig_throat(d / "fig_throat.png") + fig_thermal(d / "fig_thermal.png"))
