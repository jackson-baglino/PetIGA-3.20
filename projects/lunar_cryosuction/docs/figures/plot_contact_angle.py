"""Contact-angle diagram for the AGU 2026 talk (wall boundary condition slide).

An ice cap sits on a regolith wall Gamma. The ice is drawn as the model's own
diffuse field, phi = 1 / (1 + exp(-d / eps)) in the signed distance d to a
circular cap (the 1D equilibrium profile, docs/interface_width_conventions.md),
so the phi = 1/2 contour is the level set the contact angle is measured on.
theta is measured through the ice, as in the solver's wall condition
n . grad(phi) = cos(theta) phi (1 - phi) / eps with cos(theta) =
(gamma_as - gamma_is) / gamma_ia (docs/contact_angle/README.md); n is the
outward normal of the domain, which points into the wall.

The figure is drawn at its placed slide size, so font sizes are what the
audience sees. Colors are the deck's Slate theme; text is Franklin Gothic Book,
math is Computer Modern (as in the IguanaTex equations).

    venv_cryosuction/bin/python docs/figures/plot_contact_angle.py --out <dir>/contact_angle.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import ListedColormap, to_rgb  # noqa: E402
from matplotlib.patches import Arc, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "postprocess"))
from figcheck import check_figure  # noqa: E402

INK = "#1C2229"      # dk1
INK2 = "#2E3640"     # dk2
MUTED = "#66707C"    # accent3
CRIMSON = "#8A1F2E"  # accent1
SLATE = "#3F5F7A"    # accent2
ICE = "#9CC5DC"      # accent4

DEFAULT_FONT = Path("/Applications/Microsoft PowerPoint.app/Contents/Resources/"
                    "DFonts/Franklin Gothic Book.ttf")


def main() -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--width", type=float, default=4.2, help="placed width [in]")
    ap.add_argument("--height", type=float, default=1.6, help="placed height [in]")
    ap.add_argument("--theta", type=float, default=60.0, help="contact angle to draw [deg]")
    ap.add_argument("--font", type=Path, default=DEFAULT_FONT)
    ap.add_argument("--fontsize", type=float, default=15.0)
    ap.add_argument("--dpi", type=int, default=400)
    args = ap.parse_args()

    if args.font.exists():
        font_manager.fontManager.addfont(str(args.font))
        family = font_manager.FontProperties(fname=str(args.font)).get_name()
    else:
        family = "DejaVu Sans"
        print(f"warning: {args.font} not found, using {family}")
    plt.rcParams.update({"font.family": family, "mathtext.fontset": "cm",
                         "font.size": args.fontsize})
    fs = args.fontsize

    # geometry in inches = data units (equal aspect, axes fill the figure)
    W, H = args.width, args.height
    wall = 0.36                         # wall band thickness
    th = np.radians(args.theta)
    half = 0.95                         # half width of the cap footprint
    R = half / np.sin(th)
    xc, yc = 0.0, wall - R * np.cos(th)  # circle center (below the wall surface for theta < 90)
    eps = R / 60.0                      # thin diffuse band, close to production eps/R

    fig = plt.figure(figsize=(W, H))
    fig.patch.set_alpha(0.0)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(-W / 2, W / 2)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")

    # diffuse ice field above the wall: alpha ramps with phi
    x = np.linspace(-W / 2, W / 2, 900)
    y = np.linspace(wall, H, 400)
    X, Y = np.meshgrid(x, y)  # noqa: N806
    d = R - np.hypot(X - xc, Y - yc)
    phi = 1.0 / (1.0 + np.exp(-d / eps))
    rgba = np.zeros((256, 4))
    rgba[:, :3] = to_rgb(ICE)
    rgba[:, 3] = np.linspace(0, 1, 256)
    field = ax.imshow(phi, extent=(-W / 2, W / 2, wall, H), origin="lower",
                      cmap=ListedColormap(rgba), vmin=0, vmax=1, interpolation="bilinear", zorder=1)
    # phi = 1/2 is exactly d = 0, the circle of radius R: draw it as a line
    a = np.linspace(np.pi / 2 - th, np.pi / 2 + th, 300)
    ax.plot(xc + R * np.cos(a), yc + R * np.sin(a), color=SLATE, lw=2.0, zorder=3)

    # regolith wall
    # wall fill is background for the labels written on it (figcheck: overlap_ok)
    ax.add_patch(Rectangle((-W / 2, 0), W, wall, facecolor=MUTED, alpha=0.25,
                           edgecolor="none", zorder=2, gid="overlap_ok"))
    ax.plot([-W / 2, W / 2], [wall, wall], color=INK2, lw=1.6, zorder=4)

    # contact angle at the right contact line, measured through the ice
    cx = xc + R * np.sin(th)
    t_len = 0.62
    tx, ty = cx - t_len * np.cos(th), wall + t_len * np.sin(th)
    ax.plot([cx, tx], [wall, ty], color=INK2, lw=1.2, ls=(0, (3, 2)), zorder=4)
    r_arc = 0.34
    ax.add_patch(Arc((cx, wall), 2 * r_arc, 2 * r_arc, theta1=180 - args.theta, theta2=180,
                     color=CRIMSON, lw=2.2, zorder=5))
    am = np.radians(180 - args.theta / 2)
    ax.text(cx + (r_arc + 0.12) * np.cos(am), wall + (r_arc + 0.10) * np.sin(am), r"$\theta$",
            color=CRIMSON, fontsize=fs + 3, ha="center", va="center", zorder=6)

    # outward normal n points into the wall
    nx = W / 2 - 0.32
    ax.arrow(nx, wall, 0, -0.29, width=0.014, head_width=0.09, head_length=0.08,
             length_includes_head=True, color=INK, zorder=6)
    ax.text(nx + 0.08, wall - 0.16, r"$\mathbf{n}$", color=INK, fontsize=fs + 1,
            ha="left", va="center", zorder=6)

    # the three interface energies where they act
    am_ia = np.radians(90 + 0.62 * args.theta)
    ax.text(xc + (R + 0.07) * np.cos(am_ia), yc + (R + 0.07) * np.sin(am_ia), r"$\gamma_{ia}$",
            color=SLATE, fontsize=fs + 1, ha="right", va="bottom", zorder=6)
    ax.text(xc, wall / 2, r"$\gamma_{is}$", color=INK, fontsize=fs + 1, ha="center", va="center", zorder=6)
    ax.text((cx + nx) / 2 - 0.05, wall / 2, r"$\gamma_{as}$", color=INK, fontsize=fs + 1,
            ha="center", va="center", zorder=6)

    # labels
    ax.text(xc, wall + 0.17, "ice", color=INK, ha="center", va="center", fontsize=fs, zorder=6)
    ax.text(-W / 2 + 0.12, H - 0.16, "vapor", color=INK2, ha="left", va="top", fontsize=fs, zorder=6)
    ax.text(-W / 2 + 0.10, wall / 2, r"regolith $\Gamma$", color=INK, ha="left", va="center",
            fontsize=fs, zorder=6)

    # the ice field is the background every label sits on; check the labels
    # against all lines, arcs and arrows with it hidden, then restore it
    field.set_visible(False)
    problems = check_figure(fig, min_pt=14.0)
    field.set_visible(True)
    for p in problems:
        print("figcheck:", p)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi, transparent=True)
    print(f"font: {family}; size {W} x {H} in; theta {args.theta:g} deg; eps/R = 1/60; dpi {args.dpi}")
    print(f"figcheck (min 14 pt): {'clean' if not problems else f'{len(problems)} problem(s)'}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
