"""Plot the double-well potential F^duo(phi_i) / gamma_ia for the AGU 2026 talk.

The slide (AGU 2026 deck, slide 5) writes the well in the calibrated
physical form of docs/contact_angle/contact_angle.tex (eq. 1),

    F^duo(phi_i) = gamma_ia * phi_i^2 * (1 - phi_i)^2,

with gradient term 3 eps gamma_ia |grad phi_i|^2, so the ice-vapor interface
carries energy gamma_ia. The plot shows F^duo / gamma_ia = A phi^2 (1 - phi)^2
with A = 1 by default; the barrier at phi_i = 1/2 is A/16. (The solver's
DoubleWellDeriv uses the half well because its mobility M = 2 gamma_ia L
absorbs the factor; --prefactor 0.5 plots that scaled form.)

The figure is drawn at the size it occupies on the slide (--width, --height
in inches), so every font size below is the size the audience sees. Colors
are the deck's "Slate" theme (ppt/theme/theme1.xml). Text uses the theme
body font, Franklin Gothic Book, taken from the PowerPoint bundle; math uses
Computer Modern to match the IguanaTex equations.

    venv_cryosuction/bin/python docs/figures/plot_double_well.py --out <dir>/double_well.png
"""
from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "postprocess"))
from figcheck import check_figure  # noqa: E402

# Slate theme colors (ppt/theme/theme1.xml of the deck)
INK = "#1C2229"      # dk1, text
INK2 = "#2E3640"     # dk2, tick labels
MUTED = "#66707C"    # accent3, axes and guide lines
CRIMSON = "#8A1F2E"  # accent1, the curve (same as the equation labels)

DEFAULT_FONT = Path("/Applications/Microsoft PowerPoint.app/Contents/Resources/"
                    "DFonts/Franklin Gothic Book.ttf")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True, help="output PNG path")
    ap.add_argument("--width", type=float, default=4.12, help="placed width [in]")
    ap.add_argument("--height", type=float, default=2.89, help="placed height [in]")
    ap.add_argument("--prefactor", type=float, default=1.0,
                    help="A in F/gamma = A phi^2 (1-phi)^2 (slide: 1; solver's scaled form: 0.5)")
    ap.add_argument("--font", type=Path, default=DEFAULT_FONT, help="text font file")
    ap.add_argument("--dpi", type=int, default=400)
    ap.add_argument("--compact", action="store_true",
                    help="short placements: horizontal y label above the axis, phi_i at the axis end, "
                         "air/ice written above the minima")
    args = ap.parse_args()

    if args.font.exists():
        font_manager.fontManager.addfont(str(args.font))
        family = font_manager.FontProperties(fname=str(args.font)).get_name()
    else:
        family = "DejaVu Sans"
        print(f"warning: {args.font} not found, using {family}")

    plt.rcParams.update({
        "font.family": family,
        "mathtext.fontset": "cm",
        "font.size": 16,
        "axes.labelsize": 18,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 1.2,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
    })

    A = args.prefactor
    barrier = A / 16.0                      # value at phi = 1/2
    phi = np.linspace(-0.25, 1.25, 601)
    F = A * phi**2 * (1.0 - phi) ** 2

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    fig.patch.set_alpha(0.0)
    ax.set_facecolor("none")

    ax.plot(phi, F, color=CRIMSON, lw=3.0, zorder=3, solid_capstyle="round")
    ax.plot([0.0, 1.0], [0.0, 0.0], "o", ms=9, color=CRIMSON,
            mec="white", mew=1.5, zorder=4)

    ax.set_xlim(phi[0], phi[-1])
    ax.set_ylim(-0.08 * 1.6 * barrier, 1.6 * barrier)
    fb = Fraction(A).limit_denominator(1000) / 16
    frac = f"{fb.numerator}/{fb.denominator}" if fb.numerator != fb.denominator else "1"
    ax.set_yticks([0.0, barrier], ["0", f"${frac}$"])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if args.compact:
        ax.set_xticks([0.0, 0.5, 1.0], ["0", "0.5", "1"])
        for x0, name in ((0.0, "air"), (1.0, "ice")):
            ax.text(x0, 0.30 * barrier, name, color=INK2, fontsize=15, ha="center", va="bottom")
        ax.set_title(r"$F^{\mathrm{duo}}/\gamma_{ia}$", color=INK, fontsize=18, loc="left", pad=4)
        ax.text(1.0, -0.035, r"$\phi_i$", transform=ax.transAxes, color=INK, fontsize=18,
                ha="right", va="top")
    else:
        ax.set_xticks([0.0, 0.5, 1.0], ["0\nair", "0.5", "1\nice"])
        ax.set_xlabel(r"$\phi_i$", color=INK, labelpad=2)
        ax.set_ylabel(r"$F^{\mathrm{duo}}/\gamma_{ia}$", color=INK, labelpad=4)

    fig.tight_layout(pad=0.3)

    problems = check_figure(fig, min_pt=14.0)
    # figcheck does not test clipping: every label must lie inside the figure
    rnd = fig.canvas.get_renderer()
    fb_box = fig.bbox
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text().strip():
            bb = t.get_window_extent(rnd)
            if bb.x0 < fb_box.x0 - 0.5 or bb.y0 < fb_box.y0 - 0.5 or bb.x1 > fb_box.x1 + 0.5 or bb.y1 > fb_box.y1 + 0.5:
                problems.append(f"clipped: {t.get_text()[:30]!r}")
    for p in problems:
        print("figcheck:", p)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi, transparent=True)
    print(f"font: {family}; size {args.width} x {args.height} in; prefactor A = {A:g}; "
          f"barrier = {barrier:g} gamma_ia; dpi {args.dpi}")
    print(f"figcheck (min 14 pt): {'clean' if not problems else f'{len(problems)} problem(s)'}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
