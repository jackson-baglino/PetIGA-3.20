#!/usr/bin/env python3
"""
Publication-style schematic of the ice/air double-well potential

    F^dub(phi) = 1/2 * phi^2 * (1 - phi)^2

(the same normalization used in plotTripleWell.py; its derivative is the
f1(phi) term assembled in src/assembly.c::DoubleWellDeriv). phi = phi_i is
the ice phase fraction, so the well at phi=0 is pure air (phi_a = 1) and the
well at phi=1 is pure ice.

Produces two vector (SVG) versions, one for light slide backgrounds and one
for dark slide backgrounds, both with a transparent canvas so they drop onto
either directly.

All math (phi_i, F^dub, and the axis tick values) is set in mathtext's
Computer-Modern fontset for a consistent LaTeX look. That font's symbol
glyphs (e.g. \\phi) are NOT at standard Unicode codepoints, so they render
correctly only where the same CM math font is installed -- which Inkscape
does not ship by default. To stay portable, text is therefore baked to
vector paths (svg.fonttype="path") rather than left as live, re-typeable
text: the glyph shapes are then guaranteed correct everywhere, including
in Inkscape. Resizing (the whole figure, or a selected label) still works
normally in Inkscape; only retyping/changing the font family of a single
text run does not.
"""

import numpy as np
import matplotlib.pyplot as plt

FIGSIZE = (3.5, 4.0)        # inches
FONT_TICK = 18               # smallest font on the slide, per spec
FONT_LABEL = 22
FONT_CURVE = 22
FONT_ANNOT = 20

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",   # Computer-Modern-style math, matches LaTeX
    "svg.fonttype": "path",     # bake glyphs to paths -- see module docstring
})

THEMES = {
    "light": dict(curve="#0B4F6C", accent="#C0392B", ink="#111111"),
    "dark":  dict(curve="#7FD8FF", accent="#FF8C42", ink="#F2F2F2"),
}


def double_well(phi):
    return 0.5 * phi**2 * (1.0 - phi) ** 2


def make_plot(theme_name, out_path):
    c = THEMES[theme_name]
    phi = np.linspace(-0.18, 1.18, 600)
    F = double_well(phi)
    peak = double_well(0.5)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    ax.plot(phi, F, color=c["curve"], linewidth=2.6, solid_capstyle="round")

    # Zeros of the well: phi=0 (air) and phi=1 (ice)
    ax.plot([0, 1], [0, 0], "o", color=c["accent"], markersize=8,
             zorder=5, clip_on=False)

    ax.annotate("air", xy=(0, 0), xytext=(0, -26),
                textcoords="offset points", ha="center", va="top",
                fontsize=FONT_ANNOT, color=c["ink"], annotation_clip=False)
    ax.annotate("ice", xy=(1, 0), xytext=(0, -26),
                textcoords="offset points", ha="center", va="top",
                fontsize=FONT_ANNOT, color=c["ink"], annotation_clip=False)

    # Curve label above the central barrier
    ax.annotate(r"$F^{\mathrm{dub}}$", xy=(0.5, peak), xytext=(0, 16),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=FONT_CURVE, color=c["curve"])

    # phi_i: the ice phase fraction (assembly.c: `ice`); \phi_i throughout,
    # matching docs/model_description.md (never \varphi).
    ax.set_xlabel(r"$\phi_i$", fontsize=FONT_LABEL, color=c["ink"])
    ax.set_ylabel(r"$F^{\mathrm{dub}}(\phi_i)$",
                  fontsize=FONT_LABEL, color=c["ink"])

    ax.set_xlim(phi.min(), phi.max())
    ax.set_ylim(0.0, peak * 1.55)

    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels([r"$0$", r"$0.5$", r"$1$"],
                        fontsize=FONT_TICK, color=c["ink"])
    ax.set_yticks([])

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(c["ink"])
        ax.spines[side].set_linewidth(1.2)

    ax.tick_params(axis="x", colors=c["ink"], length=4, width=1.2)
    ax.tick_params(axis="y", length=0)

    fig.subplots_adjust(left=0.16, right=0.95, top=0.88, bottom=0.22)
    fig.savefig(out_path, transparent=True)
    plt.close(fig)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    make_plot("light", "postprocess/figures/double_well_light.svg")
    make_plot("dark", "postprocess/figures/double_well_dark.svg")
