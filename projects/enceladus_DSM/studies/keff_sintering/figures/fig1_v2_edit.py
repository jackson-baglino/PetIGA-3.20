#!/usr/bin/env python3
"""Figure 1, version 2 (advisor meeting 2026-10-15): edits applied ON TOP of the
hand-edited Figure1_mechanisms_aggregate.pdf, which stays untouched.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig1_v2_edit.py
        [--src <Figure1_mechanisms_aggregate.pdf>] [--out <dir>]

  * the w(t) arrow (axis -> upper neck surface) becomes an I-bar across the
    whole neck, upper to lower surface: w is the full waist width
    (postprocess/neck_width.py), and the bar reads as a measurement;
  * water molecules (dots) leave the small grain's surface along the
    "Vapor Transport" arrow and gather at the neck; the arrow is shortened by
    3.7 pt so its head stops below the I-bar's lower cap.

The PDF is round-tripped through Inkscape's SVG (text there is glyph outlines,
so nothing reflows). Coordinates are the SVG's (pt, y down) and were read off
the 2026-10-14 PDF: neck surfaces at y = 55.27 / 95.46 on the axis x = 219.59.
Writes Figure1_mechanisms_aggregate_v2.{pdf,png}. --dots N R (in flight,
radius pt) draws N molecules of radius R instead of the hand-placed ten, with
about N/3 deposited along the neck's lower surface, and adds --suffix to the name.
"""
from __future__ import annotations

import argparse, re, subprocess, tempfile
from pathlib import Path

import numpy as np

INKSCAPE = "/Applications/Inkscape.app/Contents/MacOS/inkscape"
SAMPLES = Path.home() / "SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign/compare/figure_samples"
X_AX, Y_TOP, Y_BOT = 219.586, 55.27, 95.46       # neck axis and surfaces
SHIFT = 3.7                                       # vapour arrow shortened by this [pt]
P0, P2, P3 = (171.031, 126.371), (222.258, 163.617), (219.6875, 99.859 + SHIFT)
MOL = "#2b4c9b"                                   # darker than the cyan initial-grain outline


def bezier(t):
    """The vapour arrow's cubic (first control point = start point)."""
    p0, p1, p2, p3 = (np.array(p) for p in (P0, P0, P2, P3))
    t = np.asarray(t)[:, None]
    pt = (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3
    dt = 3 * (1 - t) ** 2 * (p1 - p0) + 6 * (1 - t) * t * (p2 - p1) + 3 * t ** 2 * (p3 - p2)
    n = np.c_[-dt[:, 1], dt[:, 0]] / np.hypot(dt[:, 0], dt[:, 1])[:, None]
    return pt, n


def spread_dots(n, r, outline, seed=3):
    """n molecules in flight along the arrow, jittered off it, and ~n/3 more
    deposited just outside the neck's lower surface (from the outline)."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.02, 0.9, n) + rng.uniform(-0.25, 0.25, n) * 0.88 / n
    side = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    off = side * rng.uniform(0.9, 2.6, n)
    pt, nrm = bezier(np.clip(t, 0.0, 0.93))
    fly = pt + off[:, None] * nrm
    # Lower neck surface, either side of the arrow head; outward = toward the
    # pore below the neck.
    o = outline[(outline[:, 1] > 95.0) & (outline[:, 1] < 108.0) & (np.abs(outline[:, 0] - X_AX) < 13.0)]
    o = o[np.argsort(o[:, 0])]
    d = np.gradient(o, axis=0)
    nn = np.c_[-d[:, 1], d[:, 0]] / np.hypot(d[:, 0], d[:, 1])[:, None]
    nn *= np.sign(np.sum(nn * (np.array([X_AX, 125.0]) - o), axis=1))[:, None]
    m = max(2, n // 3)
    xs = np.r_[np.linspace(X_AX - 11.5, X_AX - 3.4, m // 2), np.linspace(X_AX + 3.4, X_AX + 11.5, m - m // 2)]
    dep = []
    for x in xs:
        k = int(np.argmin(np.abs(o[:, 0] - x)))
        dep.append(o[k] + (r + 0.9) * nn[k])
    return np.vstack([fly, dep])


def overlay(dots=None, outline=None):
    """SVG for the I-bar and the molecules."""
    cap, lw = 2.6, 1.0
    g = [f'<g id="v2_ibar" stroke="#000000" stroke-width="{lw}" stroke-linecap="butt" fill="none">',
         f'<path d="M {X_AX} {Y_TOP} L {X_AX} {Y_BOT}"/>',
         f'<path d="M {X_AX - cap} {Y_TOP} L {X_AX + cap} {Y_TOP}"/>',
         f'<path d="M {X_AX - cap} {Y_BOT} L {X_AX + cap} {Y_BOT}"/></g>']
    # In flight: spaced along the arrow, nudged off it alternately so they read
    # as molecules, not a dotted line; denser where they leave the surface.
    t = np.array([0.03, 0.09, 0.16, 0.25, 0.36, 0.48, 0.60, 0.71, 0.81, 0.89])
    off = np.array([2.2, -1.6, 2.0, -1.4, 2.2, -1.2, 2.2, -1.6, 2.2, -1.8])
    pt, n = bezier(t)
    xy = pt + off[:, None] * n
    # Arrived: against the neck's lower surface, either side of the arrow head.
    xy = np.vstack([xy, [[215.2, 98.4], [223.9, 98.6], [212.6, 101.6], [226.6, 101.9]]])
    rad = 1.25
    if dots is not None:
        xy, rad = spread_dots(int(dots[0]), dots[1], outline), dots[1]
    g.append(f'<g id="v2_molecules" fill="{MOL}" stroke="none">')
    g += [f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{rad:g}"/>' for x, y in xy]
    g.append("</g>")
    return "\n".join(g)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", type=Path, default=SAMPLES / "Figure1_mechanisms_aggregate.pdf")
    ap.add_argument("--out", type=Path, default=SAMPLES)
    ap.add_argument("--dots", type=float, nargs=2, default=None, metavar=("N", "R"))
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        svg = Path(td) / "f.svg"
        subprocess.run([INKSCAPE, str(a.src), "--pdf-poppler=false", "--export-type=svg",
                        f"--export-filename={svg}"], check=True, capture_output=True)
        s = svg.read_text()

        def element(pid):
            i = s.find(f'id="{pid}"')
            return s.rfind("<path", 0, i), s.find("/>", i) + 2

        # the w(t) arrow: shaft path44, head path45 -- removed
        for pid in ("path45", "path44"):
            i, j = element(pid)
            assert "219.585938" in s[i:j], f"{pid} is not the w(t) arrow"
            s = s[:i] + s[j:]
        # the vapour arrow: shaft path28 ends SHIFT lower, head path29 moves with it
        i, j = element("path28")
        s = s[:i] + s[i:j].replace("219.6875 99.859375", f"219.6875 {99.859375 + SHIFT:.6f}") + s[j:]
        i, j = element("path29")
        s = s[:i] + f'<g transform="translate(0,{SHIFT})">' + s[i:j] + "</g>" + s[j:]
        i, j = element("path25")                 # the grain-pair outline
        d = re.search(r' d="([^"]*)"', s[i:j]).group(1)
        outline = np.array([float(v) for v in re.findall(r"-?\d+\.?\d*", d)]).reshape(-1, 2)
        s = s.replace("</svg>", overlay(a.dots, outline) + "\n</svg>")
        svg.write_text(s)
        pdf = a.out / f"Figure1_mechanisms_aggregate_v2{a.suffix}.pdf"
        subprocess.run([INKSCAPE, str(svg), "--export-type=pdf", f"--export-filename={pdf}"],
                       check=True, capture_output=True)
    png = pdf.with_suffix(".png")
    subprocess.run(["mutool", "draw", "-q", "-c", "rgba", "-r", "600", "-o", str(png), str(pdf), "1"], check=True)
    print(f"wrote {pdf} / .png")


if __name__ == "__main__":
    main()
