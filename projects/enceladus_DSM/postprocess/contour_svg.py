#!/usr/bin/env python3
"""Export a phase-field level set as a vector SVG for figure work.

Why this exists: ParaView's "Export Scene -> SVG" rasterises the surface and
buries the geometry in thousands of clipped polygons, and matplotlib's SVG
output wraps every artist in its own clip path and transform. Neither is
pleasant to edit in Inkscape. This writes the contour directly: one <path>
per connected level-set loop, coordinates in micrometres, nothing else in the
file except the optional domain box and symmetry axis, each in its own named
group (Inkscape shows them as layers).

    python postprocess/contour_svg.py <run_dir> --mirror -o fig.svg

By default it takes the LAST snapshot in <run_dir>/vtkOut, contours IcePhase
at 0.5, crops to the ice bounding box, and mirrors about y = 0 for axisymmetric
runs (--mirror), which is how a 2D axisymmetric grain pair is normally drawn.

The viewBox is in micrometres and width/height are set in millimetres, so
1 user unit = 1 um and the figure opens in Inkscape at its intended print size.
Stroke weight is specified in points and converted so it comes out at that
weight on the page.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
from contourpy import contour_generator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pplib import read_vts, step_of, step_times  # noqa: E402

MM_PER_PT = 25.4 / 72.0


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def rdp(pts: np.ndarray, tol: float) -> np.ndarray:
    """Ramer-Douglas-Peucker, iterative so a 10k-vertex loop cannot blow the
    recursion limit."""
    if tol <= 0 or len(pts) < 3:
        return pts

    keep = np.zeros(len(pts), dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        seg = pts[j] - pts[i]
        norm = np.hypot(*seg)
        rel = pts[i + 1:j] - pts[i]
        if norm == 0.0:
            d = np.hypot(rel[:, 0], rel[:, 1])
        else:
            d = np.abs(rel[:, 0] * seg[1] - rel[:, 1] * seg[0]) / norm
        k = int(np.argmax(d))
        if d[k] > tol:
            k += i + 1
            keep[k] = True
            stack.append((i, k))
            stack.append((k, j))
    return pts[keep]


def close_on_axis(line: np.ndarray, mirror: bool, y_tol: float):
    """Turn one contour polyline into a closed loop.

    An axisymmetric ice body is cut by the symmetry plane, so its contour comes
    back OPEN with both endpoints sitting on y = 0. Mirroring reflects the
    profile and joins the two halves into a single genuine outline; without
    --mirror the loop is closed with a straight segment along the axis, which
    is the correct half-plane silhouette.
    """
    closed = np.allclose(line[0], line[-1])
    if closed:
        loops = [line[:-1]]
        if mirror and np.abs(line[:, 1]).min() > y_tol:
            flip = line[:-1].copy()
            flip[:, 1] *= -1.0
            loops.append(flip[::-1])
        return loops, True

    on_axis = abs(line[0, 1]) <= y_tol and abs(line[-1, 1]) <= y_tol
    if mirror and on_axis:
        flip = line[::-1].copy()
        flip[:, 1] *= -1.0
        return [np.vstack([line, flip[1:-1]])], True
    return [line], on_axis


# ---------------------------------------------------------------------------
# svg
# ---------------------------------------------------------------------------
def path_d(loop: np.ndarray, y0: float, digits: int) -> str:
    """One closed sub-path, y flipped so +y is up on the page."""
    out = []
    for k, (px, py) in enumerate(loop):
        out.append(f"{'M' if k == 0 else 'L'}{px:.{digits}f},{y0 - py:.{digits}f}")
    out.append("Z")
    return " ".join(out)


def write_svg(path, loops, bbox, args, meta):
    x0, x1, y0, y1 = bbox
    w_um, h_um = x1 - x0, y1 - y0
    width_mm = args.width_mm
    height_mm = width_mm * h_um / w_um
    um_per_mm = w_um / width_mm
    stroke_um = args.stroke_pt * MM_PER_PT * um_per_mm

    d = " ".join(path_d(lp, y1, args.digits) for lp in loops)
    fill = "none" if args.no_fill else args.fill

    body = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
        f'width="{width_mm:.4f}mm" height="{height_mm:.4f}mm" '
        f'viewBox="{x0:.{args.digits}f} 0 {w_um:.{args.digits}f} {h_um:.{args.digits}f}">',
        "  <!-- " + " | ".join(meta) + " -->",
        f'  <g id="ice" fill="{fill}" fill-rule="evenodd" stroke="{args.stroke}" '
        f'stroke-width="{stroke_um:.5g}" stroke-linejoin="round">',
        f'    <path id="phi_{args.level:g}" d="{d}"/>',
        "  </g>",
    ]
    if args.domain:
        dx0, dx1, dy0, dy1 = meta_box = args.domain_box
        body += [
            f'  <g id="domain" fill="none" stroke="{args.stroke}" '
            f'stroke-width="{stroke_um:.5g}" stroke-dasharray="{4*stroke_um:.5g},{3*stroke_um:.5g}">',
            f'    <rect x="{dx0:.{args.digits}f}" y="{y1 - dy1:.{args.digits}f}" '
            f'width="{dx1 - dx0:.{args.digits}f}" height="{dy1 - dy0:.{args.digits}f}"/>',
            "  </g>",
        ]
        del meta_box
    if args.axis:
        body += [
            f'  <g id="symmetry-axis" stroke="{args.stroke}" '
            f'stroke-width="{stroke_um:.5g}" stroke-dasharray="{6*stroke_um:.5g},{3*stroke_um:.5g},{stroke_um:.5g},{3*stroke_um:.5g}">',
            f'    <line x1="{x0:.{args.digits}f}" y1="{y1:.{args.digits}f}" '
            f'x2="{x1:.{args.digits}f}" y2="{y1:.{args.digits}f}"/>',
            "  </g>",
        ]
    body.append("</svg>")
    with open(path, "w") as fh:
        fh.write("\n".join(body) + "\n")


def write_preview(path, loops, bbox, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    x0, x1, y0, y1 = bbox
    fig, ax = plt.subplots(figsize=(7, 7 * (y1 - y0) / (x1 - x0)))
    for lp in loops:
        ax.add_patch(Polygon(lp, closed=True,
                             facecolor="none" if args.no_fill else args.fill,
                             edgecolor=args.stroke, lw=1.0))
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("z [um]")
    ax.set_ylabel("r [um]")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
def resolve_snapshot(target, step):
    if os.path.isfile(target):
        return target
    vtk_dir = os.path.join(target, "vtkOut")
    if not os.path.isdir(vtk_dir):
        vtk_dir = target
    files = sorted(glob.glob(os.path.join(vtk_dir, "solV_*.vts")), key=step_of)
    if not files:
        sys.exit(f"no solV_*.vts under {vtk_dir}")
    if step is None:
        return files[-1]
    match = [f for f in files if step_of(f) == step]
    if not match:
        sys.exit(f"step {step} not found; have {step_of(files[0])}..{step_of(files[-1])}")
    return match[0]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("target", help="run directory, vtkOut directory, or a solV_*.vts file")
    p.add_argument("-o", "--out", required=True, help="output .svg path")
    p.add_argument("--step", type=int, default=None, help="step number (default: last)")
    p.add_argument("--field", default="IcePhase")
    p.add_argument("--level", type=float, default=0.5)
    p.add_argument("--mirror", action="store_true",
                   help="reflect about y = 0 (axisymmetric runs)")
    p.add_argument("--crop", choices=["contour", "domain"], default="contour")
    p.add_argument("--margin-um", type=float, default=5.0)
    p.add_argument("--width-mm", type=float, default=120.0,
                   help="printed width of the drawing [mm]")
    p.add_argument("--stroke-pt", type=float, default=1.0,
                   help="outline weight on the page [pt]")
    p.add_argument("--stroke", default="#000000")
    p.add_argument("--fill", default="#d9d9d9")
    p.add_argument("--no-fill", action="store_true")
    p.add_argument("--simplify-um", type=float, default=0.02,
                   help="RDP tolerance [um]; 0 disables. The default is ~1/20 of "
                        "a grid cell, i.e. far below anything visible.")
    p.add_argument("--digits", type=int, default=4)
    p.add_argument("--domain", action="store_true", help="draw the domain box")
    p.add_argument("--axis", action="store_true", help="draw the symmetry axis")
    p.add_argument("--no-preview", action="store_true", help="skip the PNG preview")
    args = p.parse_args()

    fn = resolve_snapshot(args.target, args.step)
    fields, X, Y = read_vts(fn, want=[args.field])
    if args.field not in fields:
        sys.exit(f"{fn} has no point-data array {args.field!r}")

    um = 1e6
    x, y = X[0, :] * um, Y[:, 0] * um
    z = fields[args.field]
    dx = float(np.min(np.diff(x)))

    lines = contour_generator(x=x, y=y, z=z).lines(args.level)
    if not lines:
        sys.exit(f"no {args.field} = {args.level} contour in {os.path.basename(fn)}")

    loops, n_raw, open_ended = [], 0, 0
    for line in lines:
        line = np.asarray(line, dtype=float)
        n_raw += len(line)
        parts, closed_ok = close_on_axis(line, args.mirror, y_tol=0.5 * dx)
        open_ended += 0 if closed_ok else 1
        for lp in parts:
            lp = rdp(lp, args.simplify_um)
            if len(lp) >= 3:
                loops.append(lp)
    if open_ended:
        print(f"note: {open_ended} contour(s) ran off a non-axis boundary and were "
              f"closed with a straight segment")

    allpts = np.vstack(loops)
    dom = (float(x.min()), float(x.max()),
           -float(y.max()) if args.mirror else float(y.min()), float(y.max()))
    if args.crop == "domain":
        bbox = dom
    else:
        m = args.margin_um
        bbox = (allpts[:, 0].min() - m, allpts[:, 0].max() + m,
                allpts[:, 1].min() - m, allpts[:, 1].max() + m)
    args.domain_box = dom

    tmap = step_times(args.target if os.path.isdir(args.target) else
                      os.path.dirname(os.path.dirname(fn)))
    t = tmap.get(step_of(fn))
    meta = [os.path.basename(fn),
            f"step {step_of(fn)}",
            f"t = {t:.6g} s" if t is not None else "t = unknown",
            f"{args.field} = {args.level:g}",
            "mirrored about y=0" if args.mirror else "half-plane",
            f"units: um; simplify {args.simplify_um:g} um",
            f"page y = {bbox[3]:.4f} um - r"]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    write_svg(args.out, loops, bbox, args, meta)
    n_out = sum(len(lp) for lp in loops)
    print(f"{args.out}: {len(loops)} loop(s), {n_out} vertices "
          f"(from {n_raw}), bbox [{bbox[0]:.2f}, {bbox[1]:.2f}] x "
          f"[{bbox[2]:.2f}, {bbox[3]:.2f}] um")
    print("  " + " | ".join(meta))

    if not args.no_preview:
        png = os.path.splitext(args.out)[0] + "_preview.png"
        write_preview(png, loops, bbox, args)
        print(f"{png}: preview")


if __name__ == "__main__":
    main()
