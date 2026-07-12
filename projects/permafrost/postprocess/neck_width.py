#!/usr/bin/env python3
"""Extract neck width vs time from a two-grain sintering run.

Method: MINIMUM CROSS-SECTION, not curvature tracking.

For each snapshot, the ice body's vertical chord width w(x) is measured on
every mesh column (the span between the outermost phi = 0.5 crossings,
sub-cell interpolated). Restricted to the region strictly between the two
grain-center columns, w(x) has a single interior minimum — the neck plane —
and the neck width is min w(x). This definition:

  * tracks the neck plane automatically as it migrates toward the smaller
    grain (the minimum moves with it; nothing is anchored to the t = 0
    location);
  * matches the experimental convention (Molaro et al. report the neck
    WIDTH, i.e. the full waist of the bonded pair);
  * needs no derivatives of the contour. The proposed alternative — locate
    the two curvature maxima flanking the neck — identifies the same plane
    but requires second differences of a discretized contour, which is
    noise-amplifying at exactly the saddle points where the contour is
    sparsest. The waist minimum is the integral-robust version of the same
    idea. (If the rim points themselves are ever needed, e.g. for rim-radius
    measurements, add them as a separate diagnostic then.)

Caveat: if a run coarsens so far that the waist disappears (single convex
body), the interior minimum merges with a grain shoulder; the script reports
the minimum location so that regime is visible in the output.

Usage:
    python neck_width.py <run_dir> [--phi 0.5] [--out neck_width.csv]

Reads <run_dir>/vtkOut/solV_*.vts (snapshot index = timestep, -outp 1) and
<run_dir>/outp.txt (step -> time map from the monitor tables). Writes a CSV
(t, neck_width, x_neck) and a PNG plot next to it.
"""

import argparse
import base64
import glob
import re
import struct
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


def read_vts_phi(fn):
    """Return (phi[ny,nx], x[nx], y[ny]) from a solV_*.vts snapshot."""
    root = ET.parse(fn).getroot()
    grid = root.find(".//StructuredGrid")
    ext = [int(v) for v in grid.get("WholeExtent").split()]
    nx, ny = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1

    def decode(da):
        raw = base64.b64decode("".join(da.text.split()))
        n = struct.unpack("<Q", raw[:8])[0]
        return np.frombuffer(raw[8:8 + n], dtype=np.float64)

    pts = None
    for da in root.findall(".//Points/DataArray"):
        pts = decode(da).reshape(ny, nx, 3)
    phi = None
    for da in root.findall(".//PointData/DataArray"):
        if da.get("Name") == "IcePhase":
            phi = decode(da).reshape(ny, nx)
    if phi is None or pts is None:
        raise RuntimeError(f"IcePhase/Points missing in {fn}")
    return phi, pts[0, :, 0], pts[:, 0, 1]


def chord_width(col, y, level):
    """Full span between outermost `level`-crossings of phi along a column."""
    above = col >= level
    if not above.any():
        return 0.0
    idx = np.flatnonzero(above)
    lo_i, hi_i = idx[0], idx[-1]
    # sub-cell interpolation at both ends
    y_lo = y[lo_i]
    if lo_i > 0:
        f = (level - col[lo_i - 1]) / (col[lo_i] - col[lo_i - 1])
        y_lo = y[lo_i - 1] + f * (y[lo_i] - y[lo_i - 1])
    y_hi = y[hi_i]
    if hi_i < len(y) - 1:
        f = (col[hi_i] - level) / (col[hi_i] - col[hi_i + 1])
        y_hi = y[hi_i] + f * (y[hi_i + 1] - y[hi_i])
    return y_hi - y_lo


def step_times(outp):
    """Map step -> time from the monitor tables in outp.txt."""
    tmap = {}
    pat = re.compile(r"^\s+(\d+)\s+\|\s+([0-9.eE+-]+)\s+\|")
    for line in open(outp, errors="replace"):
        m = pat.match(line)
        if m and line.count("|") == 8:
            tmap[int(m.group(1))] = float(m.group(2))
    return tmap


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--phi", type=float, default=0.5, help="contour level")
    ap.add_argument("--out", type=Path, default=None, help="CSV path")
    ap.add_argument("--axisym", action="store_true",
                    help="axisymmetric r-z run: the grid's y is the radius and "
                         "the axis is y = 0, so the measured chord (axis to the "
                         "phi contour) is the neck RADIUS; report width = 2x.")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.run_dir / "vtkOut" / "solV_*.vts")))
    if not files:
        sys.exit(f"no solV_*.vts under {args.run_dir}/vtkOut")
    tmap = step_times(args.run_dir / "outp.txt")

    rows = []
    centers = None
    for fn in files:
        step = int(Path(fn).stem.split("_")[1])
        try:
            phi, x, y = read_vts_phi(fn)
        except Exception as e:
            print(f"  WARNING: skipping {Path(fn).name} ({e}) — "
                  f"empty/corrupt file, likely an incomplete transfer",
                  file=sys.stderr)
            continue
        w = np.array([chord_width(phi[:, j], y, args.phi) for j in range(len(x))])
        if args.axisym:
            # chord runs from the axis (y = 0, inside the grain) up to the
            # contour: it IS the local radius; physical width is its double.
            w = 2.0 * w

        if centers is None:
            # Grain centers = the two prominent local maxima of w(x). A local
            # max is a column beating its +-5-column neighborhood and taller
            # than 0.5*max(w); the leftmost and rightmost such columns are
            # the two grain centers. (Do NOT use the global minimum to split:
            # the grain TIPS have the smallest nonzero w and would hijack it.)
            thr = 0.5 * w.max()
            peaks = [j for j in range(5, len(w) - 5)
                     if w[j] >= thr and w[j] == w[j - 5:j + 6].max()]
            if len(peaks) < 2:
                sys.exit("could not locate two grain-center peaks in w(x)")
            centers = (peaks[0], peaks[-1])

        lo, hi = centers
        interior = np.arange(lo + 1, hi)
        interior = interior[w[interior] > 0]
        if len(interior) == 0:
            neck, xneck = 0.0, np.nan
        else:
            jn = interior[np.argmin(w[interior])]
            neck, xneck = w[jn], x[jn]
        rows.append((tmap.get(step, np.nan), neck, xneck))

    out = args.out or (args.run_dir / "neck_width.csv")
    with open(out, "w") as f:
        f.write("t_s,neck_width_m,x_neck_m\n")
        for t, nw, xn in rows:
            f.write(f"{t:.6e},{nw:.6e},{xn:.6e}\n")
    print(f"csv -> {out}")

    t = np.array([r[0] for r in rows]) / 60.0        # minutes
    nw = np.array([r[1] for r in rows]) * 1e6        # um
    print(f"neck width: {nw[0]:.2f} um (t=0) -> {nw[-1]:.2f} um "
          f"(t={t[-1]:.0f} min): {100*(nw[-1]/nw[0]-1):+.1f}%")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(t, nw, lw=2, color="#3d74d9")
    ax.set_xlabel("time [min]")
    ax.set_ylabel("neck width [µm]")
    ax.set_title("Neck width vs time (minimum cross-section of φ ≥ 0.5 body)")
    ax.grid(alpha=0.25, lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    png = Path(str(out).replace(".csv", ".png"))
    fig.savefig(png, dpi=150)
    print(f"plot -> {png}")


if __name__ == "__main__":
    main()
