"""
grain_shrinkage.py — per-grain volume and equivalent radius vs time, for the
two-grain sintering runs.

WHY THIS EXISTS
---------------
Nothing else measures grain SIZE over time. plot_mass.py integrates the whole
domain, pplib.grain_radii() reads the INITIAL radii out of the opts, and
postprocess/scratch/track_grain_mass.py splits at the domain midpoint and
reports mass, not radius. To fit the chamber humidity against Molaro et al.
(2019) we need each grain's equivalent radius per snapshot.

HOW
---
Split the domain at the neck plane and integrate the ice phase on each side.

  axisymmetric (-axisym 1, the Molaro geometry): x is the symmetry axis, y the
      radius, so the measure is dV = 2*pi*y dy dx and

          V = int int  phi_i * 2*pi*y  dy dx

      is a true 3-D volume in m^3. R_eq = (3V/(4*pi))^(1/3).

  planar 2-D: there is no third dimension, so the integral is an AREA per unit
      depth [m^2]. R_eq is then the area-equivalent disc radius sqrt(A/pi), and
      it is NOT comparable to a sphere radius. The CSV labels which one it is.

THE SPLIT PLANE comes from neck_width.py's `x_neck_m` column when a
neck_width.csv is present (it is re-located per snapshot there, so it tracks a
drifting neck), otherwise from the RADICAL PLANE of the two spheres,

    x = c0 + (d^2 + R0^2 - R1^2) / (2 d),      d = |c1 - c0|

which is where the two sphere surfaces meet. Do NOT use the midpoint between
centres: for the Molaro pair (72.5 / 101 um) that sits 14 um from the actual
contact and biases the small grain by +1.3 %, which is half the shrinkage
signal we are trying to fit. The verification gate catches exactly this.

Once the plane is within a few um of the contact the answer is insensitive to
it -- the neck disc is a ~1e-5 fraction of a grain's volume -- so tracking a
drifting neck matters much less than not being systematically off.

WHICH GRAIN IS WHICH is read from -ice_grain_cx / -ice_grain_R, not assumed:
whichever centre lies on the low-x side of the split is the "low" grain.

Usage
-----
  python postprocess/grain_shrinkage.py <run_dir>
  python postprocess/grain_shrinkage.py <run_dir> --data inputs/validation/molaro2019_fig11_T-20.csv
  python postprocess/grain_shrinkage.py <run_dir> --out /tmp/gs.csv --no-plot
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

# numpy 2.0 renamed trapz -> trapezoid; the venv here is 2.x but HPC may be 1.x.
_trapz = getattr(np, "trapezoid", None) or np.trapz

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib  # noqa: E402

# Colour by ENTITY (which grain), not by model-vs-data; validated pair.
C_LARGE, C_SMALL, C_INK = "#3d74d9", "#d1495b", "#3f434a"


def grain_centres(opts: dict):
    """(cx list [m], R list [m]) from the staged opts, or (None, None)."""
    try:
        cx = [float(v) for v in opts["-ice_grain_cx"].split(",")]
        rr = [float(v) for v in opts["-ice_grain_R"].split(",")]
    except (KeyError, ValueError):
        return None, None
    return cx, rr


def radical_plane(cx, rr) -> float | None:
    """x where the two sphere surfaces meet, along the line of centres.

    For tangent spheres (d = R0 + R1) this reduces to the contact point
    c0 + R0. The midpoint between CENTRES is a different plane entirely
    whenever the radii differ -- 14 um away for the Molaro pair.
    """
    if not cx or len(cx) < 2:
        return None
    lo, hi = (0, 1) if cx[0] <= cx[1] else (1, 0)
    d = abs(cx[hi] - cx[lo])
    if d <= 0.0:
        return None
    return cx[lo] + (d * d + rr[lo] ** 2 - rr[hi] ** 2) / (2.0 * d)


def integrate_sides(phi, X, Y, x_split: float, axisym: bool):
    """(V_low, V_high) of the ice phase either side of x = x_split.

    Trapezoid on the structured grid. Axisymmetric integrands carry the 2*pi*y
    weight, so the result is m^3; planar ones are m^2 (area per unit depth).

    The split inserts an interpolated node AT x_split and gives it to both
    sides, rather than cutting the arrays at an index. Cutting at an index
    silently drops the trapezoid segment that straddles the split, which is
    exactly where the neck material lives -- it cost 1.3 % on the small grain
    in studies/molaro_2019/verification/verify_grain_shrinkage.py. With the
    inserted node, V_low + V_high reproduces the whole-domain integral.
    """
    x = X[0, :]
    y = Y[:, 0]
    w = phi * (2.0 * math.pi * Y) if axisym else phi
    # Integrate over y first -> a line density in x, then split and integrate.
    lin = _trapz(w, y, axis=0)                       # shape (nx,)
    xs = float(min(max(x_split, x[0]), x[-1]))
    k = int(np.searchsorted(x, xs))
    ls = float(np.interp(xs, x, lin))
    lo = _trapz(np.append(lin[:k], ls), np.append(x[:k], xs))
    hi = _trapz(np.insert(lin[k:], 0, ls), np.insert(x[k:], 0, xs))
    return float(lo), float(hi)


def equiv_radius(V: float, axisym: bool) -> float:
    """Sphere-equivalent radius (axisym, m^3) or disc-equivalent (planar, m^2)."""
    if V <= 0.0:
        return 0.0
    return (3.0 * V / (4.0 * math.pi)) ** (1.0 / 3.0) if axisym else math.sqrt(V / math.pi)


def neck_planes(run_dir: Path) -> dict:
    """{step: x_neck_m} from a neighbouring neck_width.csv, or {}."""
    csv_path = run_dir / "neck_width.csv"
    if not csv_path.is_file():
        return {}
    out = {}
    with csv_path.open() as fh:
        rd = csv.DictReader(fh)
        for i, row in enumerate(rd):
            try:
                out[i] = float(row["x_neck_m"])
            except (KeyError, ValueError, TypeError):
                pass
    return out


def read_validation(path: Path):
    """[(t_s, R_large_m, R_small_m), ...] — diameters halved to radii."""
    rows = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = next(csv.reader([line]))
            rows.append((float(f[0]) * 60.0, float(f[4]) * 1e-6 / 2.0,
                         float(f[5]) * 1e-6 / 2.0))
    return rows


def make_figure(rows, data, axisym, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t_min = [r["t_s"] / 60.0 for r in rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))

    for key, col, lab in (("R_large_m", C_LARGE, "large grain"),
                          ("R_small_m", C_SMALL, "small grain")):
        R = [r[key] * 1e6 for r in rows]
        pct = 100.0 * (R[-1] / R[0] - 1.0) if R[0] else float("nan")
        ax.plot(t_min, R, "-", color=col, lw=2, zorder=3,
                label=f"model, {lab}  ({pct:+.2f} %)")

    if data:
        td = [d[0] / 60.0 for d in data]
        for idx, col, lab in ((1, C_LARGE, "large"), (2, C_SMALL, "small")):
            Rd = [d[idx] * 1e6 for d in data]
            pct = 100.0 * (Rd[-1] / Rd[0] - 1.0)
            ax.plot(td, Rd, "o", mfc="none", mec=col, mew=1.8, ms=8, zorder=4,
                    label=f"Molaro 2019, {lab}  ({pct:+.2f} %)")

    ax.set_xlabel("time  [min]")
    ax.set_ylabel(("sphere-equivalent radius" if axisym
                   else "disc-equivalent radius") + "  [$\\mu$m]")
    ax.set_title("Grain size vs time", fontsize=11, color=C_INK)
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=C_INK, labelsize=9)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    print(f"  wrote {out_png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None, help="CSV path")
    ap.add_argument("--data", type=Path, default=None,
                    help="validation CSV to overlay (default: the Molaro -20 C table)")
    ap.add_argument("--phi-field", default="IcePhase")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    run = args.run_dir
    snaps = sorted((run / "vtkOut").glob("solV_*.vts"), key=pplib.step_of)
    if not snaps:
        sys.exit(f"no vtkOut/solV_*.vts under {run} — run plot_fields.py first")

    opts = pplib.read_opts(str(run))
    axisym = str(opts.get("-axisym", "0")).strip() in ("1", "true", "True")
    cx, rr = grain_centres(opts)
    x_mid = radical_plane(cx, rr)
    planes = neck_planes(run)
    tmap = pplib.step_times(run)

    if not axisym:
        print("  NOTE: -axisym is not set; reporting AREA per unit depth and a"
              " disc-equivalent radius. Not comparable to a sphere radius.")
    if cx and len(cx) >= 2:
        lo_is_large = rr[0] > rr[1] if cx[0] < cx[1] else rr[1] > rr[0]
        print(f"  grains: cx = {cx} m, R = {rr} m"
              f"  -> the {'low' if lo_is_large else 'high'}-x grain is the large one")
        if x_mid is not None:
            print(f"  default split plane (radical) = {x_mid:.6e} m"
                  f"   [midpoint between centres would be {0.5*(cx[0]+cx[1]):.6e}]")
    else:
        lo_is_large = False
        print("  WARNING: -ice_grain_cx/-ice_grain_R not found in the staged opts;"
              " assuming the high-x grain is the large one.")

    rows = []
    for fn in snaps:
        step = pplib.step_of(fn)
        fields, X, Y = pplib.read_vts(fn, want=[args.phi_field])
        phi = fields[args.phi_field]
        x_split = planes.get(len(rows), x_mid)
        if x_split is None:
            x_split = 0.5 * (X[0, 0] + X[0, -1])
        v_lo, v_hi = integrate_sides(phi, X, Y, x_split, axisym)
        V_large, V_small = (v_lo, v_hi) if lo_is_large else (v_hi, v_lo)
        rows.append(dict(step=step, t_s=tmap.get(step, float("nan")),
                         x_neck_m=x_split,
                         V_large_m3=V_large, V_small_m3=V_small,
                         R_large_m=equiv_radius(V_large, axisym),
                         R_small_m=equiv_radius(V_small, axisym)))

    out = args.out or (run / "grain_shrinkage.csv")
    unit = "m3" if axisym else "m2_per_m"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["t_s", "R_large_m", "R_small_m",
                    f"V_large_{unit}", f"V_small_{unit}", "x_neck_m", "step"])
        for r in rows:
            w.writerow([f"{r['t_s']:.6e}", f"{r['R_large_m']:.6e}", f"{r['R_small_m']:.6e}",
                        f"{r['V_large_m3']:.6e}", f"{r['V_small_m3']:.6e}",
                        f"{r['x_neck_m']:.6e}", r["step"]])
    print(f"  wrote {out}")

    if len(rows) >= 2:
        for key, lab in (("R_large_m", "large"), ("R_small_m", "small")):
            a, b = rows[0][key], rows[-1][key]
            pct = 100.0 * (b / a - 1.0) if a else float("nan")
            print(f"  {lab:5s} grain: {a*1e6:8.3f} -> {b*1e6:8.3f} um   ({pct:+.2f} %)")
        print("  Fit the humidity to the LARGE grain: in Molaro's table only its trend")
        print("  clears the measurement noise (S/N 4.4 vs 1.1).")

    if not args.no_plot:
        data_path = args.data
        if data_path is None:
            cand = Path(__file__).resolve().parent.parent / \
                "inputs/validation/molaro2019_fig11_T-20.csv"
            data_path = cand if cand.is_file() else None
        data = read_validation(data_path) if data_path else None
        make_figure(rows, data, axisym, str(out).replace(".csv", ".png"))


if __name__ == "__main__":
    main()
