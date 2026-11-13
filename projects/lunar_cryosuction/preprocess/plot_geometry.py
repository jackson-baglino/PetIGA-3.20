#!/usr/bin/env python3
"""plot_geometry.py: a PDF preview of a geometry .opts file, the domain with its
initial ice and the full mesh.

Top panel: the pore (walls from the mesh), the initial ice drawn from the IC
options in the .opts file (-bridge_*, -wedge_band_*, -ice_grain_*), and the
wall faces named. Bottom panel: every element line of the mesh, as vector
lines, so the PDF can be zoomed in to inspect single elements. A text line
gives the domain lengths, the element counts, the element size, eps and the
degrees of freedom.

The ice drawn is the sharp shape (phi = 1/2 contour of the IC), not the
diffuse field. Opts files without a -geom_file are drawn as the rectangle
-Lx x -Ly with -Nx x -Ny elements.

    venv_cryosuction/bin/python preprocess/plot_geometry.py <geometry.opts> [--out file.pdf]

The builders (build_geometry_*.py) call render() after writing their files.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "postprocess"))
sys.path.insert(0, str(ROOT / "preprocess"))


def parse_opts(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line.startswith("-"):
            continue
        parts = line.split(None, 1)
        out[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return out


def floats(o, key):
    return [float(v) for v in o[key].split(",")] if key in o and o[key] else []


def ice_mask(o, X, Y):
    """Sharp initial ice (phi > 1/2) at points X, Y from the IC options."""
    ice = np.zeros_like(X, dtype=bool)
    cxl, rl, cxr, rr, cy = (floats(o, k) for k in ("-bridge_cxL", "-bridge_rL", "-bridge_cxR",
                                                     "-bridge_rR", "-bridge_cy"))
    ccl = [int(v) for v in floats(o, "-bridge_concaveL")] or [0] * len(cxl)
    ccr = [int(v) for v in floats(o, "-bridge_concaveR")] or [0] * len(cxl)
    for k in range(len(cxl)):        # src/initial_conditions.c BridgeField
        dl, dr = np.hypot(X - cxl[k], Y - cy[k]), np.hypot(X - cxr[k], Y - cy[k])
        left = ((dl >= rl[k]) & (X >= cxl[k])) if ccl[k] else (dl <= rl[k])   # concave: outside, beyond center
        right = ((dr >= rr[k]) & (X <= cxr[k])) if ccr[k] else (dr <= rr[k])
        ice |= left & right
    r1, r2 = floats(o, "-wedge_band_r1"), floats(o, "-wedge_band_r2")
    if r1 and "-wedge_apex_x" in o:
        rho = np.hypot(X - float(o["-wedge_apex_x"]), Y - float(o["-wedge_apex_y"]))
        for a, b in zip(r1, r2):
            if b > a:
                ice |= (rho >= a) & (rho <= b)
    gx, gy, gr = floats(o, "-ice_grain_cx"), floats(o, "-ice_grain_cy"), floats(o, "-ice_grain_R")
    ax_, ay_ = floats(o, "-ice_grain_ax") or gr, floats(o, "-ice_grain_ay") or gr
    for k in range(len(gx)):
        ice |= ((X - gx[k]) / ax_[k]) ** 2 + ((Y - gy[k]) / ay_[k]) ** 2 <= 1.0
    fx, fr, fh = floats(o, "-ice_flat_x"), floats(o, "-ice_flat_R"), floats(o, "-ice_flat_height")
    seed = float(o.get("-ice_flat_seed", 0) or 0)
    ly = float(o.get("-Ly", 0) or 0)
    for k in range(len(fx)):         # src/initial_conditions.c IceFlatField (+ -ice_flat_seed)
        xl = fx[k] - fr[k] + (seed * np.cos(2 * np.pi * Y / ly) if seed and ly else 0.0)
        ice |= (X >= xl) & (X <= fx[k] + fr[k]) & (Y <= fh[k])
    return ice


def load_surface(o, opts_path):
    """(surface or None, Lx, Ly, Nx, Ny, p)."""
    lx, ly = float(o["-Lx"]), float(o["-Ly"])
    p = int(o.get("-p", 2))
    if "-geom_file" in o:
        from igakit.io import PetIGA
        dat = ROOT / o["-geom_file"]
        if not dat.exists():
            dat = Path(opts_path).parent / o["-geom_file"]
        surf = PetIGA().read(str(dat))
        nx, ny = len(surf.breaks(0)) - 1, len(surf.breaks(1)) - 1
        return surf, lx, ly, nx, ny, p
    return None, lx, ly, int(o["-Nx"]), int(o["-Ny"]), p


def run_summary(e: dict) -> str:
    """One line of run conditions from an experiment .opts (parsed)."""
    import math
    g = [float(v) for v in e.get("-grad_temp0", "0,0,0").split(",")]
    G = -g[0]
    therm = ("isothermal" if G == 0 else f"G = {abs(G):g} K/m, warm {'top (x0)' if G > 0 else 'bottom (x1)'}")
    if e.get("-flag_BC_rhovfix", "0") in ("0", "false"):
        vap = "vapor closed"
    else:
        faces = e.get("-rhovfix_faces", "x0,x1")
        vals = [e.get("-rhovfix_lo")] + ([e.get("-rhovfix_hi")] if "," in faces else [])
        vap = "vapor fixed on " + faces + " at " + " / ".join(f"{float(v):.6g}" for v in vals if v) + " rho_vs(T0)"
    th = ""
    if "-gamma_is" in e:
        c = (float(e.get("-gamma_as", 0.3)) - float(e["-gamma_is"])) / float(e.get("-gamma_ia", 0.109))
        th = f"; wall angle {math.degrees(math.acos(max(-1, min(1, c)))):.0f}°"
    days = float(e.get("-t_final", 0)) / 86400.0
    return f"T0 {float(e.get('-temp', -20)):g} °C, {therm}; {vap}{th}; {days:g} days"


def render(opts_path, out=None, title=None, exp_path=None) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from figcheck import check_figure

    opts_path = Path(opts_path)
    o = parse_opts(opts_path)
    out = Path(out) if out else ROOT / "preprocess" / (opts_path.stem + ".pdf")
    surf, lx, ly, nx, ny, p = load_surface(o, opts_path)
    s = 1e6
    nu, nv = 2400, 400
    if surf is not None:
        P = surf(np.linspace(0, 1, nu), np.linspace(0, 1, nv))
        X, Y = P[..., 0], P[..., 1]
    else:
        X, Y = np.meshgrid(np.linspace(0, lx, nu), np.linspace(0, ly, nv), indexing="ij")
    ice = ice_mask(o, X, Y)
    rec = None
    if "-phi_file" in o:                        # -ic_type phi_file: draw the field the solver reads
        from build_ice_layout import read_phi
        pf = ROOT / o["-phi_file"]
        rec = read_phi(pf if pf.exists() else o["-phi_file"])

    aspect = ly / lx
    w_in = 7.5
    ph = max(0.9, min(3.2, (w_in - 0.9) * aspect))
    fig, ax = plt.subplots(2, 1, figsize=(w_in, 2 * ph + 1.9), constrained_layout=True)
    a, b = ax
    if rec is not None:
        a.contourf(rec[..., 0] * s, rec[..., 1] * s, rec[..., 2], levels=[0.5, 2.0], colors=["#66b3ff"])
    elif ice.any():
        a.contourf(X * s, Y * s, ice.astype(float), levels=[0.5, 1.5], colors=["#66b3ff"])
    walls_on = set(o.get("-wall_faces", "").split(","))
    for face, (xx, yy) in zip(("y0", "y1", "x0", "x1"),
                              ((X[:, 0], Y[:, 0]), (X[:, -1], Y[:, -1]), (X[0], Y[0]), (X[-1], Y[-1]))):
        if face in walls_on:                    # regolith wall (contact-angle BC): thick brown
            a.plot(xx * s, yy * s, "-", color="#8c6d46", lw=2.2, solid_capstyle="butt", clip_on=False, zorder=5)
        else:
            a.plot(xx * s, yy * s, "k-", lw=0.6)
    head = title or opts_path.name
    if exp_path is not None:
        import textwrap
        line = Path(exp_path).stem + ": " + run_summary(parse_opts(Path(exp_path))) + "; brown: regolith wall"
        head += "\n" + "\n".join(textwrap.wrap(line, 120))
    a.set_title(head, fontsize=7.5)

    segs = []
    if surf is not None:
        uu = np.linspace(0, 1, 4000)
        vv = np.linspace(0, 1, 200)
        for ub in surf.breaks(0):
            segs.append(surf([ub], vv)[0][:, :2] * s)
        for vb in surf.breaks(1):
            segs.append(surf(uu, [vb])[:, 0][:, :2] * s)
    else:
        for xb in np.linspace(0, lx, nx + 1):
            segs.append([(xb * s, 0), (xb * s, ly * s)])
        for yb in np.linspace(0, ly, ny + 1):
            segs.append([(0, yb * s), (lx * s, yb * s)])
    b.add_collection(LineCollection(segs, colors="#1f6fd0", linewidths=0.04))
    for x in (a, b):
        x.set_xlim(X.min() * s, X.max() * s)
        x.set_ylim(Y.min() * s, Y.max() * s)
        x.set_aspect("equal")
        x.set_xlabel("x [µm]")
        x.set_ylabel("y [µm]")

    eps = float(o["-eps"]) if "-eps" in o else float("nan")
    hx, hy = lx / nx, ly / ny
    gap = (Y.max(axis=1) - Y.min(axis=1)) if surf is not None else np.array([ly])
    wtxt = (f"Ly = {ly*s:.1f} µm (pore width {gap.min()*s:.1f}-{gap.max()*s:.1f} µm)"
            if gap.max() - gap.min() > 1e-9 else f"Ly = {ly*s:.1f} µm")
    dof = (nx + p) * (ny + p) * 3
    walls = o.get("-wall_faces", "none")
    b.set_title(f"mesh, every element line (zoom in): Lx = {lx*s:.1f} µm, {wtxt};\n"
                f"{nx} x {ny} elements, h = {hx*s:.3f} x {hy*s:.3f} µm (at full width), "
                f"eps = {eps*s:.4f} µm, {dof/1e6:.2f} M DOF; wall faces: {walls}", fontsize=7.5)
    bad = check_figure(fig)
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}  figcheck: {bad or 'OK'}")
    return bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("opts", nargs="+", help="geometry .opts file(s)")
    ap.add_argument("--out", default=None, help="output PDF (one opts file only)")
    ap.add_argument("--exp", default=None, help="experiment .opts: adds the run conditions to the title")
    a = ap.parse_args(argv)
    if a.out and len(a.opts) > 1:
        ap.error("--out needs a single opts file")
    bad = [render(f, a.out, exp_path=a.exp) for f in a.opts]
    return 1 if any(bad) else 0


if __name__ == "__main__":
    raise SystemExit(main())
