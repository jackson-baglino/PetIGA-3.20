"""Start frame of one type-6 (dog-bone) run, with its boundary conditions.

Everything drawn comes from the run folder, not from a sketch:
  - the ice field phi and the pore walls: the first snapshot (sol_00000.dat on
    igasol.dat, evaluated with pplib.iter_snapshots);
  - element sizes: the knot lines of the mesh in igasol.dat;
  - Nx, Ny and the boundary-condition values: the run log (outp.txt, the
    "MESH & DISCRETIZATION" and "BOUNDARY CONDITIONS" blocks);
  - eps, G, T0, contact angle, faces: the staged .opts files (pplib.read_opts),
    with the option name printed next to each label;
  - neck width and extent: the "# DOGBONE:" line of the geometry .opts;
  - lens extent: phi = 1/2 crossings of the t = 0 field, on the axis and on
    the first row of nodes off each wall.
The lower panels give T
and rho_v along the pore axis at t = 0.

Checked with postprocess/figcheck.py at the placed width (7 in); exits 1 when
a check fails.

    venv_cryosuction/bin/python docs/run_plan/t6_startframe.py <run folder> <out.png> [--sat-day 1]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "postprocess"))
import pplib  # noqa: E402
from figcheck import check_figure  # noqa: E402

UM = 1e6
FIG_W = 7.0
INK, GRAY, ROCK_F, ROCK_E = "#222222", "#555555", "#d9cbb3", "#8c7350"


def geom_opts(run: Path) -> Path:
    for f in run.glob("*.opts"):
        if "-geom_file" in f.read_text():
            return f
    raise FileNotFoundError(f"{run}: no geometry .opts")


def dogbone(run: Path) -> dict:
    m = re.search(r"#\s*DOGBONE:(.*)", geom_opts(run).read_text())
    if not m:
        raise ValueError("geometry .opts has no '# DOGBONE:' line")
    tok = m.group(1).split()
    return {tok[i]: float(tok[i + 1]) for i in range(0, len(tok), 2)}


def log_values(run: Path) -> dict:
    txt = (run / "outp.txt").read_text(errors="replace")
    out = {}
    m = re.search(r"Nx = (\d+),\s+Ny = (\d+)", txt)
    out["nx"], out["ny"] = int(m.group(1)), int(m.group(2))
    rows = {}
    for face in ("x = 0", "x = Lx", "y = 0", "y = Ly"):
        m = re.search(r"^\s*" + re.escape(face) + r"\s+\(\w+\)\s+(.*)$", txt, re.M)
        rows[face] = m.group(1) if m else ""
    out["rows"] = rows
    for face, key in (("x = 0", "x0"), ("x = Lx", "x1")):
        r = rows[face]
        mt = re.search(r"Dirichlet\s+(-?[\d.]+) C", r)
        mv = re.search(r"Dirichlet\s+([\d.eE+-]+) kg/m3", r)
        out[f"T_{key}"] = float(mt.group(1)) if mt else None
        out[f"rv_{key}"] = float(mv.group(1)) if mv else None
    m = re.search(r"theta = ([\d.]+)°", txt)
    out["theta_log"] = float(m.group(1)) if m else None
    return out


def mesh_sizes(run: Path):
    from igakit.io import PetIGA
    nrb = PetIGA().read(str(run / "igasol.dat"))
    B = nrb(nrb.breaks(0), nrb.breaks(1))
    hx = np.hypot(*(np.diff(B[..., k], axis=0) for k in (0, 1)))
    hy = np.hypot(*(np.diff(B[..., k], axis=1) for k in (0, 1)))
    return hx, hy, nrb.degree


def crossings_x(x, f):
    """x positions where f crosses 1/2 (linear interpolation)."""
    s = np.sign(f - 0.5)
    i = np.where(s[:-1] * s[1:] < 0)[0]
    return x[i] + (0.5 - f[i]) * (x[i + 1] - x[i]) / (f[i + 1] - f[i])


def rng(h, nd=2):
    lo, hi = h.min() * UM, h.max() * UM
    return f"{lo:.{nd}f}" if abs(hi - lo) < 10 ** -nd else f"{lo:.{nd}f}–{hi:.{nd}f}"


def orient(s):
    """X, Y, phi, T, rho_v of a snapshot with the along-wall index second, in um."""
    X, Y, F = s.X * UM, s.Y * UM, (s.phi, s.T, s.rhov)
    if np.ptp(Y[0, :]) > np.ptp(Y[:, 0]):
        X, Y, F = X.T, Y.T, tuple(f.T for f in F)
    return (X, Y, *F)


def main(run: Path, out: Path, sat_day: float = 1.0) -> int:
    o = pplib.read_opts(str(run))
    db, lg = dogbone(run), log_values(run)
    hx, hy, deg = mesh_sizes(run)
    it = pplib.iter_snapshots(str(run))
    X, Y, phi, T, rv = orient(next(it))
    ss = next((q for q in it if q.t >= sat_day * 86400), None)
    if ss is None:
        raise ValueError(f"no snapshot at or after {sat_day} d")
    _, _, phi_s, T_s, rv_s = orient(ss)
    jcs = phi_s.shape[0] // 2
    h0 = rv[phi.shape[0] // 2] / pplib.rho_vs(T[phi.shape[0] // 2])
    hs = rv_s[jcs] / pplib.rho_vs(T_s[jcs])
    lx, W = X.max(), db["W"] * UM
    wc, xn0, xn1, sh = db["wc"] * UM, db["x_n0"] * UM, db["x_n1"] * UM, db["shoulder"] * UM
    eps = float(o["-eps"]) * UM
    G = -float(o["-grad_temp0"].split(",")[0])
    T0 = float(o["-temp"])
    gia, gas, gis = (float(o[k]) for k in ("-gamma_ia", "-gamma_as", "-gamma_is"))
    theta = np.degrees(np.arccos((gas - gis) / gia))
    walls = o.get("-wall_faces", "")
    jc = phi.shape[0] // 2
    x_ax = crossings_x(X[jc], phi[jc])
    x_w0 = crossings_x(X[1], phi[1])
    x_w1 = crossings_x(X[-2], phi[-2])

    plt.rcParams.update({"font.family": "Arial", "font.size": 7.5})
    ex = 1.6                                     # vertical stretch of the pore panel
    ax_w = 6.25
    yl = (-95, W * ex + 34)
    xl = (-20, lx + 20)
    h_pore = ax_w * (yl[1] - yl[0]) / (xl[1] - xl[0])
    top_in, mid_in, prof_in, bot_in = 0.62, 0.95, 1.05, 0.42
    H = top_in + h_pore + mid_in + prof_in + bot_in
    fig = plt.figure(figsize=(FIG_W, H))
    ax = fig.add_axes((0.2 / FIG_W, (H - top_in - h_pore) / H, ax_w / FIG_W, h_pore / H))

    # ice field and phi = 1/2
    pc = ax.pcolormesh(X, Y * ex, phi, cmap="Blues", vmin=0, vmax=1, shading="gouraud", rasterized=True,
                       zorder=1)
    ax.contour(X, Y * ex, phi, levels=[0.5], colors=["#1f4e8c"], linewidths=0.8, zorder=3)
    def draw_walls(a):
        """Hatched regolith from the mesh boundary rows, dashed pore ends."""
        for j, sgn in ((0, -1), (-1, 1)):
            x, y = X[j], Y[j] * ex
            yb = (y.min() - 9) if sgn < 0 else (y.max() + 9)
            a.fill(np.r_[x, x[::-1]], np.r_[y, np.full_like(y, yb)], fc=ROCK_F, ec=ROCK_E, lw=0.8,
                   hatch="////", zorder=2)
        for xe, k in ((0, 0), (lx, -1)):
            a.plot([xe, xe], [Y[0, k] * ex, Y[-1, k] * ex], color=INK, lw=1.0, ls=(0, (3, 2)), zorder=3)

    draw_walls(ax)

    def hdim(a, b, yy, txt, va="top"):
        ax.annotate("", (a, yy), (b, yy), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.6,
                                                          shrinkA=0, shrinkB=0))
        ax.text((a + b) / 2, yy + (-3 if va == "top" else 3), txt, ha="center", va=va, color=GRAY)

    def vdim(x, a, b, txt, side="left"):
        ax.annotate("", (x, a * ex), (x, b * ex), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.6,
                                                                  shrinkA=0, shrinkB=0), zorder=6)
        dx = -6 if side == "left" else 6
        ax.text(x + dx, (a + b) / 2 * ex, txt, ha="right" if side == "left" else "left", va="center",
                color=GRAY, zorder=6, bbox=dict(fc="white", ec="none", pad=0.6), gid="overlap_ok")

    hdim(0, lx, -78, rf"$L_x$ = {lx:.0f} µm")
    hdim(xn0, xn1, -36, f"neck x = {xn0:.0f}–{xn1:.0f} µm")
    hdim(xn0 - sh, xn0, -36, f"{sh:.0f}")
    hdim(xn1, xn1 + sh, -36, f"{sh:.0f}")
    vdim(160, 0, W, rf"$W$ = {W:.0f} µm")
    yb0, yb1 = (W - wc) / 2, (W + wc) / 2
    vdim((xn0 + xn1) / 2, yb0, yb1, rf"$w_c$ = {wc:.0f} µm", side="right")
    if x_ax.size >= 2:
        a, b = x_ax.min(), x_ax.max()
        ax.annotate("", (a, W * ex + 24), (b, W * ex + 24),
                    arrowprops=dict(arrowstyle="<->", color="#1f4e8c", lw=0.6, shrinkA=0, shrinkB=0))
        wl = [v for v in (x_w0, x_w1) if v.size >= 2]
        wtxt = (f"; at wall {min(v.min() for v in wl):.0f}–{max(v.max() for v in wl):.0f}" if wl else "")
        ax.text(a - 6, W * ex + 24, f"lens (t = 0, $\\phi$ = 1/2): axis {a:.0f}–{b:.0f} µm{wtxt}",
                ha="right", va="center", color="#1f4e8c")
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")

    # colorbar
    cax = fig.add_axes(((0.2 + ax_w + 0.05) / FIG_W, (H - top_in - h_pore + 0.32) / H, 0.06 / FIG_W,
                        (h_pore - 0.45) / H))
    cb = fig.colorbar(pc, cax=cax, ticks=[0, 0.5, 1])
    cb.ax.tick_params(labelsize=7, length=2, pad=1)
    cax.set_title(r"$\phi$", fontsize=7.5, pad=2)

    # boundary conditions: pore ends (top block), walls (under the pore)
    rvs0 = float(pplib.rho_vs(T0))
    left = (r"$x = 0$ (warm end; dashed)" + "\n"
            + rf"$T$ = {lg['T_x0']:.3f} °C, fixed   [-Tfix_faces x0]" + "\n"
            + rf"$\rho_v$ = {lg['rv_x0']:.4e} kg m$^{{-3}}$ = $\rho_{{vs}}(T)$, fixed   "
            + f"[-rhovfix_lo {float(o['-rhovfix_lo']):.5f}]" + "\n"
            + r"$\phi$: $\partial_n\phi$ = 0   [-phi_open_bc 0]")
    right = (r"$x = L_x$ (cold end; dashed)" + "\n"
             + rf"$T$ = {lg['T_x1']:.3f} °C, fixed   [-Tfix_faces x1]" + "\n"
             + rf"$\rho_v$ = {lg['rv_x1']:.4e} kg m$^{{-3}}$ = $\rho_{{vs}}(T)$, fixed   "
             + f"[-rhovfix_hi {float(o['-rhovfix_hi']):.5f}]" + "\n"
             + r"$\phi$: $\partial_n\phi$ = 0   [-phi_open_bc 0]")
    fig.text(0.2 / FIG_W, 1 - 0.05 / H, left, ha="left", va="top", linespacing=1.3)
    fig.text(1 - 0.15 / FIG_W, 1 - 0.05 / H, right, ha="right", va="top", linespacing=1.3)
    if abs(theta - 90) < 1e-6:
        phiw = r"$\phi$: $\partial_n\phi$ = 0, $\theta$ = 90°"
    else:
        phiw = rf"$\phi$: wall energy, $\theta$ = {theta:.0f}°"
    ywall = H - top_in - h_pore
    wall_txt = (f"walls y = 0 and y = W(x) (hatched, from the mesh) [-wall_faces {walls}]: {phiw} "
                f"[-gamma_is {gis:g}, -gamma_as {gas:g}," + "\n" + f"-gamma_ia {gia:g}]; "
                r"zero flux of $T$ and $\rho_v$ (natural condition)." + "\n"
                + rf"Imposed gradient $G$ = {G:.0f} K/m, warm at $x$ = 0 [-grad_temp0]; "
                + rf"$T_0$ = {T0:.0f} °C [-temp], $\rho_{{vs}}(T_0)$ = {rvs0:.4e} kg m$^{{-3}}$.  "
                + f"Vertical scale ×{ex:g}." + "\n"
                + f"Mesh: {lg['nx']} × {lg['ny']} elements, degree {deg[0]} [igasol.dat, outp.txt]; "
                + rf"$h_x$ = {rng(hx)} µm, $h_y$ = {rng(hy)} µm; "
                + rf"$\varepsilon$ = {eps:.4f} µm [-eps]; "
                + f"{3 * (lg['nx'] + deg[0]) * (lg['ny'] + deg[1]) / 1e6:.2f} M DOF.")
    fig.text(0.2 / FIG_W, (ywall - 0.04) / H, wall_txt, ha="left", va="top", linespacing=1.3)

    # axis profiles: T at t = 0; humidity rho_v/rho_vs(T) at t = 0 and at sat_day
    pw = (FIG_W - 0.2 - 0.55 - 0.6) / 2
    yb = bot_in / H
    a1 = fig.add_axes((0.55 / FIG_W, yb, pw / FIG_W, prof_in / H))
    a1.plot(X[jc], T[jc], color=INK, lw=0.9)
    a1.plot([0, lx], (lg["T_x0"], lg["T_x1"]), "o", ms=3, color="#b2182b", clip_on=False)
    a1.set_title(r"$T$ on the axis at $t$ = 0 [°C]; dots: end values in the log", fontsize=7.5, loc="left", pad=3)
    a2 = fig.add_axes(((0.55 + pw + 0.6) / FIG_W, yb, pw / FIG_W, prof_in / H))
    a2.plot(X[jc], h0, color=GRAY, lw=0.9, ls=(0, (3, 2)), label=r"$t$ = 0")
    a2.plot(X[jcs], hs, color=INK, lw=0.9, label=rf"$t$ = {ss.t / 86400:.2f} d")
    lens = phi_s[jcs] > 0.5
    if lens.any():
        a2.axvspan(X[jcs][lens].min(), X[jcs][lens].max(), color="#c6dbef", zorder=0)
    a2.legend(fontsize=7, frameon=False, loc="lower left", handlelength=2)
    a2.set_title(r"$\rho_v/\rho_{vs}(T)$ on the axis, $T$ local; blue: ice", fontsize=7.5, loc="left", pad=3)
    for a in (a1, a2):
        a.axvspan(xn0, xn1, color="#eeeeee", zorder=0)
        a.set_xlim(0, lx)
        a.set_xlabel("x [µm]")
        a.tick_params(labelsize=7, length=2)
        a.ticklabel_format(axis="y", useOffset=False)
    bad = check_figure(fig, min_pt=7.0, placed_width_in=FIG_W)
    for b in bad:
        print("  FAIL:", b)
    fig.savefig(out, dpi=300)
    print(f"{out}: {len(bad)} problems, {FIG_W} x {H:.2f} in")
    print(f"  axis humidity: t = 0 {h0.min():.6f} to {h0.max():.6f}; t = {ss.t / 86400:.3f} d "
          f"{hs.min():.6f} to {hs.max():.6f}")
    print(f"  lens phi=1/2 on axis: {np.round(x_ax, 1)} um; near walls: {np.round(x_w0, 1)}, "
          f"{np.round(x_w1, 1)} um; theta from gammas {theta:.2f} deg (log {lg['theta_log']})")
    return 1 if bad else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--sat-day", type=float, default=1.0, help="day of the second humidity profile")
    a = ap.parse_args()
    sys.exit(main(a.run, a.out, a.sat_day))
