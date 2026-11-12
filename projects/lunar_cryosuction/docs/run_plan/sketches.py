"""Schematics for the run plan (docs/run_plan/build_runplan.py), drawn from the
real geometry and initial-condition files.

Each panel shows the pore as the solver sees it: x along the pore (warm end,
when there is a gradient, at x = 0 on the left), the walls from the mesh, the
initial ice (phi = 1/2 of the IC: the -bridge / -ice_grain options, or the
phi_file), dimensions, and a mesh line computed from the mesh itself
(elements, element-size range from the knot lines, eps, degrees of freedom).
The boundary conditions are printed above and below the panel. A panel whose
height would be under about 0.8 in is stretched vertically; the stretch is
printed on the panel ("vertical scale x2").

Every figure is checked with postprocess/figcheck.py (no overlapping text,
printed fonts >= 7 pt at the width placed in the document) before it is saved;
the script exits 1 if any check fails.

    venv_cryosuction/bin/python docs/run_plan/sketches.py <fig_dir> <scratch_geom_dir>

<scratch_geom_dir> holds the proposed geometries that are not in inputs/ yet
(built with preprocess/build_geometry_dogbone.py; commands in the run plan,
Appendix A).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "preprocess"))
sys.path.insert(0, str(ROOT / "postprocess"))
from build_ice_layout import read_phi  # noqa: E402
from figcheck import check_figure  # noqa: E402
from plot_geometry import ice_mask, load_surface, parse_opts  # noqa: E402

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "fig")
SCR = Path(sys.argv[2]) if len(sys.argv) > 2 else None
OUT.mkdir(parents=True, exist_ok=True)
GEO = ROOT / "inputs/geometry/cryosuction"
FIG_W = 6.5          # placed width in the document [in]
FAIL: list[str] = []

ICE_F, ICE_E = "#bfe0f5", "#1f6fb2"
ROCK_F, ROCK_E = "#e3d5bb", "#8c7350"
INK, GRAY = "#222222", "#555555"
FS = 7.5
plt.rcParams.update({"font.family": "Arial", "font.size": FS, "mathtext.fontset": "cm",
                     "hatch.color": ROCK_E, "hatch.linewidth": 0.5})
UM = 1e6


# ------------------------------------------------------------------ data
def load(opts: Path):
    """Pore walls, ice and mesh numbers of one geometry .opts."""
    o = parse_opts(opts)
    surf, lx, ly, nx, ny, p = load_surface(o, opts)
    if surf is not None:
        u, v = np.linspace(0, 1, 900), np.linspace(0, 1, 220)
        P = surf(u, v)
        X, Y = P[..., 0], P[..., 1]
        B = surf(surf.breaks(0), surf.breaks(1))
        hx = np.hypot(*(np.diff(B[..., k], axis=0) for k in (0, 1)))
        hy = np.hypot(*(np.diff(B[..., k], axis=1) for k in (0, 1)))
        bot = np.c_[X[:, 0], Y[:, 0]]
        top = np.c_[X[:, -1], Y[:, -1]]
    else:
        X, Y = np.meshgrid(np.linspace(0, lx, 900), np.linspace(0, ly, 450), indexing="ij")
        hx, hy = np.array([lx / nx]), np.array([ly / ny])
        bot = np.c_[[0, lx], [0, 0]]
        top = np.c_[[0, lx], [ly, ly]]
    if o.get("-ic_type") == "phi_file":
        f = ROOT / o["-phi_file"]
        a = read_phi(f)
        ice = (a[..., 0], a[..., 1], a[..., 2])
    else:
        ice = (X, Y, ice_mask(o, X, Y).astype(float))
    walls = o.get("-wall_faces", "").split(",")
    return dict(o=o, lx=lx, ly=ly, nx=nx, ny=ny, hx=hx, hy=hy, bot=bot, top=top, ice=ice, walls=walls,
                eps=float(o["-eps"]), dof=3 * (nx + p) * (ny + p))


def rng(h):
    lo, hi = h.min() * UM, h.max() * UM
    return f"{lo:.2f}" if abs(hi - lo) < 0.01 else f"{lo:.2f}–{hi:.2f}"


def mesh_line(g, br=False):
    sep = "\n" if br else ";  "
    return (f"mesh: {g['nx']} × {g['ny']} quadratic $C^1$ elements{sep}$h_x$ = {rng(g['hx'])} µm,  "
            f"$h_y$ = {rng(g['hy'])} µm;  $\\varepsilon$ = {g['eps'] * UM:.3f} µm;  "
            f"{g['dof'] / 1e6:.2f} M DOF")


# ------------------------------------------------------------------ drawing
def extent(g, ex):
    """Drawn wall thickness and the panel's x and y limits [um, y stretched by ex]."""
    lx, ly = g["lx"] * UM, g["ly"] * UM
    t = 0.10 * max(ly, 60) * ex
    ymin = min(g["bot"][:, 1].min() * UM * ex - t, 0)
    ymax = g["top"][:, 1].max() * UM * ex + t
    return t, (-0.02 * lx, 1.05 * lx), (ymin - 0.32 * (ymax - ymin) - 12, ymax + 0.20 * (ymax - ymin))


def panel(ax, g, ex, dims=(), title=None, note=None):
    lx = g["lx"] * UM
    t, xl, yl = extent(g, ex)
    for curve, sgn, face in ((g["bot"], -1, "y0"), (g["top"], 1, "y1")):
        x, y = curve[:, 0] * UM, curve[:, 1] * UM * ex
        if face in g["walls"]:
            yb = (y.min() - t) if sgn < 0 else (y.max() + t)
            ax.fill(np.r_[x, x[::-1]], np.r_[y, np.full_like(y, yb)], fc=ROCK_F, ec=ROCK_E, lw=0.9,
                    hatch="////", zorder=2)
        else:
            ax.plot(x, y, color=GRAY, lw=0.8, zorder=2)
    for xe, k in ((0, 0), (lx, -1)):
        ya, yb = g["bot"][k, 1] * UM * ex, g["top"][k, 1] * UM * ex
        if "x1" in g["walls"] and k == -1:
            ax.fill([xe, xe + t, xe + t, xe], [ya - t, ya - t, yb + t, yb + t], fc=ROCK_F, ec=ROCK_E,
                    lw=0.9, hatch="////", zorder=2)
        else:
            ax.plot([xe, xe], [ya, yb], color=INK, lw=1.0, ls=(0, (3, 2)), zorder=3)
    X, Y, F = g["ice"]
    ax.contourf(X * UM, Y * UM * ex, F, levels=[0.5, 2], colors=[ICE_F], zorder=3)
    ax.contour(X * UM, Y * UM * ex, F, levels=[0.5], colors=[ICE_E], linewidths=0.9, zorder=4)
    ymin = min(g["bot"][:, 1].min() * UM * ex - t, 0)
    ymax = g["top"][:, 1].max() * UM * ex + t
    for kind, a, b, pos, txt in dims:
        if kind == "h":            # horizontal dimension under the pore
            yy = ymin - 0.10 * (ymax - ymin) - 4
            ax.annotate("", (a, yy), (b, yy), arrowprops=dict(arrowstyle="<->", color=GRAY, lw=0.6))
            ax.text((a + b) / 2, yy - 0.03 * (ymax - ymin) - 2, txt, ha="center", va="top", color=GRAY)
        else:                      # vertical dimension at x = pos, from y a to b [um]
            ax.annotate("", (pos, a * ex), (pos, b * ex), arrowprops=dict(arrowstyle="<->", color=GRAY,
                                                                           lw=0.6), zorder=6)
            ax.text(pos - 0.02 * lx, (a + b) / 2 * ex, txt, ha="right", va="center", color=GRAY, zorder=6)
    if title:
        ax.text(0, ymax + 0.04 * (ymax - ymin), title, ha="left", va="bottom", fontweight="bold")
    if note:
        ax.text(lx, ymax + 0.04 * (ymax - ymin), note, ha="right", va="bottom", color=GRAY)
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")


def figure(name, panels, left, right, walls, ax_w=6.3, side=False):
    """panels: list of (geom, dims, title, note); left/right: face BC text; walls: wall BC text.
    side=True: panels in a left column of width ax_w, all text in a column to their right."""
    if side:
        return figure_side(name, panels, left, right, walls, ax_w)
    heights, exs = [], []
    for g, *_ in panels:
        def hgt(ex, g=g):
            _, xl, yl = extent(g, ex)
            return ax_w * (yl[1] - yl[0]) / (xl[1] - xl[0])
        ex = 1.0
        while hgt(ex) < 1.0 and ex < 8:
            ex += 1.0
        exs.append(ex)
        heights.append(hgt(ex))
    top_in, bot_in = 0.45, 0.55 + 0.20 * len(panels)
    H = top_in + sum(heights) + bot_in
    fig = plt.figure(figsize=(FIG_W, H))
    y = H - top_in
    for (g, dims, title, note), h, ex in zip(panels, heights, exs):
        y -= h
        ax = fig.add_axes((0.1 / FIG_W, y / H, ax_w / FIG_W, h / H))
        nt = "; ".join(x for x in (note, f"vertical scale ×{ex:g}" if ex > 1 else "") if x)
        panel(ax, g, ex, dims, title, nt or None)
    fig.text(0.1 / FIG_W, 1 - 0.06 / H, left, ha="left", va="top", linespacing=1.35)
    fig.text(1 - 0.1 / FIG_W, 1 - 0.06 / H, right, ha="right", va="top", linespacing=1.35)
    fig.text(0.1 / FIG_W, (bot_in - 0.08) / H, walls, ha="left", va="top", linespacing=1.35)
    for k, (g, *_) in enumerate(panels):
        pre = f"({chr(97 + k)}) " if len(panels) > 1 else ""
        fig.text(0.1 / FIG_W, (bot_in - 0.30 - 0.20 * k) / H, pre + mesh_line(g), ha="left", va="top",
                 color=GRAY)
    bad = check_figure(fig, min_pt=7.0, placed_width_in=FIG_W)
    for b in bad:
        print(f"  FAIL {name}: {b}")
    FAIL.extend(bad)
    fig.savefig(OUT / f"{name}.png", dpi=250)
    plt.close(fig)
    print(f"{name}: {len(bad)} problems, {FIG_W} x {H:.2f} in")


def figure_side(name, panels, left, right, walls, ax_w):
    heights, exs = [], []
    for g, *_ in panels:
        _, xl, yl = extent(g, 1.0)
        heights.append(ax_w * (yl[1] - yl[0]) / (xl[1] - xl[0]))
        exs.append(1.0)
    meshes = [(f"({chr(97 + k)}) " if len(panels) > 1 else "") + mesh_line(g, br=True)
              for k, (g, *_) in enumerate(panels)]
    notes = [(f"({chr(97 + k)}) " if len(panels) > 1 else "") + nt
             for k, (_, _, _, nt) in enumerate(panels) if nt]
    text = "\n\n".join(notes + [left, right, walls])
    n_lines = text.count("\n") + 1 + sum(m.count("\n") + 2 for m in meshes)
    text_h = n_lines * FS * 1.40 / 72 + 0.2
    H = max(sum(heights) + 0.1, text_h)
    fig = plt.figure(figsize=(FIG_W, H))
    y = H - 0.05 - (H - 0.1 - sum(heights)) / 2
    for k, ((g, dims, title, _), h, ex) in enumerate(zip(panels, heights, exs)):
        y -= h
        ax = fig.add_axes((0.05 / FIG_W, y / H, ax_w / FIG_W, h / H))
        tag = f"({chr(97 + k)}) " if len(panels) > 1 else ""
        panel(ax, g, ex, dims, (tag + (title or "")) or None, None)
    x0 = (ax_w + 0.25) / FIG_W
    t = fig.text(x0, 1 - 0.06 / H, text, ha="left", va="top", linespacing=1.35)
    fig.canvas.draw()
    bb = t.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.text(x0, (bb.y0 - 0.12) / H, "\n\n".join(meshes), ha="left", va="top", color=GRAY,
             linespacing=1.35)
    bad = check_figure(fig, min_pt=7.0, placed_width_in=FIG_W)
    for b in bad:
        print(f"  FAIL {name}: {b}")
    FAIL.extend(bad)
    fig.savefig(OUT / f"{name}.png", dpi=250)
    plt.close(fig)
    print(f"{name}: {len(bad)} problems, {FIG_W} x {H:.2f} in")


# ------------------------------------------------------------------ labels
PH90 = r"$\phi$: 90° face"
WALL = (r"walls (hatched): regolith; $\phi$: contact angle $\theta$ = 60° (wall energy); "
        r"zero flux of $T$ and $\rho_v$")
TFIX = r"$T$ fixed at $T_0 \mp G L_x/2$"


def fig_type1():
    g = load(GEO / "t1_box600x300um_capR50_x300um_th40_eps0.86um.opts")
    figure("type1", [(g, [("h", 0, 600, 0, r"$L_x$ = 600 µm"), ("v", 0, 300, 560, r"$L_y$ = 300 µm")],
                      None, "cap: arc radius 50 µm at x = 300 µm, 40° start" + "\n" + "(variants: radius 25, 100 µm; x = 100, 500 µm)")],
           r"$x = 0$ (warm end): " + TFIX + "\n" + r"$\rho_v$: zero flux (closed) or fixed" + "\n"
           + r"(res, ls, dh); " + PH90,
           r"$x = L_x$ (cold end): $T$ fixed;" + "\n" + r"$\rho_v$ as at $x = 0$; " + PH90,
           r"wall $y = 0$ (hatched): regolith, $\theta$ = 60°," + "\n" + r"zero flux of $T$ and $\rho_v$" + "\n"
           r"far side $y = L_y$ (gray): zero flux of $T$, $\rho_v$;" + "\n" + r"$\phi$: 90° face",
           ax_w=2.6, side=True)


def fig_type2():
    g1 = load(GEO / "t2_box600x300um_pairLW_th40_eps0.80um.opts")
    g2 = load(GEO / "t2_box900x300um_chain5rand_th40_eps0.80um.opts")
    figure("type2", [(g1, [("h", 0, 600, 0, r"$L_x$ = 600 µm")], "2A, 2B: pair",
                      "pair: 50 µm cap (warm side), 20 µm cap"),
                     (g2, [("h", 0, 900, 0, r"$L_x$ = 900 µm"), ("v", 0, 300, 860, r"$L_y$ = 300 µm")],
                      "2C: chain", "chain: five caps, 20–60 µm")],
           r"$x = 0$: " + TFIX + "\n" + r"(2B–2D) or $T_0$ (2A)" + "\n" + r"$\rho_v$: zero flux (closed); 2A-2: fixed;"
           + "\n" + PH90,
           r"$x = L_x$: $T$ fixed; $\rho_v$ as at $x = 0$;" + "\n" + PH90,
           r"wall $y = 0$: regolith, $\theta$ = 60°, zero flux" + "\n"
           + r"far side $y = L_y$: zero flux; $\phi$: 90° face", ax_w=2.6, side=True)


def fig_type3():
    g = load(GEO / "channel_2D_H125um_L500um_eps0.86um_lens150.opts")
    figure("type3", [(g, [("h", 0, 500, 0, r"$L_x$ = 500 µm"), ("v", 0, 125, 120, r"$H$ = 125 µm")],
                      None, "lens 150 µm, flat start; $H$ = 125 µm (also 60)")],
           r"$x = 0$: $T = T_0$ (isothermal)" + "\n" + r"$\rho_v = \rho_{vs}(T_0)$ (3A, 3B, 3D)" + "\n"
           + r"or zero flux (3C); " + PH90,
           r"$x = L_x$: $T = T_0$" + "\n" + r"$\rho_v = \rho_{vs}(T_0)$ (3A, 3B)" + "\n"
           + r"or zero flux (3C, 3D); " + PH90,
           r"walls (hatched): regolith; $\phi$: contact angle" + "\n" + r"$\theta$ = 30–150° by run (wall energy);"
           + "\n" + r"zero flux of $T$ and $\rho_v$", ax_w=3.3, side=True)


def fig_type4():
    g = load(GEO / "wedge_2D_om10deg_gap100um_eps0.86um.opts")
    figure("type4", [(g, [("h", 0, g["lx"] * UM, 0, r"$L_x$ ≈ 270 µm")], None,
                      r"$\omega$ = 10° shown (2.5, 5, 20°);" + "\n" + "gap 100 µm at the lens")],
           r"$x = 0$ (wide end): $T = T_0$" + "\n" + r"$\rho_v = \rho_{vs}(T_0)$; " + PH90,
           r"$x = L_x$ (narrow end): $T = T_0$" + "\n" + r"$\rho_v = \rho_{vs}(T_0)$; " + PH90,
           r"walls (hatched): regolith; $\phi$: contact angle" + "\n" + r"$\theta$ = 60° (wall energy);"
           + "\n" + r"zero flux of $T$ and $\rho_v$", ax_w=2.4, side=True)


def fig_type5():
    g = load(GEO / "t5_pore800x150um_wt40um_lens120c_eps0.80um.opts")
    figure("type5", [(g, [("h", 0, 800, 0, r"$L_x$ = 800 µm"), ("v", 0, 150, 760, "150 µm")], None,
                      r"throat $w_t$ = 40 µm (also 80) at $x$ = 400 µm; lens 120 µm at the throat"
                      + "\n" + "(also 100 µm toward either end)")],
           r"$x = 0$ (warm end): " + TFIX + "\n" + r"$\rho_v = \rho_{vs}(T(0))$ (ls) or $\rho_{vs}(T_0)$; 5B: zero flux; "
           + PH90,
           r"$x = L_x$: $T$ fixed" + "\n" + r"$\rho_v$ as at $x = 0$; " + PH90,
           WALL)


def fig_type6():
    g1 = load(GEO / "dogbone_2D_L1100um_W150um_wc40um_lens60_eps0.86um.opts")
    panels = [(g1, [("h", 0, 1100, 0, r"$L_x$ = 1100 µm"), ("v", 0, 150, 200, r"$W$ = 150 µm")], None,
               r"neck $w_c$ = 40 µm shown (80, 150), $x$ = 400–700 µm; shoulders 100 µm; lens 150 µm")]
    if SCR and (SCR / "dogbone_L720_wc40.opts").exists():
        g2 = load(SCR / "dogbone_L720_wc40.opts")
        panels.append((g2, [("h", 0, 720, 0, r"$L_x$ = 720 µm"), ("v", 0, 150, 690, r"$W$ = 150 µm")],
                       "proposed (G ≤ 30 K/m)", r"neck 150 µm long; lens 20 µm from the shoulder"))
    figure("type6", panels,
           r"$x = 0$ (warm end): " + TFIX + "\n" + r"$\rho_v = \rho_{vs}(T(0))$ (local saturation); " + PH90,
           r"$x = L_x$ (cold end): $T$ fixed" + "\n" + r"$\rho_v = \rho_{vs}(T(L_x))$; " + PH90,
           WALL.replace("60°", "60° (90° for 6B)"))


def fig_type7():
    g = load(GEO / "type7_mixed_A_2D_L2000um_eps0.86um.opts")
    figure("type7", [(g, [("h", 0, 2000, 0, r"$L_x$ = 2000 µm")], None,
                      "layout A: 3 lenses (60°), 4 caps and 3 grains (40°); 24.6 % ice")],
           r"$x = 0$ (collector): $\rho_v = h_{\mathrm{open}}\,\rho_{vs}(T_0)$, $h_{\mathrm{open}}$ = 0.99" + "\n"
           + r"$T$ fixed ($T_0$; warm for 7C); " + PH90,
           r"$x = L_x$ (dead end): $\rho_v$ zero flux; " + PH90 + "\n" + r"$T$ fixed (warm for 7D) or zero flux (7B)",
           WALL)


def fig_S():
    panels = []
    if SCR and (SCR / "S7_straight.opts").exists():
        g = load(SCR / "S7_straight.opts")
        panels.append((g, [("h", 0, 800, 0, r"$L_x$ = 800 µm"), ("v", 0, 150, 520, r"$W$ = 150 µm")],
                       "S7: face-stability scan", "lens 150 µm, 60° faces"))
    g5 = load(GEO / "throat_2D_L850um_eps0.86um_narrowL_deadend.opts")
    panels.append((g5, [("h", 0, 850, 0, r"$L_x$ = 850 µm")], None,
                   "throats 30 and 90 µm, both plugged with ice"))
    figure("S", panels,
           r"$x = 0$ (collector, warm end): $T$ fixed at $T_0 + G L_x/2$" + "\n"
           + r"$\rho_v = h_{\mathrm{open}}\,\rho_{vs}(T_0)$; " + PH90,
           r"$x = L_x$ (dead end): sealed wall," + "\n" + r"zero flux of $\rho_v$; $T$ fixed at $T_0 - G L_x/2$",
           WALL)


if __name__ == "__main__":
    for f in (fig_type1, fig_type2, fig_type3, fig_type4, fig_type5, fig_type6, fig_type7, fig_S):
        f()
    print("TOTAL problems:", len(FAIL))
    sys.exit(1 if FAIL else 0)
