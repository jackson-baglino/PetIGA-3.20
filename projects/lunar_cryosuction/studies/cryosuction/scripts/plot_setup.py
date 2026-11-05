"""Setup figures for the stage-0 cryosuction runs: geometry, boundary
conditions and initial ice, drawn from the solver's own initial state.

The ice is the solver's phi at t = 1e-9 s (one 1 ns step from the IC; nothing
moves on that scale), read from vtkOut/solV_00000.vts, so the picture is what
the runs start from, not a redrawing of the opts. Make the IC folders with:

    folder=<dir> ./lunar_cryosuction -options_file inputs/solver.opts \
        -options_file <geometry>.opts -options_file <experiment>.opts \
        -t_final 1e-9 -n_out 2 -outp 1 -pf_output 1 -dtmax 1e-9 -delt_t 1e-9
    venv_cryosuction/bin/python postprocess/plot_fields.py --dir <dir>

(-t_final 0 with output on hangs: monitoring.c advances the next output time
by t_interv = 0 forever.)

    venv_cryosuction/bin/python studies/cryosuction/scripts/plot_setup.py \
        --ic-root <folder holding chan/ throat/ regpore/ off/ int/ cap/> --out <dir>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))
from pplib import read_vts  # noqa: E402

# Reference palette (dataviz skill, references/palette.md), light mode.
ICE = "#2a78d6"
OPEN = "#eb6834"
DEAD = "#6250d6"
ROCK = "#b9b8b2"
WALL = "#52514e"
PORE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
UM = 1e6


def load(ic_root: Path, name: str):
    f = sorted((ic_root / name / "vtkOut").glob("solV_*.vts"))[0]
    fields, X, Y = read_vts(str(f), want=["IcePhase"])
    return fields["IcePhase"], X * UM, Y * UM


def panel(ax, phi, X, Y, title, faces, notes):
    """faces: dict side -> ('open'|'dead'|'wall', label) for x0, x1, y0, y1."""
    ax.set_facecolor(ROCK)
    ax.pcolormesh(X, Y, (phi > 0.5).astype(float), cmap=ListedColormap([PORE, ICE]),
                  shading="gouraud", vmin=0, vmax=1, rasterized=True)
    ax.contour(X, Y, phi, levels=[0.5], colors=[INK], linewidths=0.8)
    edges = {"x0": (X[:, 0], Y[:, 0]), "x1": (X[:, -1], Y[:, -1]),
             "y0": (X[0, :], Y[0, :]), "y1": (X[-1, :], Y[-1, :])}
    style = {"open": dict(color=OPEN, lw=4, ls=(0, (4, 2))),
             "dead": dict(color=DEAD, lw=4, ls="-"),
             "wall": dict(color=WALL, lw=1.6, ls="-"),
             "plain": dict(color=INK2, lw=1.0, ls=":")}
    for side, (kind, _label) in faces.items():
        ex, ey = edges[side]
        ax.plot(ex, ey, **style[kind], solid_capstyle="butt")
    ax.set_aspect("equal")
    ax.set_xlim(X.min() - 0.02 * np.ptp(X), X.max() + 0.02 * np.ptp(X))
    ax.set_ylim(Y.min() - 0.04 * np.ptp(X), Y.max() + 0.04 * np.ptp(X))
    ax.set_title(title, loc="left", fontsize=10.5, color=INK, fontweight="bold", pad=16)
    ax.set_xlabel("x [µm]", fontsize=8.5, color=INK2)
    ax.set_ylabel("y [µm]", fontsize=8.5, color=INK2)
    ax.tick_params(labelsize=7.5, colors=INK2)
    for s in ax.spines.values():
        s.set_visible(False)
    # Face labels above the panel at each end, so they never meet the tick labels.
    ytop = Y.max() + 0.04 * np.ptp(X)
    if "x0" in faces and faces["x0"][1]:
        ax.text(X[:, 0].mean(), ytop, faces["x0"][1], ha="left", va="bottom", fontsize=8,
                color=OPEN if faces["x0"][0] == "open" else INK, fontweight="bold")
    if "x1" in faces and faces["x1"][1]:
        ax.text(X[:, -1].mean(), ytop, faces["x1"][1], ha="right", va="bottom", fontsize=8,
                color=DEAD, fontweight="bold")
    if notes:
        ax.annotate(notes, xy=(0.0, 0.0), xycoords="axes fraction", xytext=(0, -42),
                    textcoords="offset points", fontsize=8, color=INK2, va="top")


def legend(fig, y):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h = [Patch(facecolor=ICE, edgecolor=INK, lw=0.8, label="ice (φ > 0.5, solver IC)"),
         Patch(facecolor=PORE, edgecolor=INK2, lw=0.5, label="pore: air + vapor"),
         Patch(facecolor=ROCK, label="regolith (outside the domain)"),
         Line2D([], [], color=OPEN, lw=4, ls=(0, (4, 2)), label="open end: T, ρv fixed; φ free (ice passes)"),
         Line2D([], [], color=DEAD, lw=4, label="dead pore: regolith, θ = 60°, no vapor flux"),
         Line2D([], [], color=WALL, lw=1.6, label="regolith wall: θ = 60°, no vapor flux")]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, y), ncol=3, frameon=False,
               fontsize=8.5, labelcolor=INK)


def main_figure(ic: Path, out: Path):
    fig = plt.figure(figsize=(11, 12.5))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.25, 1.25], hspace=0.70,
                          left=0.08, right=0.94, top=0.84, bottom=0.07)
    fig.suptitle("Cryosuction stage 0: geometry, boundary conditions and initial ice  (T₀ = −20 °C, α_c = 10⁻²)",
                 x=0.08, ha="left", fontsize=13, color=INK, fontweight="bold", y=0.975)
    fig.text(0.08, 0.945, "Open end at x = 0 (left) faces the dry reservoir; x = Lx (right) is a dead pore. "
             "Ice is drawn from the solver's own initial state.\nWalls carry the static contact angle θ = 60° "
             "(Young: γ_ia 0.109, γ_as 0.300, γ_is 0.2455 J/m²).\neps and mesh from comp_eps.py at each run's expected front speed.",
             fontsize=9, color=INK2, va="top")
    legend(fig, 0.915)
    faces = {"x0": ("open", "open end"), "x1": ("dead", "dead pore"),
             "y0": ("wall", ""), "y1": ("wall", "")}

    phi, X, Y = load(ic, "chan")
    panel(fig.add_subplot(gs[0]), phi, X, Y,
          "B1–B4  flat channel 125 × 375 µm, one ice bridge  (eps 0.39 µm)", faces,
          "Open end h = 0.99 of ρ_vs(T₀); every face held at T₀ (isothermal bath); 3 days.\n"
          "B1: ξ_v = 1e-3, ξ_T = 1 (M&F defaults).   B2: ξ_v = 1e-4.   B3: ξ_v = 1e-2.   B4: ξ_v = ξ_T = 1e-3.")
    phi, X, Y = load(ic, "throat")
    panel(fig.add_subplot(gs[1]), phi, X, Y,
          "C1, C3, C4  two-throat channel, 850 µm: plugs in a 30 µm and a 90 µm throat   (eps 0.86 µm)", faces,
          "Open end h = 0.99; 5 days.   C1: T₀ on every face.   C3: T fixed at both ends, +100 K/m "
          "(open end −20.04 °C, dead end −19.96 °C), walls insulated.\n"
          "C4: T fixed at the open end only; dead end and walls insulated (isolates sublimation cooling).")
    phi, X, Y = load(ic, "regpore")
    panel(fig.add_subplot(gs[2]), phi, X, Y,
          "C2  regolith pore, 814 µm, walls of 12 + 12 grains: central grain + 6 wall caps   (eps 0.86 µm)",
          faces, "Open end h = 0.998 (h = 0.99 would need eps 0.20 µm); T₀ on every face; 5 days.\n"
          "Note: the builder's 'throat_bridge' grain does not touch either wall (free disc).")
    fig.savefig(out / "stage0_setup_cryosuction.png", dpi=200, facecolor="white")
    plt.close(fig)


def verification_figure(ic: Path, out: Path):
    fig = plt.figure(figsize=(11, 6.2))
    gs = fig.add_gridspec(1, 3, wspace=0.35, left=0.07, right=0.97, top=0.72, bottom=0.30)
    fig.suptitle("Stage 0, tests A: does ice pass through the open end undisturbed?",
                 x=0.07, ha="left", fontsize=13, color=INK, fontweight="bold", y=0.97)
    fig.text(0.07, 0.915, "Pure curvature flow (phase change off): the ice area obeys dA/dt = −D·(total turning "
             "of the interface), D = eps²/τ_sub = 1.14e-15 m²/s (eps 0.64 µm, α_c 10⁻²).\nAn open face must let the cut part shrink "
             "exactly as if the domain continued. 1 rank each, 22 hours.", fontsize=9, color=INK2, va="top")
    legend(fig, 0.83)
    plain = {"x1": ("plain", ""), "y0": ("plain", ""), "y1": ("plain", "")}
    phi, X, Y = load(ic, "off")
    panel(fig.add_subplot(gs[0]), phi, X, Y, "A1 / A2  disc cut by x = 0",
          {"x0": ("open", "open (A1)  /  90° wall (A2)"), **plain},
          "Center 7.5 µm outside, R₀ = 15 µm.\nA1 open: must follow the full-disc law.\n"
          "A2 Neumann 90°: control, should miss.")
    phi, X, Y = load(ic, "int")
    panel(fig.add_subplot(gs[1]), phi, X, Y, "A3  reference disc",
          {"x0": ("open", "open"), **plain}, "Same disc fully inside:\nmeasures D (dA/dt = −2πD).")
    phi, X, Y = load(ic, "cap")
    panel(fig.add_subplot(gs[2]), phi, X, Y, "A4  wall cap at the corner",
          {"x0": ("open", "open"), "y0": ("wall", ""), "x1": ("plain", ""), "y1": ("plain", "")},
          "θ = 60° cap on wall y = 0 cut by x = 0:\ncontact line where wall and\nopen end meet.")
    fig.text(0.07, 0.03, "Dotted edges: plain 90° faces (no wall energy), sealed for vapor, T fixed at T₀. "
             "Gray regolith fill is outside the domain.", fontsize=8, color=INK2)
    fig.savefig(out / "stage0_setup_verification.png", dpi=200, facecolor="white")
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ic-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    print(f"# IC folders under {a.ic_root}; figures to {a.out}")
    main_figure(a.ic_root, a.out)
    verification_figure(a.ic_root, a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
