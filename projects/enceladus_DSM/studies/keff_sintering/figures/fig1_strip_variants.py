#!/usr/bin/env python3
"""Figure 1(b) alternatives: which instants, and which close-up, show sintering best.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig1_strip_variants.py
        <campaign dir> [--top-pdf <Figure1 v5 pdf>] [--out <dir>] [--only NAME ...]

The manuscript strip (fig1_aggregate_strip.py) shows the master run (phi 0.325,
seed 1702, -20 C) at its four saved states, 0, 10, 19 and 30 d. Those are evenly
spaced in linear time, and sintering is front-loaded, so the last three look
alike. The -20 C run saved no other states.

The same packing was also run at -5, -10, -30 and -40 C. In this model
temperature only rescales time by tau_sub (see project notes; k_eff(SSA) is T
independent), so a state of the T run at time t is the -20 C state at
t * tau_sub(-20) / tau_sub(T). Checked on the phase field itself
(2026-11-18): two runs at matched equivalent age differ in 0.01-0.06 % of
pixels (T-40 vs T-30 at 3.9 d, T-5 vs T-10 at 74 d), against ~2 % between
consecutive -20 C frames. Those runs fill in the ages the -20 C run skipped:

    T-40: 0, 1.3, 2.5, 3.9 d    T-30: 3.9, 7.3, 11.3 d   T-20: 10.2, 19.5, 30.0 d
    T-10: 25.0, 48.1, 74.2 d    T-5: 38.3, 73.6, 113.8 d    (-20 C equivalent days)

Each variant is panel (a) of the current Figure 1 (cropped, untouched, from
the hand-edited v5 PDF) over a new strip. Writes
Figure1_mechanisms_aggregate_v8_<name>.{pdf,png} and a contact sheet
Figure1_v8_variants.png. Labels are -20 C equivalent days; a variant that uses
another temperature needs a caption sentence saying so.
"""
from __future__ import annotations

import argparse, shutil, subprocess, sys, tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cmocean

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
import pplib  # noqa: E402
from pplib import step_times  # noqa: E402
from plot_keff import read_tau_sub  # noqa: E402
from plot_keff_snapshots import make_reader, snap_step, INK, MUTED, FS, FS_TINY  # noqa: E402

MM, DAY = 1 / 25.4, 86400.0
RUN = "packing_2D_phi0.325_Rave50um_LR40_seed1702_L2mm_eps1000nm_perxy_T{T}__snow_T{T}_h1.00_30d"
CUR_C, BEST_C = (700.0, 1450.0), (1740.0, 1700.0)   # BEST: most change 10->19->30 d (420 um window)
W_FIG, GAP, LEFT, TITLE, BOT = 170.0, 2.0, 1.0, 6.0, 1.0
A_CUT_MM = 51.5      # v5 PDF: panel (a) ink starts at 52.9 mm, its old (b) label ends at 50.3 mm

# name -> (frames as (T, snapshot index), window centre [um], draw the 0 d outline, note)
VARIANTS = {
    "current":   ([(-20, 0), (-20, 1), (-20, 2), (-20, 3)], CUR_C, False,
                  "the manuscript strip as it is: -20 C, linear 0/10/19/30 d"),
    "log5":      ([(-20, 0), (-40, 1), (-30, 1), (-20, 1), (-20, 3)], CUR_C, False,
                  "log spaced: 0, 1.3, 3.9, 10, 30 d (x2.9 apart); 1.3 d from -40 C, 3.9 d from -30 C"),
    "log4":      ([(-20, 0), (-30, 1), (-20, 1), (-20, 3)], CUR_C, False,
                  "log spaced, four frames: 0, 3.9, 10, 30 d; 3.9 d from -30 C"),
    "log5_best": ([(-20, 0), (-40, 1), (-30, 1), (-20, 1), (-20, 3)], BEST_C, False,
                  "log5 in the close-up with the most late change"),
    "window":    ([(-20, 0), (-20, 1), (-20, 2), (-20, 3)], BEST_C, False,
                  "-20 C only, linear 0/10/19/30 d, close-up moved to where late change is largest"),
    "outline":   ([(-20, 0), (-20, 1), (-20, 2), (-20, 3)], BEST_C, True,
                  "as 'window', with the 0 d ice outline drawn on every later frame"),
    "long":      ([(-20, 0), (-20, 1), (-20, 3), (-5, 3)], CUR_C, False,
                  "past 30 d: 0, 10, 30, 114 d (last from -5 C, its 30 d)"),
}


class Frames:
    """Snapshots of the seed-1702 packing at any temperature, with -20 C equivalent ages."""

    def __init__(self, root):
        self.root, self.cache = root, {}
        self.tau20 = read_tau_sub(root / RUN.format(T=-20))

    def get(self, T, k):
        if (T, k) not in self.cache:
            run = self.root / RUN.format(T=T)
            files, reader = make_reader(run, "sol")
            tm = step_times(str(run))
            fl, X, Y = reader(files[k], want=("IcePhase",))
            age = tm[snap_step(files[k])] * self.tau20 / read_tau_sub(run)
            self.cache[(T, k)] = (fl["IcePhase"], X[0, :] * 1e6, Y[:, 0] * 1e6, age)
        return self.cache[(T, k)]


def label(age):
    d = age / DAY
    return "0 d" if d < 0.05 else (f"{d:.1f} d" if d < 9.5 else f"{d:.0f} d")


def add_strip(fig, frames, spec, H, window=420.0):
    sel, (cx, cy), outline, _ = spec
    n = len(sel)
    cell = (W_FIG - 2 * LEFT - (n - 1) * GAP) / n
    h = window / 2
    phi0 = None
    for j, (T, k) in enumerate(sel):
        phi, x, y, age = frames.get(T, k)
        mx = (x >= cx - h) & (x <= cx + h); my = (y >= cy - h) & (y <= cy + h)
        ext = (x[mx][0], x[mx][-1], y[my][0], y[my][-1])
        ax = fig.add_axes([(LEFT + j * (cell + GAP)) / W_FIG, BOT / H, cell / W_FIG, cell / H])
        ax.imshow(phi[np.ix_(my, mx)], origin="lower", cmap=cmocean.cm.ice, vmin=0, vmax=1,
                  extent=ext, interpolation="antialiased")
        if outline and j == 0:
            phi0 = (phi[np.ix_(my, mx)], x[mx], y[my])
        elif outline:
            ax.contour(phi0[1], phi0[2], phi0[0], levels=[0.5], colors="#e8743b", linewidths=0.5)
        ax.set_xlim(ext[:2]); ax.set_ylim(ext[2:])
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6); sp.set_color(MUTED)
        ax.set_title(label(age), fontsize=FS, pad=3)
        if j == 0:
            ax.plot([cx - h + 20, cx - h + 120], [cy - h + 22, cy - h + 22], color=INK, lw=2.2,
                    solid_capstyle="butt")
            ax.text(cx - h + 70, cy - h + 32, r"100 $\mu$m", ha="center", va="bottom", fontsize=FS_TINY,
                    color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    return cell


def build(frames, name, spec, top_pdf, out):
    n = len(spec[0])
    cell = (W_FIG - 2 * LEFT - (n - 1) * GAP) / n
    hs = cell + TITLE + BOT
    fig = plt.figure(figsize=(W_FIG * MM, hs * MM))
    add_strip(fig, frames, spec, hs)
    box = subprocess.run(["mutool", "info", str(top_pdf)], capture_output=True, text=True).stdout
    nums = box.split("Mediaboxes", 1)[1].split("[", 1)[1].split("]", 1)[0].split()
    Wt, Ht = (float(nums[2]) - float(nums[0])) * 25.4 / 72, (float(nums[3]) - float(nums[1])) * 25.4 / 72
    ha = Ht - A_CUT_MM                      # panel (a) block kept from the v5 PDF
    gap = 3.0                               # as v5: 3 mm between (a) and the strip
    H = ha + gap + hs
    lab = rf"\fontsize{{{FS}}}{{{FS}}}\selectfont\textbf"
    stem = f"Figure1_mechanisms_aggregate_v8_{name}"
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        fig.savefig(td / "strip.pdf", dpi=600, transparent=True)
        plt.close(fig)
        shutil.copyfile(top_pdf, td / "a.pdf")
        (td / "f.tex").write_text(rf"""\documentclass{{article}}
\usepackage[paperwidth={W_FIG}mm,paperheight={H:.3f}mm,margin=0mm]{{geometry}}
\usepackage{{graphicx}}\pagestyle{{empty}}\setlength{{\unitlength}}{{1mm}}\setlength{{\parindent}}{{0pt}}
\begin{{document}}\begin{{picture}}({W_FIG},{H:.3f})
\put(0,{H - ha:.3f}){{\includegraphics[width={Wt:.3f}mm,trim=0 {A_CUT_MM}mm 0 0,clip]{{a.pdf}}}}
\put(0,0){{\includegraphics[width={W_FIG}mm]{{strip.pdf}}}}
\put(1,{BOT + cell + TITLE - 0.5:.3f}){{{lab}{{(b)}}}}
\end{{picture}}\end{{document}}
""")
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "f.tex"], cwd=td, capture_output=True, text=True)
        if not (td / "f.pdf").is_file():
            sys.exit("pdflatex failed:\n" + r.stdout[-1500:])
        shutil.copyfile(td / "f.pdf", out / f"{stem}.pdf")
    subprocess.run(["mutool", "draw", "-q", "-c", "rgba", "-r", "600", "-o", str(out / f"{stem}.png"),
                    str(out / f"{stem}.pdf"), "1"], check=True)
    ages = [label(frames.get(*f)[3]) for f in spec[0]]
    print(f"{stem}: {W_FIG:.0f} x {H:.1f} mm; {ages}; {spec[3]}")
    return out / f"{stem}.png"


def contact_sheet(pngs, out):
    """All variants stacked, panel (b) only, for side-by-side reading."""
    from PIL import Image
    ims = []
    for name, p in pngs:
        im = Image.open(p)
        w, hpx = im.size
        ims.append((name, im.crop((0, int(hpx * 0.42), w, hpx))))
    fig, axs = plt.subplots(len(ims), 1, figsize=(8, 2.1 * len(ims)))
    for ax, (name, im) in zip(np.atleast_1d(axs), ims):
        ax.imshow(im); ax.axis("off")
        ax.set_title(f"{name}: {VARIANTS[name][3]}", fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(out / "Figure1_v8_variants.png", dpi=200, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--top-pdf", type=Path, default=None,
                    help="Figure 1 whose panel (a) is kept (default: <root>/compare/figure_samples/"
                         "Figure1_mechanisms_aggregate_v5.pdf)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)
    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    top = a.top_pdf or a.root / "compare" / "figure_samples" / "Figure1_mechanisms_aggregate_v5.pdf"
    frames = Frames(a.root)
    pngs = [(n, build(frames, n, VARIANTS[n], top, out)) for n in (a.only or VARIANTS)]
    contact_sheet(pngs, out)
    print(f"contact sheet: {out / 'Figure1_v8_variants.png'}")


if __name__ == "__main__":
    main()
