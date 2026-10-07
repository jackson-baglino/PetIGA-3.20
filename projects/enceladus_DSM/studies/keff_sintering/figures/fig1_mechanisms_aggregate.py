#!/usr/bin/env python3
"""Figure 1 (manuscript), assembled: (a) the two sintering mechanisms on a grain
pair, (b) the same process in an aggregate.

    venv_enceladus/bin/python studies/keff_sintering/figures/fig1_mechanisms_aggregate.py
        <campaign dir> --panel-a-pdf <schematic.pdf> [--copy-to <dir>]
        <campaign dir> --pair-run <grain-pair run> [--schematic-mm 95] [--copy-to <dir>]

(a) is either a finished PDF placed as it is, at its own size and still vector
(--panel-a-pdf: the manuscript uses the hand-finished Inkscape schematic
Figure1__SinteringMechanisms/ice_phi0.5_final_mirrored.pdf, 85 mm wide, which
carries the R(t) label; assembled with pdflatex, PNG by mutool), or
postprocess/plot_mechanisms_figure.py's schematic drawn by its own build()
(--pair-run). (b) is fig1_aggregate_strip.py's
close-up of the master run (phi 0.325, seed 1702, -20 C) at four instants.
170 mm wide. Writes Figure1_mechanisms_aggregate.{pdf,png}.
"""
from __future__ import annotations

import argparse, shutil, subprocess, sys, tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
sys.path.insert(0, str(HERE))
import pplib  # noqa: E402
import plot_mechanisms_figure as mech  # noqa: E402
from fig1_aggregate_strip import add_strip  # noqa: E402
from plot_keff_snapshots import INK, FS  # noqa: E402

MM = 1 / 25.4


def assemble_pdf(a, out):
    """Panel (a) from a finished PDF, placed at its own size above the strip."""
    box = subprocess.run(["mutool", "info", str(a.panel_a_pdf)], capture_output=True, text=True).stdout
    nums = box.split("[", 1)[1].split("]", 1)[0].split()
    wa, ha = (float(nums[2]) - float(nums[0])) * 25.4 / 72, (float(nums[3]) - float(nums[1])) * 25.4 / 72
    W, gap_ab, title, bot = 170.0, 3.0, 6.0, 1.0
    cell = (W - 2 * 1.0 - 3 * 2.0) / 4
    hs = cell + title + bot
    H = ha + gap_ab + hs
    fig = plt.figure(figsize=(W * MM, hs * MM))
    add_strip(fig, a.root, a.center_um, a.window_um, W, hs, y_base=bot)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        fig.savefig(td / "strip.pdf", dpi=600, transparent=True)
        shutil.copyfile(a.panel_a_pdf, td / "a.pdf")
        lab = rf"\fontsize{{{FS}}}{{{FS}}}\selectfont\textbf"
        (td / "f.tex").write_text(rf"""\documentclass{{article}}
\usepackage[paperwidth={W}mm,paperheight={H:.3f}mm,margin=0mm]{{geometry}}
\usepackage{{graphicx}}\pagestyle{{empty}}\setlength{{\unitlength}}{{1mm}}\setlength{{\parindent}}{{0pt}}
\begin{{document}}\begin{{picture}}({W},{H:.3f})
\put({(W - wa) / 2:.3f},{H - ha:.3f}){{\includegraphics[width={wa:.3f}mm]{{a.pdf}}}}
\put(0,0){{\includegraphics[width={W}mm]{{strip.pdf}}}}
\put(1,{H - 5:.3f}){{{lab}{{(a)}}}}
\put(1,{bot + cell + title - 0.5:.3f}){{{lab}{{(b)}}}}
\end{{picture}}\end{{document}}
""")
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "f.tex"], cwd=td, capture_output=True, text=True)
        if not (td / "f.pdf").is_file():
            sys.exit("pdflatex failed:\n" + r.stdout[-1500:])
        pdf = out / "Figure1_mechanisms_aggregate.pdf"
        shutil.copyfile(td / "f.pdf", pdf)
    png = pdf.with_suffix(".png")
    subprocess.run(["mutool", "draw", "-q", "-c", "rgba", "-r", "600", "-o", str(png), str(pdf), "1"], check=True)
    if a.copy_to:
        a.copy_to.mkdir(parents=True, exist_ok=True)
        for f in (pdf, png):
            shutil.copyfile(f, a.copy_to / f.name)
    print(f"wrote {pdf} / .png ({W:.0f} x {H:.0f} mm; panel (a) {wa:.1f} x {ha:.1f} mm from {a.panel_a_pdf.name})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--pair-run", type=Path, default=None)
    ap.add_argument("--panel-a-pdf", type=Path, default=None)
    ap.add_argument("--schematic-mm", type=float, default=95.0)
    ap.add_argument("--center-um", type=float, nargs=2, default=[700.0, 1450.0])
    ap.add_argument("--window-um", type=float, default=420.0)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--copy-to", type=Path, default=None)
    a = ap.parse_args()
    plt.rcParams.update(pplib.MANUSCRIPT_RC)

    out = a.out or a.root / "compare" / "figure_samples"
    out.mkdir(parents=True, exist_ok=True)
    if a.panel_a_pdf is not None:
        return assemble_pdf(a, out)
    if a.pair_run is None:
        ap.error("give --panel-a-pdf or --pair-run")

    # the schematic with its own defaults, at the requested width
    ma = mech.build_parser().parse_args(["--run", str(a.pair_run), "--width-mm", str(a.schematic_mm)])
    fig = mech.build(a.pair_run.resolve(), ma)
    hm = fig.get_figheight() / MM
    W, gap_ab, title, bot = 170.0, 3.0, 6.0, 1.0
    cell = (W - 2 * 1.0 - 3 * 2.0) / 4
    H = hm + gap_ab + title + cell + bot
    fig.set_size_inches(W * MM, H * MM)
    fig.axes[0].set_position([(W - a.schematic_mm) / 2 / W, (H - hm) / H, a.schematic_mm / W, hm / H])
    add_strip(fig, a.root, a.center_um, a.window_um, W, H, y_base=bot)
    fig.text(1.0 / W, (H - 4.0) / H, pplib.bold("(a)"), fontsize=FS, color=INK, ha="left", va="center")
    fig.text(1.0 / W, (bot + cell + title + 0.5) / H, pplib.bold("(b)"), fontsize=FS, color=INK,
             ha="left", va="center")

    for e in ("pdf", "png"):
        f = out / f"Figure1_mechanisms_aggregate.{e}"
        fig.savefig(f, dpi=600, transparent=True)
        if a.copy_to:
            a.copy_to.mkdir(parents=True, exist_ok=True); shutil.copyfile(f, a.copy_to / f.name)
    print(f"wrote {out}/Figure1_mechanisms_aggregate.pdf/.png ({W:.0f} x {H:.0f} mm)")


if __name__ == "__main__":
    main()
