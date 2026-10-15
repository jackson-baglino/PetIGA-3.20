#!/usr/bin/env python3
"""Collect every main-text figure, and tables.tex, into ONE flat folder, for upload to Overleaf.

    python3 studies/keff_sintering/figures/collect_manuscript_figures.py [--dest <dir>]

Copies the current FigureN_<name>.pdf/.png from where each script writes it
(the campaign's compare/figure_samples/ and studies/molaro_2019/manuscript/)
into <dest> (default ~/SimulationResults/HPC_results/enceladus_DSM/ManuscriptFigures).
It rebuilds nothing: run the figure scripts first. Files in <dest> that are
not in the list below are left alone; a missing source is reported.
"""
import argparse, shutil
from pathlib import Path

HPC = Path.home() / "SimulationResults/HPC_results/enceladus_DSM"
SAMPLES = HPC / "keff_sintering_campaign/compare/figure_samples"
MOLARO = Path(__file__).resolve().parents[3] / "studies/molaro_2019/manuscript"
FIGURES = [                                   # (stem, source folder)
    ("Figure1_mechanisms_aggregate", SAMPLES),
    ("Figure2_homogenization_method", SAMPLES),
    ("Figure3_molaro_validation", MOLARO),
    ("Figure4_saturated_neck_growth", MOLARO),
    ("Figure5_gallery", SAMPLES),
    ("Figure6_keff_collapse", SAMPLES),
    ("Figure7_state_law", SAMPLES),
    ("Figure7_state_law_allphi", SAMPLES),     # variation: every porosity in panel (b)
    # v2: the advisor-meeting edits of 2026-10-15, beside the originals
    ("Figure1_mechanisms_aggregate_v2", SAMPLES),   # fig1_v2_edit.py
    ("Figure1_mechanisms_aggregate_v2_dotsA", SAMPLES),  # --dots 16 0.95
    ("Figure1_mechanisms_aggregate_v2_dotsB", SAMPLES),  # --dots 24 0.75
    ("Figure1_mechanisms_aggregate_v2_dotsC", SAMPLES),  # --dots 34 0.6
    ("Figure2_homogenization_method_v2", SAMPLES),  # fig_methods_homogenization.py --v2
    ("Figure3_molaro_validation_v2", MOLARO),       # plot_molaro_validation.py --v2
    ("Figure4_saturated_neck_growth_v2", MOLARO),   # fig_demmenie.py --v2
    ("Figure7_state_law_v2", SAMPLES),              # sample_figures.py --v2
    ("Figure7_state_law_allphi_v2", SAMPLES),
]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dest", type=Path, default=HPC / "ManuscriptFigures")
a = ap.parse_args()
a.dest.mkdir(parents=True, exist_ok=True)
for stem, src in FIGURES:
    for ext in ("pdf", "png"):
        f = src / f"{stem}.{ext}"
        if f.is_file():
            shutil.copyfile(f, a.dest / f.name)
        else:
            print(f"MISSING {f}")
# Supporting-information plots, under their own names, in <dest>/supplement/.
REPO = Path(__file__).resolve().parents[3]
CMP = HPC / "keff_sintering_campaign/compare"
SUPP = [REPO / "studies/keff_sintering/rve_convergence/rve_convergence.png",
        REPO / "studies/keff_sintering/rve_convergence/kxy_vs_L.png",
        REPO / "studies/keff_sintering/eps_sensitivity/eps_sensitivity.png",
        REPO / "studies/keff_sintering/dtmax_check/dtmax_check.png",
        REPO / "studies/keff_sintering/master_curve/master_curve.png",
        CMP / "kinks/kink_check.png", CMP / "kinks/kink_examples.png",
        CMP / "collapse/collapse_check.png",
        SAMPLES / "figA_collapse.pdf", SAMPLES / "figA_collapse.png",
        CMP / "phi_summary/T-20/phi_trends.png", CMP / "phi_summary/T-20/seeds_by_phi.png",
        CMP / "phi_summary/T-20/anisotropy_time.png", CMP / "phi_summary/T-20/ssa_by_phi.png"]
(a.dest / "supplement").mkdir(exist_ok=True)
for f in SUPP:
    if f.is_file():
        shutil.copyfile(f, a.dest / "supplement" / f.name)
    else:
        print(f"MISSING {f}")
# The tables file rides along, so one folder holds everything the manuscript
# takes in. The repo copy (studies/keff_sintering/tables/tables.tex) is the
# source; edit that one, never the copy.
TABLES = REPO / "studies/keff_sintering/tables/tables.tex"
if TABLES.is_file():
    shutil.copyfile(TABLES, a.dest / TABLES.name)
else:
    print(f"MISSING {TABLES}")
print(f"{a.dest}:")
for f in sorted(x for x in a.dest.rglob("*") if x.is_file()):
    print(f"  {str(f.relative_to(a.dest)):44s} {f.stat().st_size / 1e6:6.2f} MB  {__import__('time').strftime('%m-%d %H:%M', __import__('time').localtime(f.stat().st_mtime))}")
