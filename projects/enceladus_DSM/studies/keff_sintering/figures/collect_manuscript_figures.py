#!/usr/bin/env python3
"""Collect every main-text figure into ONE flat folder, for upload to Overleaf.

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
    ("Figure7_state_law_alt", SAMPLES),
    ("Figure8_timescales", SAMPLES),
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
print(f"{a.dest}:")
for f in sorted(a.dest.iterdir()):
    print(f"  {f.name:42s} {f.stat().st_size / 1e6:6.2f} MB")
