#!/usr/bin/env python3
"""Collect every main-text figure, and tables.tex, for upload to Overleaf: ONE
subfolder per figure iteration, each holding ALL of Figures 1-7.

    python3 studies/keff_sintering/figures/collect_manuscript_figures.py [--dest <dir>]

Each iteration folder <dest>/<iteration>/ gets the figures current at that
iteration under their CANONICAL names (FigureN_<name>.pdf/.png, no _vN
suffix), so the folder alone says which version they are and Overleaf paths
never change. Alternatives that were offered alongside keep their own names.
NOTES.txt in each folder says what changed. Sources are where each script
writes (the campaign's compare/figure_samples/ and studies/molaro_2019/manuscript/);
nothing is rebuilt -- run the figure scripts first. tables.tex and
supplement/ stay at the top of <dest>. An iteration folder that already exists
is NEVER rewritten (the user prunes them by hand) unless named in --rebuild.
"omit" drops a canonical figure from a round (e.g. the all-phi Figure 7). Loose FigureN files left at the top by
the old flat layout are removed (they are copies; the sources stay).
"""
import argparse, shutil
from pathlib import Path

HPC = Path.home() / "SimulationResults/HPC_results/enceladus_DSM"
SAMPLES = HPC / "keff_sintering_campaign/compare/figure_samples"
MOLARO = Path(__file__).resolve().parents[3] / "studies/molaro_2019/manuscript"
CANON = [("Figure1_mechanisms_aggregate", SAMPLES), ("Figure2_homogenization_method", SAMPLES),
         ("Figure3_molaro_validation", MOLARO), ("Figure4_saturated_neck_growth", MOLARO),
         ("Figure5_gallery", SAMPLES), ("Figure6_keff_collapse", SAMPLES),
         ("Figure7_state_law", SAMPLES), ("Figure7_state_law_allphi", SAMPLES)]
# iteration -> (suffix of the version used per canonical stem ("" = original),
#               extra alternatives kept under their own names, notes)
ITERATIONS = {
    "v1_2026-10-15_original": ({}, [], "The figures before the 2026-10-15 advisor meeting."),
    "v2_2026-10-15": (
        {"Figure1_mechanisms_aggregate": "_v2", "Figure2_homogenization_method": "_v2",
         "Figure3_molaro_validation": "_v2", "Figure4_saturated_neck_growth": "_v2",
         "Figure7_state_law": "_v2", "Figure7_state_law_allphi": "_v2"},
        ["Figure1_mechanisms_aggregate_v2_dotsA", "Figure1_mechanisms_aggregate_v2_dotsB",
         "Figure1_mechanisms_aggregate_v2_dotsC"],
        "Advisor meeting: Fig 1 molecules + w(t) I-bar; Fig 2 scale bar in ice, (c) = T(x,y); "
        "Figs 3-4 'simulation', 'theory', 'Data from Molaro et al.', w(t) I-bars, lighter sigma map "
        "(Fig 3); Fig 7(b) as a scatter. dotsA/B/C: denser, smaller molecules."),
    "v3_2026-10-16": (
        {"Figure1_mechanisms_aggregate": "_v3a", "Figure2_homogenization_method": "_v3",
         "Figure3_molaro_validation": "_v3", "Figure4_saturated_neck_growth": "_v2",
         "Figure7_state_law": "_v2", "Figure7_state_law_allphi": "_v2"},
        ["Figure1_mechanisms_aggregate_v3b"],
        "Fig 1 molecule band from the grain surface to the neck, grey outline (v3b: more, smaller); "
        "Fig 2(c) no isotherms, RdYlBu_r; Fig 3 original balance map with sigma below the bar grey."),
    "v4_2026-10-16": (
        {"Figure1_mechanisms_aggregate": "_v4", "Figure2_homogenization_method": "_v4_bluered",
         "Figure3_molaro_validation": "_v4", "Figure4_saturated_neck_growth": "_v2",
         "Figure7_state_law": "_v2", "Figure7_state_law_allphi": "_v2"},
        ["Figure2_homogenization_method_v4_amp", "Figure3_molaro_validation_v4_fullrange"],
        "Fig 1 molecules on the arrow's trajectory, crowded where they leave the grain and reattach "
        "at the neck; Fig 2(c) sequential two-colour map (alternative: _v4_amp, one hue); Fig 3 back "
        "to the original sigma map, v2 labels kept (proposal: _v4_fullrange, sigma bar over the pore's "
        "full range)."),
    "v5_2026-10-16": (
        {"Figure1_mechanisms_aggregate": "_v5", "Figure2_homogenization_method": "_v4_bluered",
         "Figure3_molaro_validation": "_v4_fullrange_trim", "Figure4_saturated_neck_growth": "_v4_amp",
         "Figure7_state_law": "_v2", "omit": ["Figure7_state_law_allphi"]},
        [],
        "Finalised: Fig 2 (two-colour T map) and Fig 3 (full-range sigma bar). Fig 1 molecules on both "
        "sides of the arrow, 'Vapor Transport' moved left to make room. Fig 3 simulation curves stop at "
        "each series' last data point. Fig 4 sigma on cmocean amp from sigma = 0 (white, always on the "
        "bar) to the max, gamma = 3 power norm so the depletion halo at the neck shows (the far field is "
        "flat). Fig 7 without the all-phi variant."),
    "v6_2026-10-20": (
        {"Figure1_mechanisms_aggregate": "_v5", "Figure2_homogenization_method": "_v4_bluered",
         "Figure3_molaro_validation": "_v4_fullrange_trim", "Figure4_saturated_neck_growth": "_v4_amp_fit",
         "Figure7_state_law": "_v2", "omit": ["Figure7_state_law_allphi"]},
        [],
        "Fig 4 fitted lines restored: the free power law (a = 0.29) and the one-third theory curve were "
        "missing from v5 (a tick-list variable overwrote the fit window in fig_demmenie.py). Nothing else "
        "changed from v5."),
    "v7_2026-10-20": (
        {"Figure1_mechanisms_aggregate": "_v5", "Figure2_homogenization_method": "_v4_bluered",
         "Figure3_molaro_validation": "_v4_fullrange_trim", "Figure4_saturated_neck_growth": "_v4_amp_fit_s2",
         "Figure7_state_law": "_v2", "omit": ["Figure7_state_law_allphi"]},
        [],
        "Fig 4 opens on the first resolved sample (1 h, w = 16.8 um): the first two output samples (0 and "
        "1.3 s) are not drawn, since the 1 h gap to the next one is unresolved and a line across it is not "
        "a measurement. Instant 1 is that snapshot (it already has a neck). Its few crevice cells below "
        "sigma = 0 show in blue under a 'below 0' arrow on the sigma bar. Fits unchanged (a = 0.289)."),
}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dest", type=Path, default=HPC / "ManuscriptFigures")
ap.add_argument("--rebuild", nargs="*", default=[], help="iteration folders to rewrite even if present")
a = ap.parse_args()
a.dest.mkdir(parents=True, exist_ok=True)
for it, (use, extras, notes) in ITERATIONS.items():
    d = a.dest / it
    if d.is_dir() and it not in a.rebuild:
        continue
    d.mkdir(exist_ok=True)
    jobs = [(f"{stem}{use.get(stem, '')}", src, stem) for stem, src in CANON if stem not in use.get("omit", [])]
    jobs += [(x, next(src for stem, src in CANON if x.startswith(stem)), x) for x in extras]
    for name, src, out in jobs:
        for ext in ("pdf", "png"):
            f = src / f"{name}.{ext}"
            if f.is_file():
                shutil.copyfile(f, d / f"{out}.{ext}")
            else:
                print(f"MISSING {f}")
    (d / "NOTES.txt").write_text(notes + "\n")
for f in a.dest.glob("Figure*.*"):          # the old flat layout
    f.unlink()
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
