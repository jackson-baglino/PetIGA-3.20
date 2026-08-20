#!/usr/bin/env bash
# =============================================================================
# run_postprocess.sh — Run all post-processing on a completed simulation folder
#
# This script is copied into every run folder by stage_output_folder(), so it
# can be run standalone after transferring HPC results to a local machine.
#
# Usage (from inside the run folder):
#   bash postprocess/run_postprocess.sh
#
# Or with an explicit run folder path:
#   bash /path/to/run_folder/postprocess/run_postprocess.sh /path/to/run_folder
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve the run folder: either the argument, or two levels above this script
# (postprocess/ sits one level below the run folder root)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${1:-$(cd "$SCRIPT_DIR/.." && pwd)}"

if [ ! -d "$RUN_DIR" ]; then
    echo "❌ Run folder not found: $RUN_DIR"
    exit 1
fi

POSTPROCESS_DIR="$RUN_DIR/postprocess"
if [ ! -d "$POSTPROCESS_DIR" ]; then
    echo "❌ postprocess/ directory not found inside run folder: $RUN_DIR"
    exit 1
fi

echo ""
echo "========================================================================="
echo "  Phase-field post-processing"
echo "  Run folder : $RUN_DIR"
echo "========================================================================="

# ---------------------------------------------------------------------------
# Detect Python
# ---------------------------------------------------------------------------
PYTHON=$(command -v python3 2>/dev/null || command -v python 2>/dev/null || echo "")
if [[ -z "$PYTHON" ]]; then
    echo "❌ python3 not found — cannot run post-processing."
    exit 1
fi

# ---------------------------------------------------------------------------
# Detect simulation dimension from any .opts file in the run folder
# ---------------------------------------------------------------------------
dim=2
for f in "$RUN_DIR"/*.opts; do
    if [[ -f "$f" ]]; then
        d=$(awk '$1=="-dim"{print $2}' "$f" | head -n1)
        [[ -n "$d" ]] && dim=$d && break
    fi
done
echo "Detected dim = $dim"

overall_exit=0

# Every figure goes here, so the run folder root stays solver output only.
PLOTS="$RUN_DIR/plots"
mkdir -p "$PLOTS"

# ---------------------------------------------------------------------------
# run_step <label> <script> [args...]
#
# Runs one post-processing script and folds its REAL exit status into
# overall_exit. The previous version did `cmd | sed` then `$?`, which reads
# sed's status -- always 0 -- so every failure was silently counted as a pass.
# PIPESTATUS[0] is the actual one.
# ---------------------------------------------------------------------------
run_step() {
    local label="$1" script="$2"; shift 2
    [[ -f "$script" ]] || { echo "⚠️  missing $(basename "$script") — skipping $label"; return; }
    echo ""
    echo "--- $label ---"
    set +e
    "$PYTHON" "$script" "$@" 2>&1 | sed 's/^/  /'
    local rc=${PIPESTATUS[0]}
    set -e
    (( rc != 0 )) && { echo "  ⚠️  $(basename "$script") exited $rc"; (( overall_exit += rc )); }
    return 0
}

# ---------------------------------------------------------------------------
# VTK conversion (2D and 3D only — not meaningful for 1D)
#
# This is a format conversion, not a figure, so it writes to vtkOut/ and is
# deliberately NOT routed into plots/. plot_fields_highres.py is likewise not
# run here: it is slow and only wanted on demand.
# ---------------------------------------------------------------------------
if [[ "$dim" != "1" ]]; then
    if [[ -f "$RUN_DIR/igasol.dat" ]]; then
        mkdir -p "$RUN_DIR/vtkOut"
        run_step "VTK conversion" "$POSTPROCESS_DIR/plot_fields.py" --dir "$RUN_DIR"
    else
        echo "⚠️  igasol.dat not found — skipping VTK conversion."
    fi
fi

# ---------------------------------------------------------------------------
# Scalar time-series from SSA_evo.dat
# ---------------------------------------------------------------------------
if [[ -f "$RUN_DIR/SSA_evo.dat" ]]; then
    run_step "Porosity + interface density" \
        "$POSTPROCESS_DIR/plot_porosity.py" --dir "$RUN_DIR" --save "$PLOTS/porosity.png"
    run_step "Ice-air surface area" \
        "$POSTPROCESS_DIR/plot_ssa.py" --dir "$RUN_DIR" --save "$PLOTS/ssa.png"
else
    echo "⚠️  SSA_evo.dat not found — skipping porosity and surface-area plots."
fi

# ---------------------------------------------------------------------------
# Time step diagnostic (outp.txt)
# ---------------------------------------------------------------------------
if [[ -f "$RUN_DIR/outp.txt" ]]; then
    run_step "Time step diagnostic" \
        "$POSTPROCESS_DIR/plot_timestep.py" --dir "$RUN_DIR" --save "$PLOTS/timestep.png"
fi

# ---------------------------------------------------------------------------
# Phase mass vs. time (igasol.dat + sol_*.dat required)
# ---------------------------------------------------------------------------
if [[ -f "$RUN_DIR/igasol.dat" ]] && ls "$RUN_DIR"/sol_*.dat &>/dev/null 2>&1; then
    run_step "Phase mass vs. time" \
        "$POSTPROCESS_DIR/plot_mass.py" --dir "$RUN_DIR" \
        --save "$PLOTS/mass.png" --per-phase-dir "$PLOTS/mass"
fi

# ---------------------------------------------------------------------------
# Two-grain sintering: neck width, per-grain radii, and the Molaro comparison.
#
# Only for runs that actually have a grain PAIR -- these scripts locate a neck
# between two centres, so they are meaningless on a packing or a single grain.
# Detected from -ice_grain_cx having exactly two entries.
# ---------------------------------------------------------------------------
n_grains=$(awk '$1=="-ice_grain_cx"{n=split($2,a,","); print n; exit}' \
           "$RUN_DIR"/*.opts 2>/dev/null | head -n1)
axisym=$(awk '$1=="-axisym"{print $2; exit}' "$RUN_DIR"/*.opts 2>/dev/null | head -n1)

if [[ "${n_grains:-0}" == "2" ]] && compgen -G "$RUN_DIR/vtkOut/solV_*.vts" >/dev/null; then
    ax_flag=""; [[ "${axisym:-0}" == "1" ]] && ax_flag="--axisym"
    run_step "Neck width" \
        "$POSTPROCESS_DIR/neck_width.py" "$RUN_DIR" $ax_flag
    run_step "Per-grain radii" \
        "$POSTPROCESS_DIR/grain_shrinkage.py" "$RUN_DIR"
    # Needs neck_width.csv from the step above; skips itself with a clear
    # message if that failed (e.g. a truncated outp.txt gives no step->time map).
    if [[ -f "$RUN_DIR/neck_width.csv" ]]; then
        run_step "Neck growth vs Molaro 2019" \
            "$POSTPROCESS_DIR/plot_neck_vs_molaro.py" "$RUN_DIR"
    fi
fi

echo ""
echo "========================================================================="
if [[ "$overall_exit" -ne 0 ]]; then
    echo "  ⚠️  Post-processing completed with errors (exit sum $overall_exit)"
else
    echo "  ✅ Post-processing complete"
fi
echo "  Results in: $RUN_DIR"
echo "  Figures in: $PLOTS"
if compgen -G "$PLOTS/*.png" >/dev/null; then
    (cd "$PLOTS" && ls *.png) | sed 's/^/      /'
fi
echo "========================================================================="
echo ""

exit "$overall_exit"
