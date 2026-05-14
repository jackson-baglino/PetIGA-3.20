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
echo "  Permafrost post-processing"
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

# ---------------------------------------------------------------------------
# VTK conversion (all dimensions)
# ---------------------------------------------------------------------------
if [[ -f "$RUN_DIR/igasol.dat" ]]; then
    echo ""
    echo "--- VTK conversion ---"
    mkdir -p "$RUN_DIR/vtkOut"
    set +e
    "$PYTHON" "$POSTPROCESS_DIR/plotpermafrost.py" --dir "$RUN_DIR" \
        2>&1 | sed 's/^/  /'
    (( overall_exit += $? )) || true
    set -e
else
    echo "⚠️  igasol.dat not found — skipping VTK conversion."
fi

# ---------------------------------------------------------------------------
# 1D-specific plots
# ---------------------------------------------------------------------------
if [[ "$dim" == "1" ]]; then
    echo ""
    echo "--- 1D phase field profiles ---"
    set +e
    "$PYTHON" "$POSTPROCESS_DIR/plot1D_profiles.py" \
        --dir "$RUN_DIR" --out-dir "$RUN_DIR" \
        2>&1 | sed 's/^/  /'
    (( overall_exit += $? )) || true

    echo ""
    echo "--- 1D derived quantities ---"
    "$PYTHON" "$POSTPROCESS_DIR/plot1D_profiles.py" \
        --dir "$RUN_DIR" --derived --save "$RUN_DIR/derived.png" \
        2>&1 | sed 's/^/  /'
    (( overall_exit += $? )) || true
    set -e
fi

# ---------------------------------------------------------------------------
# Scalar time-series (SSA_evo.dat)
# ---------------------------------------------------------------------------
if [[ -f "$RUN_DIR/SSA_evo.dat" ]]; then
    echo ""
    echo "--- Scalar time-series ---"
    set +e
    "$PYTHON" "$POSTPROCESS_DIR/plot_scalars.py" \
        --file "$RUN_DIR/SSA_evo.dat" --save "$RUN_DIR/scalars.png" \
        2>&1 | sed 's/^/  /'
    (( overall_exit += $? )) || true
    set -e
fi

# ---------------------------------------------------------------------------
# Time step diagnostic (outp.txt)
# ---------------------------------------------------------------------------
if [[ -f "$RUN_DIR/outp.txt" ]]; then
    echo ""
    echo "--- Time step diagnostic ---"
    set +e
    "$PYTHON" "$POSTPROCESS_DIR/plot_timestep.py" \
        --dir "$RUN_DIR" --save "$RUN_DIR/timestep.png" \
        2>&1 | sed 's/^/  /'
    (( overall_exit += $? )) || true
    set -e
fi

# ---------------------------------------------------------------------------
# Phase mass vs. time (igasol.dat + sol_*.dat required)
# ---------------------------------------------------------------------------
if [[ -f "$RUN_DIR/igasol.dat" ]] && ls "$RUN_DIR"/sol_*.dat &>/dev/null 2>&1; then
    echo ""
    echo "--- Phase mass vs. time ---"
    set +e
    "$PYTHON" "$POSTPROCESS_DIR/plot_mass.py" \
        --dir "$RUN_DIR" --save "$RUN_DIR/mass.png" \
        2>&1 | sed 's/^/  /'
    (( overall_exit += $? )) || true
    set -e
fi

echo ""
echo "========================================================================="
if [[ "$overall_exit" -ne 0 ]]; then
    echo "  ⚠️  Post-processing completed with errors (exit sum $overall_exit)"
else
    echo "  ✅ Post-processing complete"
fi
echo "  Results in: $RUN_DIR"
echo "========================================================================="
echo ""

exit "$overall_exit"
