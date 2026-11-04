#!/usr/bin/env bash
# =============================================================================
# run_batch_postprocess.sh — local post-processing for a downloaded HPC batch.
#
# Submitted by scripts/HPC/submit_batch.sh and copied into the batch parent
# folder, so the workflow is:
#
#   1. On HPC: ./scripts/HPC/submit_batch.sh --tag mytag --tests "..."
#   2. Wait for jobs to finish.
#   3. Locally: rsync the whole batch parent down.
#   4. Locally: bash <downloaded-parent>/run_batch_postprocess.sh
#
# The script auto-detects every per-test subfolder of the parent directory
# (anything that contains an `igasol.dat`) and runs the standard sweep into
# that run's own `plots/` folder:
#   - plot_porosity.py   porosity + interface density   -> plots/porosity.png
#   - plot_ssa.py        ice-air surface area           -> plots/ssa.png
#   - plot_mass.py       phase mass vs time             -> plots/mass{.png,/}
#   - plot_timestep.py   adaptive dt history            -> plots/timestep.png
#   - plot_fields.py     VTK conversion (dim >= 2)      -> vtkOut/
#   - wedge_gt_velocity.py  Gibbs-Thomson v_n check (wedge geometries only)
#                                                      -> plots/wedge_gt_velocity.png
#
# plot_fields.py is a format conversion rather than a figure, so it stays in
# vtkOut/. plot_fields_highres.py is not in the sweep: it is slow and only
# wanted on demand.
#
# 1D runs get no field plots -- plot_1d_profiles.py is in postprocess/scratch/.
#
# It uses the `postprocess/` directory that submit_batch.sh staged at the
# batch root, so no source tree is required.
#
# Usage:
#   bash run_batch_postprocess.sh                      # from the batch folder
#   bash run_batch_postprocess.sh /path/to/batch       # from anywhere
# =============================================================================
set -uo pipefail

if [[ $# -ge 1 ]]; then
    BATCH_DIR="$(cd "$1" && pwd)"
else
    BATCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Locate the postprocess/ directory. Prefer the one staged inside the batch
# (so the user can run this without the project repo). Fall back to a
# postprocess/ in the parent of run_batch_postprocess.sh if the staged copy
# is missing.
if   [[ -d "$BATCH_DIR/postprocess" ]]; then
    POSTPROCESS="$BATCH_DIR/postprocess"
elif [[ -d "$(dirname "$BATCH_DIR")/postprocess" ]]; then
    POSTPROCESS="$(dirname "$BATCH_DIR")/postprocess"
else
    echo "❌ Could not find postprocess/ next to or inside $BATCH_DIR"
    echo "   submit_batch.sh stages it automatically; re-download or re-run."
    exit 1
fi

PYTHON="$(command -v python3 || command -v python || true)"
if [[ -z "$PYTHON" ]]; then
    echo "❌ Neither python3 nor python found on PATH."
    exit 1
fi

echo "============================================================"
echo "  Batch post-processing"
echo "  Batch dir   : $BATCH_DIR"
echo "  postprocess : $POSTPROCESS"
echo "  python      : $PYTHON"
echo "============================================================"

# Iterate over every per-test subdir (one containing igasol.dat = a real run).
# Sort for deterministic order.
n_total=0
n_done=0
n_skipped=0
n_failed=0

for run in "$BATCH_DIR"/*/; do
    run="${run%/}"
    name="$(basename "$run")"

    # Skip the staging directories created by submit_batch.sh
    case "$name" in
        inputs_snapshot|src_snapshot|postprocess) continue ;;
    esac

    if [[ ! -f "$run/igasol.dat" ]]; then
        echo ""
        echo "⚠ $name — no igasol.dat, skipping"
        ((n_skipped++))
        continue
    fi

    ((n_total++))
    echo ""
    echo "------------------------------------------------------------"
    echo "▶ $name"
    echo "------------------------------------------------------------"

    # Pull dim from the staged geometry opts file inside the run
    dim=$(awk '$1=="-dim"{print $2; exit}' "$run"/*.opts 2>/dev/null | head -n1)
    dim=${dim:-2}

    failed_this=0

    # Every figure for this run goes in its own plots/ subfolder.
    plots="$run/plots"
    mkdir -p "$plots"

    # Run one script, folding its REAL status in. `cmd | sed` then `$?` reads
    # sed's status (always 0), so the previous `if ! cmd | sed` never fired.
    step() {
        local script="$POSTPROCESS/$1"; shift
        [[ -f "$script" ]] || return 0
        echo "  $(basename "$script") ..."
        "$PYTHON" "$script" "$@" 2>&1 | sed 's/^/    /'
        local rc=${PIPESTATUS[0]}
        (( rc != 0 )) && { echo "  ⚠ $(basename "$script") exited $rc"; failed_this=1; }
        return 0
    }

    # ── 1. scalar time-series from SSA_evo.dat ──────────────────────────
    if [[ -f "$run/SSA_evo.dat" ]]; then
        step plot_porosity.py --dir "$run" --save "$plots/porosity.png"
        step plot_ssa.py      --dir "$run" --save "$plots/ssa.png"
    fi

    # ── 2. phase mass ───────────────────────────────────────────────────
    step plot_mass.py --dir "$run" --save "$plots/mass.png" \
                      --per-phase-dir "$plots/mass"

    # ── 3. time step diagnostic ─────────────────────────────────────────
    if [[ -f "$run/outp.txt" ]]; then
        step plot_timestep.py --dir "$run" --save "$plots/timestep.png"
    fi

    # ── 4. Gibbs-Thomson check on the wedge menisci. Only meaningful for
    #      the apex-centred band geometry, so it is gated on the apex being
    #      declared; every other geometry falls through silently. Needs the
    #      .vts, hence after the conversion below reruns -- on a first pass
    #      over a fresh batch it is skipped and picked up on the second. ──
    if grep -qs '^-wedge_apex_x' "$run"/*.opts && \
       compgen -G "$run/vtkOut/solV_*.vts" >/dev/null; then
        step wedge_gt_velocity.py --dir "$run" \
             --save "$plots/wedge_gt_velocity.png"
    fi

    # ── 5. VTK conversion (2D/3D only). A format conversion, not a figure,
    #      so it stays in vtkOut/ rather than plots/. ────────────────────
    if [[ "$dim" != "1" ]]; then
        step plot_fields.py --dir "$run"
    else
        echo "  (dim 1 — no 1D plotter in the current sweep; see"
        echo "   postprocess/scratch/plot_1d_profiles.py)"
    fi

    if (( failed_this )); then
        ((n_failed++))
    else
        ((n_done++))
        echo "  ✅ $name done"
    fi
done

echo ""
echo "============================================================"
echo "  Batch post-processing summary"
echo "  Total runs : $n_total"
echo "  Done       : $n_done"
echo "  Failed     : $n_failed"
echo "  Skipped    : $n_skipped (no igasol.dat)"
echo "============================================================"
