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
# (anything that contains an `igasol.dat`) and runs:
#   - postprocess/plot_mass.py             (always; produces mass.png +
#                                            mass_plots/{total,ice,sediment,
#                                            vapor,change_loglog}.png)
#   - postprocess/plotpermafrost.py        (dim >= 2: VTK conversion)
#   - postprocess/plot1D_profiles.py       (dim == 1: phase-field profiles)
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

    # ── 1. plot_mass.py (always) ────────────────────────────────────────
    if [[ -f "$POSTPROCESS/plot_mass.py" ]]; then
        echo "  plot_mass.py ..."
        if ! "$PYTHON" "$POSTPROCESS/plot_mass.py" --dir "$run" 2>&1 | sed 's/^/    /'; then
            echo "  ⚠ plot_mass.py failed"
            failed_this=1
        fi
    fi

    # ── 2. dim-specific field plots ─────────────────────────────────────
    if [[ "$dim" == "1" ]]; then
        if [[ -f "$POSTPROCESS/plot1D_profiles.py" ]]; then
            echo "  plot1D_profiles.py ..."
            if ! "$PYTHON" "$POSTPROCESS/plot1D_profiles.py" --dir "$run" 2>&1 | sed 's/^/    /'; then
                echo "  ⚠ plot1D_profiles.py failed"
                failed_this=1
            fi
        fi
    else
        if [[ -f "$POSTPROCESS/plotpermafrost.py" ]]; then
            echo "  plotpermafrost.py (VTK) ..."
            if ! "$PYTHON" "$POSTPROCESS/plotpermafrost.py" --dir "$run" 2>&1 | sed 's/^/    /'; then
                echo "  ⚠ plotpermafrost.py failed"
                failed_this=1
            fi
        fi
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
