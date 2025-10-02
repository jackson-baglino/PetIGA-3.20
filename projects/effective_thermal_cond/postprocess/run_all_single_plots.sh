#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Script: run_all_single_plots.sh
# Purpose: Loop over all DSM result folders in a parent directory and
#          call plot_single_run_k_eff.py for each one.
# Usage:   ./run_all_single_plots.sh /path/to/parent_dir
# ----------------------------------------------------------------------

set -euo pipefail

PARENT_DIR="${1:-.}"   # default to current dir if not provided
SCRIPT_DIR="$(dirname "$0")"
PLOT_SCRIPT="${SCRIPT_DIR}/plot_single_run_k_eff.py"

# The plot script doesn't need the executable bit; we'll run it via python3.
# Check for existence (-f) instead of executability (-x).
if [[ ! -f "$PLOT_SCRIPT" ]]; then
  echo "[ERROR] Cannot find plot_single_run_k_eff.py at $PLOT_SCRIPT"
  echo "        (SCRIPT_DIR resolved to: $SCRIPT_DIR)"
  exit 1
fi

echo "[INFO] Scanning parent directory: $PARENT_DIR"
echo "[INFO] Using plot script: $PLOT_SCRIPT"

# loop over subdirectories
for folder in "$PARENT_DIR"/*/; do
  # skip if not a directory
  [[ -d "$folder" ]] || continue

  echo "=== Processing $(basename "$folder") ==="
  python3 "$PLOT_SCRIPT" "$folder" --theme both || {
    echo "[WARN] Plotting failed for $folder, skipping."
  }
done

echo "[INFO] All done."