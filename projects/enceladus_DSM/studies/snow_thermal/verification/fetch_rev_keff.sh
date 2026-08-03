#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# fetch_rev_keff.sh — pull ONLY the k_eff time series off the cluster.
#
# The rev_sweep output is far too large to download: the 3 mm window alone is a
# 16.7 M-DOF field written every output step. But everything plot_rev_tensor.py
# needs is the per-run k_eff.csv (a few kB) plus the .opts files that record
# what was run. This script rsyncs exactly those and nothing else, so the whole
# sweep comes down in well under a megabyte.
#
# Usage:
#   ./fetch_rev_keff.sh [remote_base] [local_dest]
#
#   remote_base  default /resnick/scratch/jbaglino/enceladus_DSM/rev_sweep
#   local_dest   default <this dir>/raw_rev
#
# The directory layout is mirrored, so the local tree is
#   raw_rev/rev_L<L>mm_T-20_rev/<timestamp>_<exp>_<tag>_job<id>/k_eff.csv
# which is exactly what plot_rev_tensor.py --results raw_rev expects.
#
# SSH to login.hpc.caltech.edu is interactive (2FA), so run this yourself and
# answer the prompt; it cannot be run unattended.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REMOTE_HOST="${REMOTE_HOST:-hpc}"
REMOTE_BASE="${1:-/resnick/scratch/jbaglino/enceladus_DSM/rev_sweep}"
LOCAL_DEST="${2:-$HERE/raw_rev}"

mkdir -p "$LOCAL_DEST"

echo "=== fetching k_eff series only ==="
echo "  from : $REMOTE_HOST:$REMOTE_BASE"
echo "  to   : $LOCAL_DEST"
echo

# --include of the directory levels is REQUIRED before the file patterns:
# rsync prunes a directory the moment it matches an --exclude, so the trailing
# '--exclude=*' would stop the descent at the top level and nothing would ever
# be tested against the k_eff.csv rule.
rsync -avz --prune-empty-dirs \
      --include='*/' \
      --include='k_eff.csv' \
      --include='*.opts' \
      --include='metadata.json' \
      --exclude='*' \
      "$REMOTE_HOST:$REMOTE_BASE/" "$LOCAL_DEST/"

echo
n=$(find "$LOCAL_DEST" -name k_eff.csv | wc -l | tr -d ' ')
sz=$(du -sh "$LOCAL_DEST" | cut -f1)
echo "got $n k_eff.csv files, $sz total"
echo
echo "Now plot with:"
echo
echo "  venv_enceladus/bin/python3 studies/snow_thermal/verification/plot_rev_tensor.py \\"
echo "      --results $LOCAL_DEST"
