#!/usr/bin/env bash
# =============================================================================
# supersede_runs.sh — move the existing run folders of a stage file aside, so
# the runs can be resubmitted. Run on the HPC.
#
#   ./scripts/HPC/supersede_runs.sh <stage file>            # list what would move
#   ./scripts/HPC/supersede_runs.sh <stage file> --yes      # move them
#
# submit_keff_production.sh never writes over a run folder. To REDO a run, its
# folder goes to <campaign>/stages/superseded__<date>/<run>/ -- kept, not
# deleted -- and the name is free again. Nothing is moved while a job of that
# name is queued or running.
# =============================================================================
set -euo pipefail
CAMPAIGN_DIR="${CAMPAIGN_DIR:-/resnick/groups/rubyfu/jbaglino/simulation_outputs/enceladus_DSM/keff_sintering_campaign}"
stage="${1:-}"; go="${2:-}"
[[ -f "$stage" ]] || { echo "usage: $0 <stage file> [--yes]" >&2; exit 1; }
dest="$CAMPAIGN_DIR/stages/superseded__$(date +%Y-%m-%d)"
n=0
while IFS= read -r line; do
    line="${line%%#*}"; line="${line#"${line%%[![:space:]]*}"}"; line="${line%"${line##*[![:space:]]}"}"
    [[ -z "$line" || "$line" != *:* ]] && continue
    IFS=':' read -r geom exp _ <<< "$line"
    run="${geom}__${exp}"
    [[ -d "$CAMPAIGN_DIR/$run" ]] || { echo "  (no folder)  $run"; continue; }
    if command -v squeue >/dev/null 2>&1 && [[ -n "$(squeue -h -u "${USER:-$(id -un)}" -n "$run" -o %i 2>/dev/null)" ]]; then
        echo "  IN QUEUE, left alone: $run"; continue
    fi
    if [[ "$go" == "--yes" ]]; then
        mkdir -p "$dest"
        mv "$CAMPAIGN_DIR/$run" "$dest/$run"
        echo "  moved  $run"
    else
        echo "  would move  $run"
    fi
    n=$((n + 1))
done < "$stage"
if [[ "$go" == "--yes" ]]; then echo "$n folder(s) moved to $dest"
else echo "$n folder(s) would move to $dest -- rerun with --yes"; fi
