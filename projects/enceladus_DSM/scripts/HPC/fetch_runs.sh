#!/usr/bin/env bash
# =============================================================================
# fetch_runs.sh — pull runs from the cluster in two stages.
#
#   ./scripts/HPC/fetch_runs.sh --tables <remote-glob> <local-dir>
#   ./scripts/HPC/fetch_runs.sh --fields <files.txt> <remote-root> <local-dir>
#
# Uses the `hpc:` host alias from ~/.ssh/config, same as hpcget.
#
# WHY NOT hpcget. Two reasons, both specific rather than stylistic:
#
#   * hpcget passes --ignore-existing, which SKIPS any file already present
#     locally. That is right for bulk archival pulls and wrong here: a table
#     re-pulled after a rerun would be silently kept at its old contents. This
#     script omits it, so rsync's normal size/mtime check decides.
#   * hpcget takes no filters. At L/R_ave = 64 a snapshot is ~460 MB and a run
#     writes one per step, so an unfiltered pull is ~170 GB per run. The point
#     of staging is to never ask for those.
#
# STAGE 1 (--tables) is a few MB per run and is usually all that is needed:
# k_eff.csv, SSA_evo.dat, outp.txt, igasol.dat, *.opts. The REV comparison runs
# entirely off these, because both ensembles carry the same k_eff bias and the
# comparison is therefore valid without re-measuring either.
#
# STAGE 2 (--fields) takes the list from postprocess/plan_download.py, which
# picks one snapshot per k_eff sample -- ~55 files instead of ~370, ~2% of the
# bytes, with every k_eff row still backed by a field. Only needed to recompute
# k_eff with the sharp or tensorial coefficient.
#
# Structure is preserved per run: each remote run directory lands in its own
# local subdirectory, so the analysis tools (which discover runs by content)
# find them.
# =============================================================================
set -euo pipefail

usage() { sed -n '2,32p' "${BASH_SOURCE[0]}"; exit "${1:-0}"; }
[ $# -eq 0 ] && usage 1

TABLES=(--include='*/' --include='k_eff.csv' --include='SSA_evo.dat'
        --include='outp.txt' --include='igasol.dat' --include='*.opts'
        --exclude='*')

case "${1:-}" in
  --tables)
    shift
    [ $# -eq 2 ] || usage 1
    remote="$1"; local_dir="$2"
    mkdir -p "$local_dir"
    # Expand the remote glob ON THE CLUSTER, then pull each match into its own
    # subdirectory. A single rsync with several sources would merge them all
    # into one directory and lose which run is which.
    # No mapfile: macOS ships bash 3.2, which does not have it. (Shipped that
    # bug once already in resume_batch.sh.) A here-string keeps the loop in the
    # current shell so set -e still applies to a failed rsync.
    listing=$(ssh hpc "ls -d ${remote} 2>/dev/null" || true)
    [ -z "$listing" ] && { echo "no remote match for: $remote" >&2; exit 1; }
    echo "=== $(printf '%s\n' "$listing" | wc -l | tr -d ' ') run(s) -> $local_dir ==="
    while IFS= read -r d; do
      [ -n "$d" ] || continue
      name=$(basename "$d")
      echo "  $name"
      rsync -avhP "${TABLES[@]}" "hpc:$d/" "$local_dir/$name/"
    done <<< "$listing"
    echo ""
    echo "Next: venv_enceladus/bin/python postprocess/plan_download.py $local_dir --out files.txt"
    ;;
  --fields)
    shift
    [ $# -eq 3 ] || usage 1
    list="$1"; remote_root="$2"; local_dir="$3"
    [ -f "$list" ] || { echo "no such list: $list" >&2; exit 1; }
    n=$(wc -l < "$list" | tr -d ' ')
    echo "=== $n snapshot(s) from $remote_root -> $local_dir ==="
    rsync -avhP --files-from="$list" "hpc:$remote_root/" "$local_dir/"
    ;;
  -h|--help) usage 0 ;;
  *) usage 1 ;;
esac
