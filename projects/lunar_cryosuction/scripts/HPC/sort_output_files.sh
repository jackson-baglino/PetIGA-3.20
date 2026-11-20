#!/usr/bin/env bash
# sort_output_files.sh: move the old flat SLURM logs (output_files/<geom>__<exp>.o<id> / .e<id>)
# into output_files/<batch folder name>/, the layout submit_batch.sh uses since 2026-11-20.
# A log goes to the batch whose folder holds a run folder of the same name. Nothing is deleted;
# logs with no matching batch, or with more than one, stay where they are and are listed.
#
#   ./scripts/HPC/sort_output_files.sh [--apply] [output_files dir] [results root]
#
# Default is a dry run (prints the moves). Defaults: output_files/ and
# /resnick/groups/rubyfu/jbaglino/simulation_outputs/lunar_cryosuction. Run from the project root.
set -euo pipefail
apply=0
[[ "${1:-}" == "--apply" ]] && { apply=1; shift; }
of=${1:-output_files}
root=${2:-/resnick/groups/rubyfu/jbaglino/simulation_outputs/lunar_cryosuction}
[[ -d "$of" && -d "$root" ]] || { echo "❌ Need an output_files dir ($of) and a results root ($root)"; exit 1; }

# run folder name -> batch name(s), from every batch_* folder under the root
map=$(find "$root" -mindepth 2 -maxdepth 2 -type d -path '*/batch_*/*__*' \
        ! -name 'inputs_snapshot*' ! -name 'src_snapshot*' | awk -F/ '{print $NF "\t" $(NF-1)}')

n_move=0; n_none=0; n_many=0
while IFS= read -r f; do
  base=$(basename "$f"); stem=${base%.[oe][0-9]*}
  batches=$(awk -F'\t' -v s="$stem" '$1==s {print $2}' <<<"$map" | sort -u)
  nb=$(grep -c . <<<"$batches" || true)
  if (( nb == 1 )); then
    echo "move  $base -> $batches/"
    if (( apply )); then mkdir -p "$of/$batches"; mv -n "$f" "$of/$batches/"; fi
    n_move=$((n_move + 1))
  elif (( nb == 0 )); then
    echo "keep  $base  (no batch has a run folder of this name)"; n_none=$((n_none + 1))
  else
    echo "keep  $base  (in several batches: $(tr '\n' ' ' <<<"$batches"))"; n_many=$((n_many + 1))
  fi
done < <(find "$of" -maxdepth 1 -type f -name '*__*.[oe][0-9]*' | sort)

echo "---"
echo "$( (( apply )) && echo moved || echo 'would move'): $n_move   kept, no batch: $n_none   kept, ambiguous: $n_many"
(( apply )) || echo "(dry run; add --apply to move)"
