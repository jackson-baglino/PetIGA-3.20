#!/usr/bin/env bash
# =============================================================================
# fetch_stage.sh — download one stage of the k_eff campaign, from your Mac.
#
#   ./scripts/HPC/fetch_stage.sh studies/keff_sintering/batch3b_T-20.txt
#   ./scripts/HPC/fetch_stage.sh <stage file> --full        # + snapshots, vtk, everything
#   ./scripts/HPC/fetch_stage.sh <stage file> --full-run <name-substring>
#   ./scripts/HPC/fetch_stage.sh <stage file> --dry-run     # list, transfer nothing
#
# Every manuscript run lives in ONE folder on the cluster,
#   $REMOTE = .../simulation_outputs/enceladus_DSM/keff_sintering_campaign/<geom>__<exp>/
# (scripts/HPC/submit_keff_production.sh), and the stage file in the repo says
# which runs a stage holds. So the include list is built HERE, from the stage
# file, and the whole stage comes down in ONE rsync -- one 2FA prompt -- into
# the same layout locally:
#   $LOCAL = ~/SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign/
# Re-running is cheap: rsync only moves what changed.
#
# Default is the TABLES: k_eff*.csv, SSA_evo.dat, outp*.txt, igasol.dat, *.opts,
# phi_bounds.csv, cost_job*.txt, the SLURM .o/.e, plus the stage's records in
# stages/<stage>__*/. That is everything the k_eff analysis, the health check
# and the cost check read (a few MB per run). --full adds the snapshots
# (~9 GB/run at L/R 40); --full-run does that for matching runs only.
# =============================================================================
set -euo pipefail

REMOTE_HOST="${REMOTE_HOST:-hpc}"
REMOTE="${REMOTE:-/resnick/groups/rubyfu/jbaglino/simulation_outputs/enceladus_DSM/keff_sintering_campaign}"
LOCAL="${LOCAL:-$HOME/SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign}"

stage="" ; full=0 ; dry=0 ; fullpat=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --full) full=1; shift ;;
        --full-run) fullpat="$2"; shift 2 ;;
        --dry-run|-n) dry=1; shift ;;
        -h|--help) sed -n '2,26p' "$0"; exit 0 ;;
        -*) echo "unknown option: $1" >&2; exit 1 ;;
        *) stage="$1"; shift ;;
    esac
done
[[ -f "$stage" ]] || { echo "❌ stage file not found: ${stage:-<none>}" >&2; exit 1; }
stage_name="$(basename "$stage" .txt)"

filt=$(mktemp)
trap 'rm -f "$filt"' EXIT
# the stage's own records
printf '%s\n' "+ /stages/" "+ /stages/${stage_name}__*/***" >> "$filt"
n=0
while IFS= read -r line; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"; line="${line%"${line##*[![:space:]]}"}"
    [[ -z "$line" || "$line" != *:* ]] && continue
    IFS=':' read -r geom exp _ <<< "$line"
    run="${geom}__${exp}"
    n=$((n + 1))
    if (( full )) || [[ -n "$fullpat" && "$run" == *"$fullpat"* ]]; then
        echo "+ /${run}/***" >> "$filt"
    else
        echo "+ /${run}/" >> "$filt"
        for p in 'k_eff*.csv' SSA_evo.dat 'outp*.txt' igasol.dat '*.opts' phi_bounds.csv \
                 'cost_job*.txt' '*.o[0-9]*' '*.e[0-9]*'; do
            echo "+ /${run}/${p}" >> "$filt"
        done
    fi
done < "$stage"
echo "- *" >> "$filt"

echo "stage    : $stage_name ($n runs)$( ((full)) && echo ', FULL' )${fullpat:+, full for *$fullpat*}"
echo "from     : $REMOTE_HOST:$REMOTE/"
echo "to       : $LOCAL/"
mkdir -p "$LOCAL"
args=(-avhP --prune-empty-dirs --filter="merge $filt")
(( dry )) && args+=(-n)
rsync "${args[@]}" "$REMOTE_HOST:$REMOTE/" "$LOCAL/"
