#!/usr/bin/env bash
# =============================================================================
# fetch_stage.sh — download one stage of the k_eff campaign, from your Mac.
#
#   ./scripts/HPC/fetch_stage.sh studies/keff_sintering/batch3b_T-20.txt
#   ./scripts/HPC/fetch_stage.sh <stage file> --full        # + snapshots, vtk, everything
#   ./scripts/HPC/fetch_stage.sh <stage file> --full-run <name-substring>
#   ./scripts/HPC/fetch_stage.sh <stage file> --dry-run     # list, transfer nothing
#   ./scripts/HPC/fetch_stage.sh <stage A> <stage B> ... --full-run <s1> --full-run <s2>
#                       # several stages and several full runs, still ONE rsync
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

REMOTE_HOST="${REMOTE_HOST-hpc}"     # REMOTE_HOST= (empty) reads a local REMOTE, for testing
REMOTE="${REMOTE:-/resnick/groups/rubyfu/jbaglino/simulation_outputs/enceladus_DSM/keff_sintering_campaign}"
LOCAL="${LOCAL:-$HOME/SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign}"

stages=() ; full=0 ; dry=0 ; fullpats=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --full) full=1; shift ;;
        --full-run) fullpats+=("$2"); shift 2 ;;
        --dry-run|-n) dry=1; shift ;;
        -h|--help) sed -n '2,28p' "$0"; exit 0 ;;
        -*) echo "unknown option: $1" >&2; exit 1 ;;
        *) stages+=("$1"); shift ;;
    esac
done
(( ${#stages[@]} )) || { echo "❌ no stage file given" >&2; exit 1; }
for stage in "${stages[@]}"; do
    [[ -f "$stage" ]] || { echo "❌ stage file not found: $stage" >&2; exit 1; }
done

filt=$(mktemp)
trap 'rm -f "$filt"' EXIT
echo "+ /stages/" >> "$filt"
n=0 ; nfull=0 ; names=()
for stage in "${stages[@]}"; do
stage_name="$(basename "$stage" .txt)"
names+=("$stage_name")
# the stage's own records
echo "+ /stages/${stage_name}__*/***" >> "$filt"
while IFS= read -r line; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"; line="${line%"${line##*[![:space:]]}"}"
    [[ -z "$line" || "$line" != *:* ]] && continue
    IFS=':' read -r geom exp _ <<< "$line"
    run="${geom}__${exp}"
    n=$((n + 1))
    want_full=$full
    for fp in ${fullpats[@]+"${fullpats[@]}"}; do [[ "$run" == *"$fp"* ]] && want_full=1; done
    if (( want_full )); then
        echo "+ /${run}/***" >> "$filt"; nfull=$((nfull + 1))
    else
        echo "+ /${run}/" >> "$filt"
        for p in 'k_eff*.csv' SSA_evo.dat 'outp*.txt' igasol.dat '*.opts' phi_bounds.csv \
                 'cost_job*.txt' '*.o[0-9]*' '*.e[0-9]*'; do
            echo "+ /${run}/${p}" >> "$filt"
        done
    fi
done < "$stage"
done
echo "- *" >> "$filt"

echo "stages   : ${names[*]} ($n runs, $nfull with full snapshots)"
for fp in ${fullpats[@]+"${fullpats[@]}"}; do
    grep -q "$fp" "$filt" || echo "⚠  --full-run '$fp' matches no run in these stages" >&2
done
echo "from     : $REMOTE_HOST:$REMOTE/"
echo "to       : $LOCAL/"
mkdir -p "$LOCAL"
args=(-avhP --prune-empty-dirs --filter="merge $filt")
(( dry )) && args+=(-n)
rsync "${args[@]}" "${REMOTE_HOST:+$REMOTE_HOST:}$REMOTE/" "$LOCAL/"
