#!/usr/bin/env bash
# =============================================================================
# resume_batch.sh — restart every run in a killed batch from its last good
# snapshot, without typing any of it in.
#
# WHY THIS EXISTS. A resume needs, per run, the geometry name, the experiment
# name and the snapshot path. Those are long, there is one set per seed, and a
# wrong PAIRING is silent: IGAReadVec checks the vector LENGTH, so it catches a
# wrong mesh, but a snapshot from seed 2 loaded into seed 3's geometry has the
# same length and loads happily. Everything needed is already in the batch
# directory layout,
#
#     $SCRATCH/enceladus_DSM/batch_<ts>[_<tag>]/<geom>__<exp>/sol_NNNNN.dat
#
# so this derives all three from it and never pairs them by hand.
#
# THE SECOND-TO-LAST SNAPSHOT, NOT THE LAST. SSA_evo.dat is flushed every step
# so it is always intact, but a snapshot the scheduler interrupted mid-write
# can be truncated. One extra step is cheap insurance. --last overrides if you
# have checked the final file yourself.
#
# The clock and the time step are NOT passed: the solver reads both from the
# snapshot's own SSA_evo.dat (see RestartStateFromLog). Anything you do pass
# via --extra-opts still wins.
#
# USAGE
#   ./scripts/HPC/resume_batch.sh <batch_dir> [--last] [--dry-run] [--tag T]
#                                 [--extra-opts "..."] [-- <sbatch flags>]
#
# Submits through submit_enceladus.sh, which sbatch's each job. It does NOT
# call run_enceladus.sh: that is the script sbatch executes, and running it
# directly puts the solver on the login node with no allocation.
#
# Defaults to the same k_eff sampling the pilot used.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

KEFF_OPTS="-keff 1 -keff_step0 1 -keff_t_interv 5.0e4 -keff_ksp_type cg -keff_pc_type gamg"

batch=""; use_last=0; dry=0; extra=""; tag="resume"; sbatch_extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        --last)        use_last=1; shift ;;
        --dry-run)     dry=1; shift ;;
        --extra-opts)  extra="$2"; shift 2 ;;
        --tag)         tag="$2"; shift 2 ;;
        --)            shift; sbatch_extra=("$@"); break ;;
        -h|--help)     sed -n '2,32p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)             batch="$1"; shift ;;
    esac
done
[ -n "$batch" ] && [ -d "$batch" ] || { echo "usage: $0 <batch_dir> [--last] [--dry-run]" >&2; exit 2; }

# Refuse to "submit" on a machine with no scheduler: without this the loop
# would happily fall through to running jobs wherever it was invoked.
if [ "$dry" -eq 0 ] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not found -- this submits SLURM jobs and must run on a" >&2
    echo "       cluster login node. Use --dry-run to inspect the commands." >&2
    exit 1
fi

echo "batch: $batch"
echo "tag:   $tag"
echo ""
n=0
for run in "$batch"/*__*/; do
    [ -d "$run" ] || continue
    run="${run%/}"                     # no trailing slash: it doubles up in -initial_cond
    name=$(basename "$run")
    geom="${name%%__*}"
    exp="${name#*__}"

    # find, not ls: under `set -o pipefail` a glob that matches nothing makes
    # ls exit non-zero and takes the whole script down, so a run that died
    # before its first snapshot would abort the resume of every later one.
    # Portable to bash 3.2 as well (no mapfile, no negative indexing).
    nsnap=$(find "$run" -maxdepth 1 -name 'sol_*.dat' | wc -l | tr -d ' ')
    if [ "$nsnap" -eq 0 ]; then
        echo "  SKIP $name — no sol_*.dat written yet"; continue
    fi
    if [ "$use_last" -eq 1 ] || [ "$nsnap" -eq 1 ]; then
        snap=$(find "$run" -maxdepth 1 -name 'sol_*.dat' | sort | tail -1)
    else
        snap=$(find "$run" -maxdepth 1 -name 'sol_*.dat' | sort | tail -2 | head -1)
    fi
    step=$(basename "$snap" .dat | sed 's/sol_0*//')
    # Report what the solver will recover, so a bad pick is visible BEFORE
    # burning an allocation on it.
    tdt=$(awk -v s="${step:-0}" '$4==s {t=$3; d=$5} END{print (t==""?"?":t), (d==""?"?":d)}' \
          "$run/SSA_evo.dat" 2>/dev/null || echo "? ?")
    t=${tdt%% *}; dt=${tdt##* }
    echo "  $name"
    echo "      snapshot  $(basename "$snap")  of $nsnap"
    echo "      resumes   t = ${t:-?} s   dt = ${dt:-?} s"

    # submit_enceladus.sh, NOT run_enceladus.sh. run_enceladus.sh is the
    # script sbatch EXECUTES; invoking it directly runs the solver on whatever
    # node you are sitting on -- the login node -- with no allocation. The two
    # names are one word apart and the failure is silent until someone notices
    # the head node is pinned. submit_enceladus.sh sizes the rank count and
    # sbatch's it.
    #
    # Argument order matters: <geom> <exp> [tag] [sbatch flags] -- [solver flags]
    cmd=(./scripts/HPC/submit_enceladus.sh "$geom" "$exp" "$tag")
    [ "${#sbatch_extra[@]}" -gt 0 ] && cmd+=("${sbatch_extra[@]}")
    cmd+=(-- $KEFF_OPTS -initial_cond "$snap" $extra)

    if [ "$dry" -eq 1 ]; then
        printf '      $ '; printf '%q ' "${cmd[@]}"; echo
    else
        "${cmd[@]}"
    fi
    n=$((n+1))
    echo ""
done
echo "$n run(s) $( [ "$dry" -eq 1 ] && echo 'would be resumed (dry run)' || echo resumed )"
