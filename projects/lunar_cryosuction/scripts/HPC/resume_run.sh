#!/bin/bash
# resume_run.sh: continue a stopped run (time limit, node failure) in its own
# folder, from its last sol_N.dat, as one new SLURM job.
#
# Usage (from the project root on the HPC, after `git pull`):
#   ./scripts/HPC/resume_run.sh <run_folder> [--t-final-days D] [--dry-run] [-- <sbatch options>]
# e.g.
#   ./scripts/HPC/resume_run.sh \
#     /resnick/groups/rubyfu/jbaglino/simulation_outputs/lunar_cryosuction/batch_X/geom__exp \
#     -- --time=2-00:00:00
#
# The run folder must be a batch run folder (<geom>__<exp>, as written by
# submit_batch.sh) or a single-run folder holding exactly one geometry and one
# experiment .opts besides solver.opts. The job runs scripts/HPC/run_cryosuction.sh
# with RESUME_DIR set; see "Resume mode" there for what it does to the folder
# (log backup and trim, -initial_cond ... -restart_append, outp.txt appended).
#
# --t-final-days D extends a run (also one that already finished) to D days:
# the experiment .opts in the run folder is copied to <exp>.opts.before_extend_<ts>
# and gets -t_final and -n_out lines appended (later lines win, in PETSc and in
# pplib.read_opts). n_out is scaled so the snapshot interval t_final/(n_out-1)
# stays the same, so the second leg continues the first leg's output grid.
# --dry-run prints what would be done and stops before changing anything,
# building or submitting.
#
# The job gets the same name and the same rank count as the first leg
# (TARGET_DOFS_PER_CORE as given, default from scripts/lib/alloc.sh); a SLURM
# job ID cannot be reused, so the new leg has a new one.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/run_cryosuction.sh"
source "$PROJECT_ROOT/scripts/lib/alloc.sh"

[[ $# -ge 1 ]] || { sed -n '2,28p' "$0"; exit 1; }
run_dir="$(cd "$1" && pwd)"; shift
sbatch_extra=(); t_final_days=""; dry=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --t-final-days) [[ $# -ge 2 ]] || { echo "❌ --t-final-days needs a value"; exit 1; }
                        t_final_days="$2"; shift 2 ;;
        --dry-run) dry=1; shift ;;
        --) shift; sbatch_extra=("$@"); break ;;
        *) echo "❌ Unknown argument: $1"; exit 1 ;;
    esac
done

# --- identify geometry and experiment ------------------------------------
base="$(basename "$run_dir")"
if [[ "$base" == *__* ]]; then
    geom="${base%%__*}"; exp="${base#*__}"
else
    mapfile -t others < <(find "$run_dir" -maxdepth 1 -name '*.opts' ! -name solver.opts -exec basename {} .opts \;)
    [[ ${#others[@]} -eq 2 ]] || { echo "❌ Expected 2 non-solver .opts in $run_dir, found ${#others[@]}"; exit 1; }
    if grep -q '^-geom_file\|^-Nx' "$run_dir/${others[0]}.opts"; then geom="${others[0]}"; exp="${others[1]}"
    else geom="${others[1]}"; exp="${others[0]}"; fi
fi
for f in solver.opts "$geom.opts" "$exp.opts" SSA_evo.dat; do
    [[ -f "$run_dir/$f" ]] || { echo "❌ Missing $run_dir/$f"; exit 1; }
done
# Last COMPLETE snapshot (same size as sol_00000.dat); the job itself sets
# shorter ones aside (run_cryosuction.sh, prepare_resume).
full="$(stat -c %s "$run_dir/sol_00000.dat" 2>/dev/null || stat -f %z "$run_dir/sol_00000.dat")"
last=""
while IFS= read -r f; do
    sz="$(stat -c %s "$f" 2>/dev/null || stat -f %z "$f")"
    if [[ "$sz" == "$full" ]]; then last="$f"; else echo "⚠️  $(basename "$f") is incomplete ($sz of $full bytes); it will be set aside"; fi
done < <(find "$run_dir" -maxdepth 1 -name 'sol_*.dat' | sort)
[[ -n "$last" ]] || { echo "❌ No complete sol_*.dat in $run_dir"; exit 1; }

# Already finished? Compare the last logged time with -t_final.
t_final="$(awk '$1=="-t_final"{print $2}' "$run_dir/$exp.opts" | tail -n1)"
t_last="$(awk '!/^#/ {t=$3} END {print t}' "$run_dir/SSA_evo.dat")"
echo "Run folder : $run_dir"
echo "Geometry   : $geom"
echo "Experiment : $exp"
echo "Last snap  : $(basename "$last")"
echo "Logged t   : $t_last s of t_final $t_final s"
if [[ -n "$t_final_days" ]]; then
    n_out="$(awk '$1=="-n_out"{print $2}' "$run_dir/$exp.opts" | tail -n1)"
    [[ -n "$n_out" ]] || { echo "❌ No -n_out in $exp.opts; cannot keep the snapshot interval"; exit 1; }
    new_t="$(awk -v d="$t_final_days" 'BEGIN {printf "%.6e", d * 86400}')"
    awk -v a="$new_t" -v b="$t_last" 'BEGIN {exit !(a > b)}' \
        || { echo "❌ --t-final-days $t_final_days ($new_t s) is not past the logged t $t_last s"; exit 1; }
    new_n="$(awk -v n="$n_out" -v a="$new_t" -v b="$t_final" 'BEGIN {printf "%d", (n - 1) * a / b + 1.5}')"
    echo "Extend     : t_final $t_final -> $new_t s ($t_final_days d), n_out $n_out -> $new_n (same snapshot interval)"
    if [[ $dry -eq 0 ]]; then
        ts="$(date +%Y-%m-%d__%H.%M.%S)"
        cp -p "$run_dir/$exp.opts" "$run_dir/$exp.opts.before_extend_$ts"
        {
            echo ""
            echo "# Extended by resume_run.sh $ts: t_final $t_final -> $new_t s ($t_final_days d),"
            echo "# n_out $n_out -> $new_n (same snapshot interval). Original: $exp.opts.before_extend_$ts"
            echo "-t_final $new_t"
            echo "-n_out $new_n"
        } >> "$run_dir/$exp.opts"
    fi
    t_final="$new_t"
fi
if awk -v a="$t_last" -v b="$t_final" 'BEGIN {exit !(a >= b)}'; then
    echo "❌ The run already reached t_final; nothing to resume."; exit 1
fi

# --- allocation: same rule as submit_batch.sh ------------------------------
dof=$(awk '$1=="-dof"{print $2}' "$run_dir/solver.opts" | head -n1); dof=${dof:-3}
if grep -q "^-geom_file" "$run_dir/$geom.opts"; then
    read -r nx ny nz <<< "$(awk '$1=="#" && $2=="DOF_GRID:"{print $3, $4, $5}' "$run_dir/$geom.opts" | head -n1)"
    gf="$(awk '$1=="-geom_file"{print $2}' "$run_dir/$geom.opts" | head -n1)"
    [[ -f "$PROJECT_ROOT/$gf" ]] || { echo "❌ Mesh file not found: $PROJECT_ROOT/$gf"; exit 1; }
else
    nx=$(awk '$1=="-Nx"{print $2}' "$run_dir/$geom.opts" | head -n1)
    ny=$(awk '$1=="-Ny"{print $2}' "$run_dir/$geom.opts" | head -n1)
    nz=$(awk '$1=="-Nz"{print $2}' "$run_dir/$geom.opts" | head -n1)
fi
nx=${nx:-1}; ny=${ny:-1}; nz=${nz:-1}
total_dofs=$((dof * nx * ny * nz))
nprocs=$(( (total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE ))
(( nprocs < 1 )) && nprocs=1
tpn=$nprocs; (( tpn > MAX_TASKS_PER_NODE )) && tpn=$MAX_TASKS_PER_NODE
nnodes=$(( (nprocs + tpn - 1) / tpn ))
echo "Allocation : DoFs=$total_dofs nprocs=$nprocs nodes=$nnodes tasks/node=$tpn"

# --- build on the submission host, then submit ------------------------------
[[ $dry -eq 0 ]] || { echo "(dry run: nothing changed, built or submitted)"; exit 0; }
cd "$PROJECT_ROOT"
make all || { echo "❌ Build failed."; exit 1; }
sbatch --job-name="${geom}__${exp}" \
       --nodes="$nnodes" \
       --ntasks="$nprocs" \
       --ntasks-per-node="$tpn" \
       --mem-per-cpu="$(mem_per_cpu "$total_dofs" "$nprocs")" \
       --export=ALL,SKIP_COMPILE=1,RESUME_DIR="$run_dir" \
       "${sbatch_extra[@]}" \
       "$RUN_SCRIPT" "$geom" "$exp" "resume"
