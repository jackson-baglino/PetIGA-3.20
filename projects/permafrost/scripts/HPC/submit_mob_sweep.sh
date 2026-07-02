#!/usr/bin/env bash
# =============================================================================
# submit_mob_sweep.sh — Fan out one SLURM job per mob_sub value.
#
# All jobs share a single timestamped parent folder; each gets its own
# mob_<value>/ subfolder within it. Compiles once on the login node, then
# submits inline batch scripts that skip recompilation.
#
# Usage (from project root):
#   ./scripts/HPC/submit_mob_sweep.sh [experiment] [tag] [-- sbatch_overrides...]
#
#   experiment   Experiment opts name (default: 30day_T-5_h1.00_GTphys)
#   tag          Optional label appended to the batch folder name
#
# Examples:
#   ./scripts/HPC/submit_mob_sweep.sh
#   ./scripts/HPC/submit_mob_sweep.sh 5yr_T-5_h1.00_GTphys longrun
#   ./scripts/HPC/submit_mob_sweep.sh 30day_T-5_h1.00_GTphys -- --time=0-06:00:00
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$SCRIPT_DIR/../.." && pwd)"
EXEC="$PROJ/permafrost"
SOLVER_OPTS="$PROJ/inputs/solver.opts"

GEOM="2D_two_ice_grains_boundary"
GEOM_OPTS="$PROJ/inputs/geometry/${GEOM}.opts"

EXP="${1:-30day_T-5_h1.00_GTphys}"; shift || true
TAG="${1:-}"; [[ "${1:-}" == "--" ]] && TAG="" || shift || true

sbatch_extra=()
if [[ "${1:-}" == "--" ]]; then shift; sbatch_extra=("$@"); fi

EXP_OPTS="$PROJ/inputs/experiment/${EXP}.opts"

# ---------------------------------------------------------------------------
# mob_sub values to sweep (baseline → 100× reduction)
#
# "kp_default" is a special sentinel: the geometry opts -mob_sub override is
# stripped so the executable computes mob_sub via K&P matching at the actual
# eps in the opts file. This is the physically-correct but untuned reference.
# ---------------------------------------------------------------------------
MOB_VALUES=(
    kp_default   # K&P-computed at current eps (untuned reference)
    6.5404e-10   # frozen at safety=0.5 K&P values; baseline
    1.3e-10      #  5× reduction relative to safety=0.5 baseline
    6.5e-11      # 10× reduction
    1.3e-11      # 50× reduction
    6.5e-12      # 100× reduction
)

# ---------------------------------------------------------------------------
# Resource sizing from geometry opts
# ---------------------------------------------------------------------------
Nx=$(awk '$1=="-Nx"{print $2}' "$GEOM_OPTS" | head -n1); Nx=${Nx:-1}
Ny=$(awk '$1=="-Ny"{print $2}' "$GEOM_OPTS" | head -n1); Ny=${Ny:-1}
Nz=$(awk '$1=="-Nz"{print $2}' "$GEOM_OPTS" | head -n1); Nz=${Nz:-1}
total_dofs=$(( 4 * Nx * Ny * Nz ))
TARGET=10000; MAX_PER_NODE=32
nprocs=$(( (total_dofs + TARGET - 1) / TARGET ))
(( nprocs < 1 )) && nprocs=1
tasks_per_node=$nprocs
(( tasks_per_node > MAX_PER_NODE )) && tasks_per_node=$MAX_PER_NODE
nnodes=$(( (nprocs + tasks_per_node - 1) / tasks_per_node ))
(( nnodes < 1 )) && nnodes=1

# ---------------------------------------------------------------------------
# Compile once on the login node
# ---------------------------------------------------------------------------
echo "--- Compiling ---"
cd "$PROJ"
make clean && make all
echo "✅ Build complete."

# ---------------------------------------------------------------------------
# Shared parent folder
# ---------------------------------------------------------------------------
TS=$(date +%Y-%m-%d__%H.%M.%S)
BATCH_PARENT="$SCRATCH/permafrost/${GEOM}/${TS}_mob_sweep${TAG:+_$TAG}"
mkdir -p "$BATCH_PARENT"

echo ""
echo "============================================================"
echo "  mob_sub sweep — $GEOM × $EXP"
echo "  Parent dir  : $BATCH_PARENT"
echo "  DOFs        : $total_dofs  →  nprocs=$nprocs  nodes=$nnodes"
echo "  Values      : ${MOB_VALUES[*]}"
echo "============================================================"

# ---------------------------------------------------------------------------
# One sbatch job per mob_sub value — each gets its own subfolder
# ---------------------------------------------------------------------------
for mob in "${MOB_VALUES[@]}"; do
    out_dir="$BATCH_PARENT/mob_${mob}"
    mkdir -p "$out_dir"
    cp "$SOLVER_OPTS" "$out_dir/"
    cp "$EXP_OPTS"    "$out_dir/"

    geom_copy="$out_dir/$(basename "$GEOM_OPTS")"
    cp "$GEOM_OPTS" "$geom_copy"

    # kp_default: strip the frozen -mob_sub line so the executable uses K&P
    # matching at the current eps; all other cases pass -mob_sub explicitly.
    if [[ "$mob" == "kp_default" ]]; then
        sed -i '/^-mob_sub/d' "$geom_copy"
        mob_arg=""
        run_geom="$geom_copy"
    else
        mob_arg="-mob_sub $mob"
        run_geom="$GEOM_OPTS"
    fi

    sbatch \
        --job-name="pf_mob_${mob}" \
        --nodes="$nnodes" \
        --ntasks="$nprocs" \
        --ntasks-per-node="$tasks_per_node" \
        --cpus-per-task=1 \
        --partition=expansion \
        --account=rubyfu \
        --time=0-12:00:00 \
        --mem-per-cpu=2G \
        --output="${out_dir}/slurm-%j.out" \
        --error="${out_dir}/slurm-%j.err" \
        --mail-user=jbaglino@caltech.edu \
        --mail-type=END,FAIL \
        --constraint='icelake|skylake|cascadelake' \
        "${sbatch_extra[@]}" \
        --wrap="export folder=$out_dir && srun -n $nprocs $EXEC \
            -options_file $SOLVER_OPTS \
            -options_file $run_geom \
            -options_file $EXP_OPTS \
            -output_path  $out_dir \
            $mob_arg \
            2>&1 | tee $out_dir/outp.txt"

    echo "  → submitted mob_sub=$mob  →  $out_dir"
done

echo ""
echo "  All jobs submitted. Check with: squeue -u \$USER"
echo "  Results: $BATCH_PARENT"
