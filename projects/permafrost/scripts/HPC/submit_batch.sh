#!/usr/bin/env bash
# =============================================================================
# submit_batch.sh — Submit a batch of permafrost simulations to SLURM,
# with all outputs going into a single timestamped parent folder for easy
# bulk download.
#
# Mirrors the local scripts/Studio/run_batch_tests.sh interface so the same
# test specs work on both. Each test is submitted as an independent sbatch
# job (so they run in parallel on the cluster), but they all write into:
#
#   $SCRATCH/permafrost/batch_<timestamp>[_<tag>]/<geom>__<exp>/
#
# Each sbatch call:
#   - sizes its allocation to the geometry's grid (TARGET_DOFS_PER_CORE)
#   - calls run_permafrost.sh as the actual SLURM script
#   - sets BATCH_OUT_DIR so run_permafrost.sh writes into the shared parent
#   - sets SKIP_COMPILE=1 since this script builds once on the submission host
#
# Usage (run from project root):
#   ./scripts/HPC/submit_batch.sh --tag mytag \
#       --tests "1D_separated_grains:1day_T-20_h1.00,2D_separated_grains:30day_T-5_h1.00"
#
#   ./scripts/HPC/submit_batch.sh --tag mytag --tests-file tests.txt
#
# --extra-opts forwards a single quoted string of permafrost CLI flags to
# EVERY fanned-out job (appended after the three -options_file flags, same
# as submit_permafrost.sh's `-- ...` convention, so they override anything
# set in the opts files):
#   ./scripts/HPC/submit_batch.sh --tag mytag --tests "..." --extra-opts "-beta_sub0 1.4e3"
#
# Extra sbatch flags can be appended after --:
#   ./scripts/HPC/submit_batch.sh --tag mytag --tests "..." -- --time=0-04:00:00
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/run_permafrost.sh"
INPUTS_DIR="$PROJECT_ROOT/inputs"
GEOMETRY_DIR="$INPUTS_DIR/geometry"
EXPERIMENT_DIR="$INPUTS_DIR/experiment"
SOLVER_OPTS="$INPUTS_DIR/solver.opts"

# Resource-sizing parameters. FOUR copies of TARGET_DOFS_PER_CORE exist -- keep
# them in sync:
#   scripts/Studio/run_permafrost.sh :: compute_optimal_nprocs
#   scripts/HPC/run_permafrost.sh    :: compute_optimal_nprocs
#   scripts/HPC/submit_permafrost.sh
#   scripts/HPC/submit_batch.sh      :: here
# (The other three said "three copies" and omitted this file, which is exactly
# how it was left behind at the old 10000 when the rest moved to 40000 on
# 2026-07-12 -- a 4x over-allocation, compounded to 5.3x by the hardcoded
# dof=4 in compute_alloc. Fixed 2026-07-15.)
TARGET_DOFS_PER_CORE=40000
MAX_TASKS_PER_NODE=32

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------
tag=""
tests_arg=""
tests_file=""
sbatch_extra=()
extra_opts=()

usage() {
    sed -n '2,28p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)         tag="$2"; shift 2 ;;
        --tests)       tests_arg="$2"; shift 2 ;;
        --tests-file)  tests_file="$2"; shift 2 ;;
        --extra-opts)  read -ra extra_opts <<< "$2"; shift 2 ;;
        --)            shift; sbatch_extra=("$@"); break ;;
        -h|--help)     usage ;;
        *)             echo "Unknown argument: $1"; usage ;;
    esac
done

# Build the list of tests
TESTS=()
if [[ -n "$tests_arg" ]]; then
    IFS=',' read -ra TESTS <<< "$tests_arg"
    # Trim leading/trailing whitespace from each entry (so `--tests "a, b,c"`
    # and indented backslash-continuation lines both work).
    for i in "${!TESTS[@]}"; do
        s="${TESTS[$i]}"
        s="${s#"${s%%[![:space:]]*}"}"   # strip leading whitespace
        s="${s%"${s##*[![:space:]]}"}"   # strip trailing whitespace
        TESTS[$i]="$s"
    done
elif [[ -n "$tests_file" ]]; then
    if [[ ! -f "$tests_file" ]]; then
        echo "❌ Tests file not found: $tests_file"
        exit 1
    fi
    # One "geom:exp" per line, ignore blank lines and # comments. Robustly
    # strip leading/trailing whitespace, surrounding quotes, and trailing
    # commas so the file format is forgiving (e.g. accidentally indented
    # heredocs, copy-pasted backslash continuations, quoted EOF markers).
    while IFS= read -r line; do
        line="${line%%#*}"                          # strip # comments
        line="${line%$'\r'}"                        # strip CR (Windows line endings)
        line="${line#"${line%%[![:space:]]*}"}"     # strip leading whitespace
        line="${line%"${line##*[![:space:]]}"}"     # strip trailing whitespace
        line="${line%,}"                            # strip trailing comma
        # Strip a single pair of surrounding quotes if present
        if [[ "$line" =~ ^\"(.*)\"$ ]]; then line="${BASH_REMATCH[1]}"; fi
        if [[ "$line" =~ ^\'(.*)\'$ ]]; then line="${BASH_REMATCH[1]}"; fi
        # Skip anything that isn't of the form geom:exp
        [[ -z "$line" ]] && continue
        [[ "$line" != *:* ]] && continue
        TESTS+=("$line")
    done < "$tests_file"
else
    echo "❌ Must supply --tests \"g1:e1,g2:e2,...\" or --tests-file <file>"
    usage
fi

if [[ ${#TESTS[@]} -eq 0 ]]; then
    echo "❌ No tests specified."
    exit 1
fi

cd "$PROJECT_ROOT"

# ---------------------------------------------------------------------------
# Build once (sets SKIP_COMPILE=1 for each fanned-out job — see compile_code
# in run_permafrost.sh for the race-condition rationale).
# ---------------------------------------------------------------------------
echo ""
echo "--- Building permafrost on submission host ---"
if ! make all; then
    echo "❌ Build failed. Fix the build before submitting jobs."
    exit 1
fi
if [[ ! -x ./permafrost ]]; then
    echo "❌ ./permafrost still missing after make all."
    exit 1
fi
echo "✅ Build complete."

# ---------------------------------------------------------------------------
# Create the shared parent batch folder (under $SCRATCH on HPC,
# $PROJECT_ROOT/scratch as fallback for local testing).
# ---------------------------------------------------------------------------
TS=$(date +%Y-%m-%d__%H.%M.%S)
batch_name="batch_${TS}${tag:+_$tag}"

if [[ -d "${SCRATCH:-}" ]]; then
    BATCH_PARENT="$SCRATCH/permafrost/$batch_name"
else
    BATCH_PARENT="$PROJECT_ROOT/scratch/$batch_name"
fi
mkdir -p "$BATCH_PARENT"

echo "============================================================"
echo "  Permafrost batch submission"
echo "  Tag         : ${tag:-<none>}"
echo "  Tests       : ${#TESTS[@]}"
echo "  Parent dir  : $BATCH_PARENT"
echo "  Target      : ${TARGET_DOFS_PER_CORE} DoFs/core (max ${MAX_TASKS_PER_NODE}/node)"
echo "============================================================"

# ---------------------------------------------------------------------------
# Stage shared assets at the batch level for reproducibility + easy download
# ---------------------------------------------------------------------------
mkdir -p "$BATCH_PARENT/inputs_snapshot/geometry"
mkdir -p "$BATCH_PARENT/inputs_snapshot/experiment"
mkdir -p "$BATCH_PARENT/src_snapshot"
cp "$INPUTS_DIR/solver.opts"           "$BATCH_PARENT/inputs_snapshot/"   2>/dev/null || true
cp -r "$GEOMETRY_DIR"/*.opts           "$BATCH_PARENT/inputs_snapshot/geometry/"   2>/dev/null || true
cp -r "$EXPERIMENT_DIR"/*.opts         "$BATCH_PARENT/inputs_snapshot/experiment/" 2>/dev/null || true
for ext in c h; do
    cp "$PROJECT_ROOT/src/"*.$ext     "$BATCH_PARENT/src_snapshot/"      2>/dev/null || true
done
cp -r "$PROJECT_ROOT/include"          "$BATCH_PARENT/src_snapshot/"      2>/dev/null || true
cp    "$PROJECT_ROOT/makefile"         "$BATCH_PARENT/src_snapshot/"      2>/dev/null || true
cp    "$PROJECT_ROOT/postprocess"      -r  "$BATCH_PARENT/"               2>/dev/null || true
cp    "${BASH_SOURCE[0]}"              "$BATCH_PARENT/submit_batch.sh"

# Copy the local-postprocessing helper (for after the user downloads the batch)
if [[ -f "$PROJECT_ROOT/postprocess/run_batch_postprocess.sh" ]]; then
    cp "$PROJECT_ROOT/postprocess/run_batch_postprocess.sh" "$BATCH_PARENT/"
fi

# ---------------------------------------------------------------------------
# Per-geometry allocation sizer.
# Echoes "nprocs nnodes tasks_per_node total_dofs"
# ---------------------------------------------------------------------------
compute_alloc() {
    local geom_file="$1"
    local nx ny nz dof

    # dof from solver.opts, NOT hardcoded. This used to be a literal 4 while
    # solver.opts sets -dof 3, inflating every allocation by 4/3 on top of the
    # stale DoFs/core target. The other three sizers all read it from the file.
    dof=$(awk '$1=="-dof"{print $2}' "$SOLVER_OPTS" 2>/dev/null | head -n1)
    [[ -z "${dof:-}" ]] && dof=4

    # -geom_file meshes override -Nx/-Ny/-Nz; read the grid from the
    # "# DOF_GRID: nx ny [nz]" comment, matching submit_permafrost.sh.
    if grep -q "^-geom_file" "$geom_file"; then
        read -r nx ny nz <<< "$(awk '$1=="#" && $2=="DOF_GRID:"{print $3, $4, $5}' "$geom_file" | head -n1)"
    else
        nx=$(awk '$1=="-Nx"{print $2}' "$geom_file" | head -n1)
        ny=$(awk '$1=="-Ny"{print $2}' "$geom_file" | head -n1)
        nz=$(awk '$1=="-Nz"{print $2}' "$geom_file" | head -n1)
    fi
    nx=${nx:-1}; ny=${ny:-1}; nz=${nz:-1}

    local total_dofs=$((dof * nx * ny * nz))
    local nprocs=$(( (total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE ))
    (( nprocs < 1 )) && nprocs=1

    local tasks_per_node=$nprocs
    (( tasks_per_node > MAX_TASKS_PER_NODE )) && tasks_per_node=$MAX_TASKS_PER_NODE

    local nnodes=$(( (nprocs + tasks_per_node - 1) / tasks_per_node ))
    (( nnodes < 1 )) && nnodes=1

    echo "$nprocs $nnodes $tasks_per_node $total_dofs"
}

# ---------------------------------------------------------------------------
# Submit one job: sbatch run_permafrost.sh <geom> <exp> <tag>
# with BATCH_OUT_DIR pointing at the shared parent so all jobs end up there.
# ---------------------------------------------------------------------------
N_SUBMITTED=0
N_SKIPPED=0
submit_one() {
    local spec="$1"
    local geom="${spec%%:*}"
    local exp="${spec##*:}"
    if [[ "$geom" == "$exp" || -z "$geom" || -z "$exp" ]]; then
        echo "⚠ Invalid test spec (expected geom:exp): $spec"
        ((N_SKIPPED++)) || true
        return
    fi

    local geom_file="$GEOMETRY_DIR/${geom}.opts"
    local exp_file="$EXPERIMENT_DIR/${exp}.opts"
    if [[ ! -f "$geom_file" ]]; then
        echo "⚠ Skipping $spec — geometry file not found: $geom_file"
        ((N_SKIPPED++)) || true
        return
    fi
    if [[ ! -f "$exp_file" ]]; then
        echo "⚠ Skipping $spec — experiment file not found: $exp_file"
        ((N_SKIPPED++)) || true
        return
    fi

    local job_name="${geom}__${exp}"
    local nprocs nnodes tasks_per_node total_dofs
    read -r nprocs nnodes tasks_per_node total_dofs < <(compute_alloc "$geom_file")

    printf "→ %-45s DoFs=%-8d nprocs=%-3d nodes=%-2d tasks/node=%d\n" \
        "$job_name" "$total_dofs" "$nprocs" "$nnodes" "$tasks_per_node"

    sbatch --job-name="$job_name" \
           --nodes="$nnodes" \
           --ntasks="$nprocs" \
           --ntasks-per-node="$tasks_per_node" \
           --export=ALL,SKIP_COMPILE=1,BATCH_OUT_DIR="$BATCH_PARENT" \
           "${sbatch_extra[@]}" \
           "$RUN_SCRIPT" "$geom" "$exp" "$tag" "${extra_opts[@]}"
    ((N_SUBMITTED++)) || true
}

# ---------------------------------------------------------------------------
# Fan out
# ---------------------------------------------------------------------------
for spec in "${TESTS[@]}"; do
    submit_one "$spec"
done

echo ""
echo "============================================================"
echo "  Parsed     : ${#TESTS[@]} test specs"
echo "  Submitted  : $N_SUBMITTED jobs to SLURM"
echo "  Skipped    : $N_SKIPPED (file-not-found or malformed)"
echo "  Parent dir : $BATCH_PARENT"
echo "  Check with: squeue -u \$USER"
echo "  Once all jobs are done, download the parent dir, then run:"
echo "    bash $BATCH_PARENT/run_batch_postprocess.sh"
echo "============================================================"
