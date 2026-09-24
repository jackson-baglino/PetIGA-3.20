#!/usr/bin/env bash
# =============================================================================
# submit_batch.sh — Submit a batch of enceladus_dsm simulations to SLURM,
# with all outputs going into a single timestamped parent folder for easy
# bulk download.
#
# Mirrors the local scripts/Studio/run_batch_tests.sh interface so the same
# test specs work on both. Each test is submitted as an independent sbatch
# job (so they run in parallel on the cluster), but they all write into:
#
#   $SCRATCH/enceladus_DSM/batch_<timestamp>[_<tag>]/<geom>__<exp>/
#
# Each sbatch call:
#   - sizes its allocation to the geometry's grid (TARGET_DOFS_PER_CORE)
#   - calls run_enceladus.sh as the actual SLURM script
#   - sets BATCH_OUT_DIR so run_enceladus.sh writes into the shared parent
#   - sets SKIP_COMPILE=1 since this script builds once on the submission host
#
# Usage (run from project root):
#   ./scripts/HPC/submit_batch.sh --tag mytag \
#       --tests "1D_separated_grains:base_T-20_h1.00_1d,2D_separated_grains:base_T-5_h1.00_30d"
#
#   ./scripts/HPC/submit_batch.sh --tag mytag --tests-file tests.txt
#
# --extra-opts forwards a single quoted string of enceladus_dsm CLI flags to
# EVERY fanned-out job (appended after the three -options_file flags, same
# as submit_enceladus.sh's `-- ...` convention, so they override anything
# set in the opts files):
#   ./scripts/HPC/submit_batch.sh --tag mytag --tests "..." --extra-opts "-beta_sub0 1.4e3"
#
# PER-JOB options: a test spec may carry a THIRD field, geom:exp:<opts>, which
# is appended after --extra-opts for that job only. --extra-opts goes to every
# job, so it cannot carry anything that differs between them -- a per-run
# -keff_replay directory being the case this was added for.
#
# Use --tests-file for these: one spec per line, so the opts may contain spaces
# and colons. (--tests splits on commas, so a third field there must not.)
#
#   packing_2D_..._seed1_...:snow_T-20_h1.00_30d:--label tensor -keff_replay /path/seed1
#
# --label <name> inside that third field is consumed HERE, not passed to the
# solver. It suffixes the job name and the output subfolder, which is what lets
# two jobs share a geometry and an experiment -- the same run replayed under two
# conductivity laws would otherwise both write to <geom>__<exp>/ and clobber
# each other. Omitted, it defaults to the spec's position, j01, j02, ...
#
# Extra sbatch flags can be appended after --:
#   ./scripts/HPC/submit_batch.sh --tag mytag --tests "..." -- --time=0-04:00:00
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/run_enceladus.sh"
INPUTS_DIR="$PROJECT_ROOT/inputs"
GEOMETRY_DIR="$INPUTS_DIR/geometry"
EXPERIMENT_DIR="$INPUTS_DIR/experiment"
SOLVER_OPTS="$INPUTS_DIR/solver.opts"

# Resource-sizing parameters. FOUR copies of TARGET_DOFS_PER_CORE exist -- keep
# them in sync:
#   scripts/Studio/run_enceladus.sh :: compute_optimal_nprocs
#   scripts/HPC/run_enceladus.sh    :: compute_optimal_nprocs
# TARGET_DOFS_PER_CORE and MAX_TASKS_PER_NODE are sourced from
# scripts/lib/alloc.sh (single source of truth; see rationale there).
source "$PROJECT_ROOT/scripts/lib/provenance.sh"
source "$PROJECT_ROOT/scripts/lib/alloc.sh"
source "$PROJECT_ROOT/scripts/lib/opts.sh"

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------
tag=""
tests_arg=""
tests_file=""
sbatch_extra=()
extra_opts=()

usage() {
    sed -n '2,44p' "$0"
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
# in run_enceladus.sh for the race-condition rationale).
# ---------------------------------------------------------------------------
echo ""
echo "--- Building enceladus_dsm on submission host ---"
if ! make all; then
    echo "❌ Build failed. Fix the build before submitting jobs."
    exit 1
fi
if [[ ! -x ./enceladus_dsm ]]; then
    echo "❌ ./enceladus_dsm still missing after make all."
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
    BATCH_PARENT="$SCRATCH/enceladus_DSM/$batch_name"
else
    BATCH_PARENT="$PROJECT_ROOT/scratch/$batch_name"
fi
mkdir -p "$BATCH_PARENT"

echo "============================================================"
echo "  Enceladus DSM batch submission"
echo "  Tag         : ${tag:-<none>}"
print_repo_provenance "$PROJECT_ROOT"
echo "  Tests       : ${#TESTS[@]}"
echo "  Parent dir  : $BATCH_PARENT"
echo "  Target      : ${TARGET_DOFS_PER_CORE} DoFs/core (max ${MAX_TASKS_PER_NODE}/node)"
echo "============================================================"

# ---------------------------------------------------------------------------
# Stage shared assets at the batch level for reproducibility + easy download
# ---------------------------------------------------------------------------
mkdir -p "$BATCH_PARENT/src_snapshot"
# The whole inputs/ tree, named `inputs` rather than `inputs_snapshot`, and
# minus scratch/. Two reasons for both choices:
#   * run_batch_measure.sh runs $BATCH_PARENT/postprocess/*.py, and those
#     locate the experimental series as Path(__file__).parent.parent /
#     "inputs/validation/...". Under the old name that lookup missed and the
#     figures came out with no data on them.
#   * the piecemeal copy only took solver.opts + geometry/ + experiment/, so
#     validation/ and the inputs README were never in the snapshot at all.
# scratch/ is 70 MB of retired files pending deletion; the rest is ~4 MB.
# (`inputs_snapshot` is still skipped by the batch iterators, so older batches
# keep working.)
cp -r "$INPUTS_DIR"                    "$BATCH_PARENT/inputs"             2>/dev/null || true
rm -rf "$BATCH_PARENT/inputs/scratch"
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
    local geom_file="$1"; shift
    # Remaining args: the job's solver flags (--extra-opts + per-job), which
    # are appended after the opts files and so win at solve time. The
    # allocation must read the same values, or a per-job -Nx 4096 on a
    # geometry file that says 256 is sized for 256.
    local job_opts=("$@")
    local nx ny nz dof

    # dof from solver.opts, NOT hardcoded. This used to be a literal 4 while
    # solver.opts sets -dof 3, inflating every allocation by 4/3 on top of the
    # stale DoFs/core target. The other three sizers all read it from the file.
    dof=$(awk '$1=="-dof"{print $2}' "$SOLVER_OPTS" 2>/dev/null | head -n1)
    [[ -z "${dof:-}" ]] && dof=3

    # -geom_file meshes override -Nx/-Ny/-Nz; read the grid from the
    # "# DOF_GRID: nx ny [nz]" comment, matching submit_enceladus.sh.
    if grep -q "^-geom_file" "$geom_file"; then
        read -r nx ny nz <<< "$(awk '$1=="#" && $2=="DOF_GRID:"{print $3, $4, $5}' "$geom_file" | head -n1)"
    else
        nx=$(awk '$1=="-Nx"{print $2}' "$geom_file" | head -n1)
        ny=$(awk '$1=="-Ny"{print $2}' "$geom_file" | head -n1)
        nz=$(awk '$1=="-Nz"{print $2}' "$geom_file" | head -n1)
    fi
    local i keff_only=0
    for ((i = 0; i < ${#job_opts[@]}; i++)); do
        case "${job_opts[$i]}" in
            -Nx) nx="${job_opts[$((i+1))]:-$nx}" ;;
            -Ny) ny="${job_opts[$((i+1))]:-$ny}" ;;
            -Nz) nz="${job_opts[$((i+1))]:-$nz}" ;;
            -keff_only) [[ "${job_opts[$((i+1))]:-1}" != 0 ]] && keff_only=1 ;;
        esac
    done
    # -keff_only never builds the phase-field Jacobian: the cost is the scalar
    # corrector solve, one unknown per node, so size on that.
    (( keff_only )) && dof=1
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
# Submit one job: sbatch run_enceladus.sh <geom> <exp> <tag>
# with BATCH_OUT_DIR pointing at the shared parent so all jobs end up there.
# ---------------------------------------------------------------------------
N_SUBMITTED=0
N_SKIPPED=0
submit_one() {
    local spec="$1"
    local idx="$2"
    # geom:exp[:per-job opts]. The third field is the REMAINDER of the line,
    # colons included, so a -keff_replay path with a colon in it survives.
    local geom exp perjob_str
    IFS=':' read -r geom exp perjob_str <<< "$spec"
    if [[ -z "$geom" || -z "$exp" ]]; then
        echo "⚠ Invalid test spec (expected geom:exp[:opts]): $spec"
        ((N_SKIPPED++)) || true
        return
    fi

    # Split the per-job options and pull out --label, which is ours, not the
    # solver's: it disambiguates the job name and the output subfolder when
    # two jobs share a geometry and an experiment (e.g. the same run replayed
    # under two conductivity laws). Without it they would both land in
    # $BATCH_OUT_DIR/<geom>__<exp>/ and overwrite each other's staged inputs.
    local perjob=() label=""
    if [[ -n "${perjob_str:-}" ]]; then
        read -ra perjob <<< "$perjob_str"
        local keep=() i=0
        while [[ $i -lt ${#perjob[@]} ]]; do
            if [[ "${perjob[$i]}" == "--label" ]]; then
                label="${perjob[$((i+1))]:-}"
                i=$((i+2))
            else
                keep+=("${perjob[$i]}")
                i=$((i+1))
            fi
        done
        perjob=(${keep[@]+"${keep[@]}"})
        [[ -z "$label" ]] && label="j$(printf '%02d' "$idx")"
    fi

    # Resolve through the shared helper (scripts/lib/opts.sh): .opts live in
    # per-family subdirectories while the spec carries a bare name.
    local geom_file exp_file
    if ! geom_file="$(require_opts "$GEOMETRY_DIR" "$geom" geometry)"; then
        echo "⚠ Skipping $spec — geometry not found"
        ((N_SKIPPED++)) || true
        return
    fi
    if ! exp_file="$(require_opts "$EXPERIMENT_DIR" "$exp" experiment)"; then
        echo "⚠ Skipping $spec — experiment not found"
        ((N_SKIPPED++)) || true
        return
    fi

    local job_name="${geom}__${exp}${label:+__${label}}"
    local nprocs nnodes tasks_per_node total_dofs
    read -r nprocs nnodes tasks_per_node total_dofs < <(compute_alloc "$geom_file" ${extra_opts[@]+"${extra_opts[@]}"} ${perjob[@]+"${perjob[@]}"})

    printf "→ %-45s DoFs=%-8d nprocs=%-3d nodes=%-2d tasks/node=%d\n" \
        "$job_name" "$total_dofs" "$nprocs" "$nnodes" "$tasks_per_node"
    [[ ${#perjob[@]} -gt 0 ]] && printf "    per-job opts: %s\n" "${perjob[*]}"

    sbatch --job-name="$job_name" \
           --nodes="$nnodes" \
           --ntasks="$nprocs" \
           --ntasks-per-node="$tasks_per_node" \
           --export=ALL,SKIP_COMPILE=1,BATCH_OUT_DIR="$BATCH_PARENT",BATCH_JOB_LABEL="$label" \
           ${sbatch_extra[@]+"${sbatch_extra[@]}"} \
           "$RUN_SCRIPT" "$geom" "$exp" "$tag" \
           ${extra_opts[@]+"${extra_opts[@]}"} ${perjob[@]+"${perjob[@]}"}
    ((N_SUBMITTED++)) || true
}

# ---------------------------------------------------------------------------
# Fan out
# ---------------------------------------------------------------------------
spec_idx=0
for spec in "${TESTS[@]}"; do
    spec_idx=$((spec_idx + 1))
    submit_one "$spec" "$spec_idx"
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
