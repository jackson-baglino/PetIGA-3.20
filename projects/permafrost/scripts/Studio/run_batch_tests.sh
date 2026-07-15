#!/usr/bin/env bash
# =============================================================================
# run_batch_tests.sh — Run a curated set of tests sequentially into one folder
#
# Location: $PETIGA_DIR/projects/permafrost/scripts/Studio/run_batch_tests.sh
#
# All runs go under a single timestamped parent folder:
#   $RESULTS_BASE/batch_<timestamp>/<test_name>/
#
# Each test is composed from three opts files:
#   inputs/solver.opts                  (numerical defaults — always)
#   inputs/geometry/<geom>.opts         (mesh, domain, IC, geometry-specific tweaks)
#   inputs/experiment/<exp>.opts        (t_final, T, humidity, ...)
#
# Test list entries are "<geometry>:<experiment>" (without .opts extension).
# The per-test output directory name is "<geometry>__<experiment>".
#
# Usage:
#   ./scripts/Studio/run_batch_tests.sh                                  # default set
#   ./scripts/Studio/run_batch_tests.sh --skip-1d                        # 2D-only
#   ./scripts/Studio/run_batch_tests.sh --tag mylabel
#   ./scripts/Studio/run_batch_tests.sh --tests "2D_touching_grains:1day_T-20_h1.00,2D_single_ice:1day_T-20_h0.95"
#
# Or run a single one-off:
#   ./scripts/Studio/run_batch_tests.sh --geom 2D_touching_grains --exp 1day_T-20_h1.00
# =============================================================================

set -uo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PETIGA_DIR}/projects/permafrost"
EXEC="$PROJECT_ROOT/permafrost"
INPUTS_DIR="$PROJECT_ROOT/inputs"
SOLVER_OPTS="$INPUTS_DIR/solver.opts"
GEOMETRY_DIR="$INPUTS_DIR/geometry"
EXPERIMENT_DIR="$INPUTS_DIR/experiment"
POSTPROCESS="$PROJECT_ROOT/postprocess"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
RESULTS_BASE="/Users/jacksonbaglino/SimulationResults/permafrost/scratch"

MAX_LOCAL_CORES=12
# One of FIVE copies -- keep in sync (Studio/run, Studio/run_batch_tests,
# HPC/run, HPC/submit, HPC/submit_batch). This file and HPC/submit_batch.sh
# were both left at the old 10000 when the rest moved to 40000 on 2026-07-12;
# fixed 2026-07-15. Locally the hw-core cap hid it; on HPC it did not.
TARGET_DOFS_PER_CORE=40000

if [ ! -f "$PROJECT_ROOT/makefile" ] && [ ! -f "$PROJECT_ROOT/Makefile" ]; then
    echo "❌ Could not find makefile at $PROJECT_ROOT"
    exit 1
fi
if [ ! -f "$SOLVER_OPTS" ]; then
    echo "❌ Could not find solver opts at $SOLVER_OPTS"
    exit 1
fi

PYTHON=$(command -v python3 || command -v python || true)

# ---------------------------------------------------------------------------
# Default test list. Format: "<geometry>:<experiment>"
# (geometry and experiment names map to inputs/geometry/<g>.opts and
#  inputs/experiment/<e>.opts).
# ---------------------------------------------------------------------------
DEFAULT_TESTS=(
    # ---- 1D tests at h=0.95, T=-20 ----
    "1D_ice_slab:1day_T-20_h0.95"
    "1D_ice_slab_hires:1day_T-20_h0.95"
    "1D_single_ice:1day_T-20_h0.95"
    "1D_single_ice_hires:1day_T-20_h0.95"
    "1D_ice_sed_pair:1day_T-20_h0.95"
    "1D_ice_sed_pair_hires:1day_T-20_h0.95"
    "1D_separated_grains:1day_T-20_h0.95"
    "1D_separated_grains_hires:1day_T-20_h0.95"
    "1D_touching_grains:1day_T-20_h0.95"
    "1D_touching_grains_hires:1day_T-20_h0.95"

    # ---- 2D tests at saturation (multi-grain physics — Ostwald ripening, sintering) ----
    "2D_ice_sed_pair:1day_T-20_h1.00"
    "2D_ice_sed_pair_hires:1day_T-20_h1.00"
    "2D_touching_grains:1day_T-20_h1.00"
    "2D_touching_grains_hires:1day_T-20_h1.00"
    "2D_separated_grains:1day_T-20_h1.00"
    "2D_separated_grains_hires:1day_T-20_h1.00"
    "2D_single_sed:1day_T-20_h1.00"

    # ---- 2D tests at undersaturation (bulk sublimation driving force) ----
    "2D_single_ice:1day_T-20_h0.95"
    "2D_single_ice_hires:1day_T-20_h0.95"
    "2D_ice_slab:1day_T-20_h0.95"
    "2D_ice_slab_hires:1day_T-20_h0.95"
)

# ---------------------------------------------------------------------------
# Arg parsing
# ---------------------------------------------------------------------------
tag=""
skip_1d=false
custom_tests=""
oneshot_geom=""
oneshot_exp=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)
            tag="$2"; shift 2 ;;
        --skip-1d)
            skip_1d=true; shift ;;
        --tests)
            custom_tests="$2"; shift 2 ;;
        --geom|--geometry)
            oneshot_geom="$2"; shift 2 ;;
        --exp|--experiment)
            oneshot_exp="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,28p' "$0"; exit 0 ;;
        *)
            echo "Unknown argument: $1"
            exit 1 ;;
    esac
done

# Build the actual test list
TESTS=()
if [[ -n "$oneshot_geom" || -n "$oneshot_exp" ]]; then
    if [[ -z "$oneshot_geom" || -z "$oneshot_exp" ]]; then
        echo "❌ --geom and --exp must be given together"
        exit 1
    fi
    TESTS+=("${oneshot_geom}:${oneshot_exp}")
elif [[ -n "$custom_tests" ]]; then
    IFS=',' read -ra TESTS <<< "$custom_tests"
    # Trim leading/trailing whitespace from each entry so indented backslash
    # continuations and `"a, b, c"` style spaces don't break the geom:exp split.
    for _i in "${!TESTS[@]}"; do
        _s="${TESTS[$_i]}"
        _s="${_s#"${_s%%[![:space:]]*}"}"
        _s="${_s%"${_s##*[![:space:]]}"}"
        TESTS[$_i]="$_s"
    done
    unset _i _s
else
    for t in "${DEFAULT_TESTS[@]}"; do
        if $skip_1d && [[ "$t" == 1D_* ]]; then
            continue
        fi
        TESTS+=("$t")
    done
fi

if [[ ${#TESTS[@]} -eq 0 ]]; then
    echo "❌ No tests selected"
    exit 1
fi

# ---------------------------------------------------------------------------
# Set up batch folder
# ---------------------------------------------------------------------------
TS=$(date +%Y-%m-%d__%H.%M.%S)
batch_name="batch_${TS}${tag:+_$tag}"
BATCH_DIR="$RESULTS_BASE/$batch_name"
mkdir -p "$BATCH_DIR"

SUMMARY="$BATCH_DIR/SUMMARY.txt"
{
    echo "Batch run: $batch_name"
    echo "Started:   $(date)"
    echo "Project:   $PROJECT_ROOT"
    echo "Tests:     ${#TESTS[@]}"
    echo ""
    echo "STATUS legend:"
    echo "  OK        Simulation completed normally."
    echo "  OK?       Completed, but t_final < 60 s (likely a leftover debug value)."
    echo "  FAIL(N)   Simulation aborted with exit code N."
    echo "  MISSING   .opts file not found."
    echo ""
    printf "%-45s | %-9s | %-7s | %-9s\n" "TEST" "STATUS" "NPROCS" "TIME (s)"
    printf "%-45s-+-%-9s-+-%-7s-+-%-9s\n" "$(printf -- '-%.0s' {1..45})" "---------" "-------" "---------"
} > "$SUMMARY"

# ---------------------------------------------------------------------------
# Compile once
# ---------------------------------------------------------------------------
echo "=========================================================================="
echo "  Batch run: $batch_name"
echo "  Output:    $BATCH_DIR"
echo "  Tests:     ${#TESTS[@]}"
echo "=========================================================================="
echo ""
echo "--- Compiling ---"
cd "$PROJECT_ROOT" && make
if [ ! -f "$EXEC" ]; then
    echo "❌ Build failed: $EXEC not found"
    exit 1
fi
echo "✅ Build complete."

# Stage shared assets at the batch level (source + the three input dirs)
mkdir -p "$BATCH_DIR/src_snapshot"
cp "$SOLVER_OPTS" "$BATCH_DIR/" 2>/dev/null || true
mkdir -p "$BATCH_DIR/inputs_snapshot"
cp -r "$GEOMETRY_DIR"   "$BATCH_DIR/inputs_snapshot/" 2>/dev/null || true
cp -r "$EXPERIMENT_DIR" "$BATCH_DIR/inputs_snapshot/" 2>/dev/null || true
for ext in c h; do
    cp "$PROJECT_ROOT/src/"*.$ext "$BATCH_DIR/src_snapshot/" 2>/dev/null || true
done
cp -r "$PROJECT_ROOT/include" "$BATCH_DIR/src_snapshot/" 2>/dev/null || true
cp "$PROJECT_ROOT/makefile" "$BATCH_DIR/src_snapshot/" 2>/dev/null || true
cp "${BASH_SOURCE[0]}" "$BATCH_DIR/run_batch_tests.sh"

# ---------------------------------------------------------------------------
# nprocs helper
# ---------------------------------------------------------------------------
choose_nprocs() {
    local opts="$1"
    local Nx Ny Nz total dofs cap hw dof
    Nx=$(awk '$1=="-Nx"{print $2}' "$opts" | head -n1)
    Ny=$(awk '$1=="-Ny"{print $2}' "$opts" | head -n1)
    Nz=$(awk '$1=="-Nz"{print $2}' "$opts" | head -n1)
    [[ -z "${Nx:-}" ]] && Nx=1
    [[ -z "${Ny:-}" ]] && Ny=1
    [[ -z "${Nz:-}" ]] && Nz=1
    # dof from solver.opts, not a hardcoded 4 (solver.opts sets -dof 3).
    dof=$(awk '$1=="-dof"{print $2}' "$SOLVER_OPTS" 2>/dev/null | head -n1)
    [[ -z "${dof:-}" ]] && dof=4
    dofs=$((dof * Nx * Ny * Nz))
    local n=$(((dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE))
    (( n < 1 )) && n=1
    hw=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo $MAX_LOCAL_CORES)
    cap=$(( hw < MAX_LOCAL_CORES ? hw : MAX_LOCAL_CORES ))
    (( n > cap )) && n=$cap
    echo "$n"
}

# ---------------------------------------------------------------------------
# Per-test runner. spec format: "<geometry>:<experiment>"
# ---------------------------------------------------------------------------
run_one_test() {
    local spec="$1"
    local geom="${spec%%:*}"
    local exp="${spec##*:}"
    if [[ "$geom" == "$exp" || -z "$geom" || -z "$exp" ]]; then
        echo "⚠️  Invalid test spec (expected geom:exp): $spec"
        printf "%-45s | %-9s | %-7s | %-9s\n" "$spec" "BAD_SPEC" "-" "-" >> "$SUMMARY"
        return
    fi

    local geom_path="$GEOMETRY_DIR/${geom}.opts"
    local exp_path="$EXPERIMENT_DIR/${exp}.opts"
    local test_name="${geom}__${exp}"
    local test_out="$BATCH_DIR/$test_name"

    if [ ! -f "$geom_path" ]; then
        echo "⚠️  Missing geometry: $geom_path"
        printf "%-45s | %-9s | %-7s | %-9s\n" "$test_name" "MISSING_G" "-" "-" >> "$SUMMARY"
        return
    fi
    if [ ! -f "$exp_path" ]; then
        echo "⚠️  Missing experiment: $exp_path"
        printf "%-45s | %-9s | %-7s | %-9s\n" "$test_name" "MISSING_E" "-" "-" >> "$SUMMARY"
        return
    fi

    mkdir -p "$test_out"
    cp "$SOLVER_OPTS" "$test_out/"
    cp "$geom_path"   "$test_out/"
    cp "$exp_path"    "$test_out/"

    local nprocs
    nprocs=$(choose_nprocs "$geom_path")

    echo ""
    echo "=========================================================================="
    echo "  ▶  $geom × $exp  (nprocs=$nprocs)"
    echo "     → $test_out"
    echo "=========================================================================="

    local start end elapsed exit_code status
    start=$(date +%s)

    # The simulation reads its output path from the 'folder' env var (used by
    # monitoring.c) AND from -output_path.  Export both for safety.
    export folder="$test_out"

    set +e
    mpiexec -np "$nprocs" "$EXEC"             \
        -options_file "$SOLVER_OPTS"          \
        -options_file "$geom_path"            \
        -options_file "$exp_path"             \
        -output_path  "$test_out"             \
        2>&1 | tee "$test_out/outp.txt"
    exit_code=${PIPESTATUS[0]}
    set -e

    end=$(date +%s)
    elapsed=$((end - start))

    if [ "$exit_code" -eq 0 ]; then
        status="OK"
        # Sanity check: if the run "succeeded" suspiciously fast, the
        # experiment opts file probably has a tiny t_final left over from
        # debugging. Read t_final from the experiment opts.
        local tfinal
        tfinal=$(awk '$1=="-t_final"{print $2}' "$exp_path" | head -n1)
        if [ -n "${tfinal:-}" ]; then
            # Flag if t_final is unreasonably small (< 60 s of simulated time)
            if awk "BEGIN{exit !($tfinal < 60)}"; then
                status="OK?"
            fi
        fi
    else
        status="FAIL($exit_code)"
    fi

    # Post-processing per dim — read from geometry opts
    local dim
    dim=$(awk '$1 == "-dim" { print $2 }' "$geom_path" | head -n1)
    dim=${dim:-2}

    set +e
    if [[ "$dim" != "1" ]]; then
        # 2D/3D: VTK conversion. Call plotpermafrost.py directly from the
        # project's postprocess/ — earlier this went through
        # run_plotpermafrost.sh, which expects a per-run staged postprocess/
        # subdirectory and silently exits 1 if it doesn't find one.
        if [ -n "$PYTHON" ] && [ -f "$POSTPROCESS/plotpermafrost.py" ]; then
            echo "  Post-processing (VTK)..."
            "$PYTHON" "$POSTPROCESS/plotpermafrost.py" \
                --dir "$test_out" 2>&1 | sed 's/^/    /' || true
        fi
    else
        # 1D: phase profile plots
        if [ -n "$PYTHON" ] && [ -f "$POSTPROCESS/plot1D_profiles.py" ]; then
            echo "  Post-processing (1D profiles)..."
            "$PYTHON" "$POSTPROCESS/plot1D_profiles.py" \
                --dir "$test_out" --out-dir "$test_out" 2>&1 | sed 's/^/    /' || true
        fi
    fi

    # Time-step diagnostic (any dim)
    if [ -n "$PYTHON" ] && [ -f "$POSTPROCESS/plot_timestep.py" ] && [ -f "$test_out/outp.txt" ]; then
        "$PYTHON" "$POSTPROCESS/plot_timestep.py" \
            --dir "$test_out" --save "$test_out/timestep.png" 2>&1 | sed 's/^/    /' || true
    fi

    # Phase mass plot (any dim, if snapshots exist)
    if [ -n "$PYTHON" ] && [ -f "$POSTPROCESS/plot_mass.py" ] && ls "$test_out"/sol_*.dat &>/dev/null; then
        "$PYTHON" "$POSTPROCESS/plot_mass.py" \
            --dir "$test_out" --save "$test_out/mass.png" 2>&1 | sed 's/^/    /' || true
    fi

    # Phase-field free-energy plot (any dim, if snapshots exist)
    if [ -n "$PYTHON" ] && [ -f "$POSTPROCESS/plot_energy.py" ] && ls "$test_out"/sol_*.dat &>/dev/null; then
        "$PYTHON" "$POSTPROCESS/plot_energy.py" \
            --dir "$test_out" --save "$test_out/energy.png" 2>&1 | sed 's/^/    /' || true
    fi
    set -e

    printf "%-45s | %-9s | %-7d | %-9d\n" "$test_name" "$status" "$nprocs" "$elapsed" >> "$SUMMARY"
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
batch_start=$(date +%s)
for opts in "${TESTS[@]}"; do
    run_one_test "$opts"
done
batch_end=$(date +%s)
batch_elapsed=$((batch_end - batch_start))

{
    echo ""
    echo "Finished: $(date)"
    echo "Wall-clock: ${batch_elapsed} s ($(echo "scale=2; $batch_elapsed / 60" | bc) min)"
} >> "$SUMMARY"

echo ""
echo "=========================================================================="
echo "  Batch complete in $(echo "scale=1; $batch_elapsed / 60" | bc) min"
echo "  Summary:"
echo "=========================================================================="
cat "$SUMMARY"
echo ""
echo "Results dir: $BATCH_DIR"
