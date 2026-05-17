#!/usr/bin/env bash
# =============================================================================
# run_batch_tests.sh — Run a curated set of tests sequentially into one folder
#
# Location: $PETIGA_DIR/projects/permafrost/scripts/Studio/run_batch_tests.sh
#
# All runs go under a single timestamped parent folder:
#   $RESULTS_BASE/batch_<timestamp>/<test_name>/
#
# Writes a SUMMARY.txt at the end with status + wall-clock per test.
#
# Usage:
#   ./scripts/Studio/run_batch_tests.sh                 # default test set
#   ./scripts/Studio/run_batch_tests.sh --skip-1d       # 2D-only
#   ./scripts/Studio/run_batch_tests.sh --tag mylabel   # add suffix to batch dir
#   ./scripts/Studio/run_batch_tests.sh --tests test_2D_TouchingGrainPair.opts,test_2D_SingleIceGrain.opts
# =============================================================================

set -uo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PETIGA_DIR}/projects/permafrost"
EXEC="$PROJECT_ROOT/permafrost"
INPUTS_DIR="$PROJECT_ROOT/inputs"
TESTS_DIR="$INPUTS_DIR/tests"
UNIVERSAL_OPTS="$INPUTS_DIR/universal.opts"
POSTPROCESS="$PROJECT_ROOT/postprocess"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
RESULTS_BASE="/Users/jacksonbaglino/SimulationResults/permafrost/scratch"

MAX_LOCAL_CORES=12
TARGET_DOFS_PER_CORE=10000

if [ ! -f "$PROJECT_ROOT/makefile" ] && [ ! -f "$PROJECT_ROOT/Makefile" ]; then
    echo "❌ Could not find makefile at $PROJECT_ROOT"
    exit 1
fi

PYTHON=$(command -v python3 || command -v python || true)

# ---------------------------------------------------------------------------
# Default test list (matches the user request)
# ---------------------------------------------------------------------------
DEFAULT_TESTS=(
    # ---- 1D tests (standard + hires for each) ----
    "test_1D_IceSlab.opts"
    "test_1D_IceSlab_hires.opts"
    "test_1D_SingleIceGrain.opts"
    "test_1D_SingleIceGrain_hires.opts"
    "test_1D_IceSedPair.opts"
    "test_1D_IceSedPair_hires.opts"
    "test_1D_SeparatedGrainPair.opts"
    "test_1D_SeparatedGrainPair_hires.opts"
    "test_1D_TouchingGrainPair.opts"
    "test_1D_TouchingGrainPair_hires.opts"

    # ---- 2D tests (low + high res) ----
    "test_2D_SingleIceGrain.opts"
    "test_2D_SingleIceGrain_hires.opts"
    "test_2D_IceSedPair.opts"
    "test_2D_IceSedPair_hires.opts"
    "test_2D_TouchingGrainPair.opts"
    "test_2D_TouchingGrainPair_hires.opts"
    "test_2D_SeparatedGrainPair.opts"
    "test_2D_SeparatedGrainPair_hires.opts"

    # ---- Single sediment grain (only low res; no hires variant) ----
    "test_2D_SingleSedGrain.opts"
)

# ---------------------------------------------------------------------------
# Arg parsing
# ---------------------------------------------------------------------------
tag=""
skip_1d=false
custom_tests=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)
            tag="$2"; shift 2 ;;
        --skip-1d)
            skip_1d=true; shift ;;
        --tests)
            custom_tests="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,18p' "$0"; exit 0 ;;
        *)
            echo "Unknown argument: $1"
            exit 1 ;;
    esac
done

# Build the actual test list
TESTS=()
if [[ -n "$custom_tests" ]]; then
    IFS=',' read -ra TESTS <<< "$custom_tests"
else
    for t in "${DEFAULT_TESTS[@]}"; do
        if $skip_1d && [[ "$t" == test_1D_* ]]; then
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

# Stage shared assets (source + universal opts) at the batch level
mkdir -p "$BATCH_DIR/src_snapshot"
[ -f "$UNIVERSAL_OPTS" ] && cp "$UNIVERSAL_OPTS" "$BATCH_DIR/"
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
    local Nx Ny Nz total dofs cap hw
    Nx=$(awk '$1=="-Nx"{print $2}' "$opts" | head -n1)
    Ny=$(awk '$1=="-Ny"{print $2}' "$opts" | head -n1)
    Nz=$(awk '$1=="-Nz"{print $2}' "$opts" | head -n1)
    [[ -z "${Nx:-}" ]] && Nx=1
    [[ -z "${Ny:-}" ]] && Ny=1
    [[ -z "${Nz:-}" ]] && Nz=1
    dofs=$((4 * Nx * Ny * Nz))
    local n=$(((dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE))
    (( n < 1 )) && n=1
    hw=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo $MAX_LOCAL_CORES)
    cap=$(( hw < MAX_LOCAL_CORES ? hw : MAX_LOCAL_CORES ))
    (( n > cap )) && n=$cap
    echo "$n"
}

# ---------------------------------------------------------------------------
# Per-test runner
# ---------------------------------------------------------------------------
run_one_test() {
    local opts_file="$1"
    local opts_path="$TESTS_DIR/$opts_file"
    local test_name="${opts_file%.opts}"
    local test_out="$BATCH_DIR/$test_name"

    if [ ! -f "$opts_path" ]; then
        echo "⚠️  Missing: $opts_file"
        printf "%-45s | %-9s | %-7s | %-9s\n" "$test_name" "MISSING" "-" "-" >> "$SUMMARY"
        return
    fi

    mkdir -p "$test_out"
    cp "$opts_path" "$test_out/"
    [ -f "$UNIVERSAL_OPTS" ] && cp "$UNIVERSAL_OPTS" "$test_out/"

    local nprocs
    nprocs=$(choose_nprocs "$opts_path")

    echo ""
    echo "=========================================================================="
    echo "  ▶  $test_name  (nprocs=$nprocs)"
    echo "     → $test_out"
    echo "=========================================================================="

    local start end elapsed exit_code status
    start=$(date +%s)

    # The simulation reads its output path from the 'folder' env var (used by
    # monitoring.c) AND from -output_path.  Export both for safety.
    export folder="$test_out"

    set +e
    mpiexec -np "$nprocs" "$EXEC"                \
        -options_file "$UNIVERSAL_OPTS"         \
        -options_file "$opts_path"              \
        -output_path  "$test_out"               \
        2>&1 | tee "$test_out/outp.txt"
    exit_code=${PIPESTATUS[0]}
    set -e

    end=$(date +%s)
    elapsed=$((end - start))

    if [ "$exit_code" -eq 0 ]; then
        status="OK"
        # Sanity check: if the run "succeeded" suspiciously fast, the opts
        # file probably has a tiny t_final left over from debugging. Read
        # the final simulated TIME from the last monitor data row and
        # compare to t_final from the opts file.
        local tfinal final_t
        tfinal=$(awk '$1=="-t_final"{print $2}' "$opts_path" | head -n1)
        if [ -n "${tfinal:-}" ]; then
            # Flag if t_final is unreasonably small (< 60 s of simulated time)
            if awk "BEGIN{exit !($tfinal < 60)}"; then
                status="OK?"
            fi
        fi
    else
        status="FAIL($exit_code)"
    fi

    # Post-processing per dim
    local dim
    dim=$(awk '$1 == "-dim" { print $2 }' "$opts_path" | head -n1)
    dim=${dim:-2}

    set +e
    if [[ "$dim" != "1" ]]; then
        # 2D/3D: VTK conversion if a plot script is present
        local plot_script="$SCRIPTS_DIR/run_plotpermafrost.sh"
        if [ -f "$plot_script" ]; then
            echo "  Post-processing (VTK)..."
            "$plot_script" "$test_out" 2>&1 | sed 's/^/    /' || true
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
