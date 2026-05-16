#!/usr/bin/env bash
# =============================================================================
# run_permafrost.sh — Build, run, and post-process a permafrost simulation
#
# Location: $PETIGA_DIR/projects/permafrost/scripts/Studio/run_permafrost.sh
#
# Usage:
#   ./scripts/Studio/run_permafrost.sh <PETSc_options_file> [tag]
#
# The output folder is auto-derived from the opts file's ic_type and basename:
#   $RESULTS_BASE/<ic_category>/<opts_basename>[_tag]_<timestamp>/
#
# Example:
#   ./scripts/Studio/run_permafrost.sh inputs/tests/test_2D_TouchingGrainPair.opts
#   ./scripts/Studio/run_permafrost.sh inputs/tests/test_1D_IceSlab.opts sweep_a
# =============================================================================

set -uo pipefail

# ---------------------------------------------------------------------------
# Global state
# ---------------------------------------------------------------------------
folder=""
name=""
sim_exit=0
NPROCS=1
MAX_LOCAL_CORES=12   # physical cores on this Mac

# ---------------------------------------------------------------------------
# Resolve project root — always two levels above this script's location
# Script lives at: <project_root>/scripts/Studio/run_permafrost.sh
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PROJECT_ROOT="${PETIGA_DIR}/projects/permafrost"

# Validate that we resolved the project root correctly
if [ ! -f "$PROJECT_ROOT/makefile" ] && [ ! -f "$PROJECT_ROOT/Makefile" ]; then
    echo "❌ Error: Could not find makefile in resolved project root: $PROJECT_ROOT"
    echo "   Script location: $SCRIPT_DIR"
    echo "   Expected project root two levels up from script."
    exit 1
fi

# Convenience aliases
EXEC="$PROJECT_ROOT/permafrost"
SRC_DIR="$PROJECT_ROOT/src"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
INPUTS_DIR="$PROJECT_ROOT/inputs"
RESULTS_BASE="/Users/jacksonbaglino/SimulationResults/permafrost/scratch"
UNIVERSAL_OPTS="$INPUTS_DIR/universal.opts"

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    echo ""
    echo "Usage:"
    echo "  ./scripts/Studio/run_permafrost.sh <PETSc_options_file> [tag]"
    echo ""
    echo "Arguments:"
    echo "  PETSc_options_file   Path to PETSc .opts file (relative to project root)"
    echo "  tag                  Optional label appended to the run folder name"
    echo ""
    echo "The output folder is derived automatically from the opts file:"
    echo "  \$RESULTS_BASE/<ic_type_category>/<opts_basename>[_tag]_<timestamp>/"
    echo ""
    echo "Example:"
    echo "  ./scripts/Studio/run_permafrost.sh ./inputs/tests/test_2D_TouchingGrainPair.opts"
    echo "  ./scripts/Studio/run_permafrost.sh ./inputs/tests/test_1D_IceSlab.opts sweep_a"
    echo ""
}

if [ "$#" -lt 1 ]; then
    echo "❌ Error: Missing required arguments."
    usage
    exit 1
fi

params_file="$1"
title="${2:-}"

# Resolve params_file relative to project root if not absolute
if [[ "$params_file" != /* ]]; then
    params_file="$PROJECT_ROOT/$params_file"
fi

if [ ! -f "$params_file" ]; then
    echo "❌ Error: Options file not found: $params_file"
    exit 1
fi

trap 'echo "❌ Script error on line $LINENO"' ERR

# ---------------------------------------------------------------------------
# Simulation time parameters
# ---------------------------------------------------------------------------
t_final=$(echo "60 * 24 * 3600" | bc -l)   # 60 days in seconds
n_out=100

if (( $(echo "$t_final <= 0" | bc -l) )); then
    echo "❌ Error: t_final must be > 0. Got $t_final"
    exit 1
fi

# ---------------------------------------------------------------------------
# compute_optimal_nprocs
# Reads Nx/Ny/Nz from the opts file, targets ~10 000 DOFs/core, and caps at
# MAX_LOCAL_CORES (12).  Sets the global NPROCS.
# ---------------------------------------------------------------------------
compute_optimal_nprocs() {
    local Nx Ny Nz total_dofs
    local TARGET_DOFS_PER_CORE=10000

    Nx=$(awk '$1=="-Nx"{print $2}' "$params_file" | head -n1)
    Ny=$(awk '$1=="-Ny"{print $2}' "$params_file" | head -n1)
    Nz=$(awk '$1=="-Nz"{print $2}' "$params_file" | head -n1)

    [[ -z "${Nx:-}" ]] && Nx=1
    [[ -z "${Ny:-}" ]] && Ny=1
    [[ -z "${Nz:-}" ]] && Nz=1

    total_dofs=$((4 * Nx * Ny * Nz))
    NPROCS=$(((total_dofs + TARGET_DOFS_PER_CORE - 1) / TARGET_DOFS_PER_CORE))
    (( NPROCS < 1 )) && NPROCS=1

    # On macOS use sysctl; fall back to nproc for Linux
    local hw_cores
    hw_cores=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo $MAX_LOCAL_CORES)
    local cap=$(( hw_cores < MAX_LOCAL_CORES ? hw_cores : MAX_LOCAL_CORES ))
    (( NPROCS > cap )) && NPROCS=$cap

    echo "------------------------------------------------------------"
    echo "Grid from opts: Nx=${Nx}, Ny=${Ny}, Nz=${Nz}"
    echo "Total DoFs:     4 × Nx × Ny × Nz = ${total_dofs}"
    echo "Target/core:    ${TARGET_DOFS_PER_CORE} DoFs"
    echo "Chosen NPROCS:  ${NPROCS}  (cap: ${cap} logical cores)"
    echo "------------------------------------------------------------"
}

# ---------------------------------------------------------------------------
# compile_code
# Runs make from the project root where the makefile lives
# ---------------------------------------------------------------------------
compile_code() {
    echo ""
    echo "--- Compiling ---"
    echo "Project root: $PROJECT_ROOT"

    cd "$PROJECT_ROOT"
    make

    if [ ! -f "$EXEC" ]; then
        echo "❌ Executable not found after compilation: $EXEC"
        exit 1
    fi
    echo "✅ Compilation successful."
}

# ---------------------------------------------------------------------------
# derive_ic_subfolder
# Reads -ic_type from the opts file and maps it to a clean category name
# ---------------------------------------------------------------------------
derive_ic_subfolder() {
    local ic
    ic=$(awk '$1=="-ic_type"{print $2}' "$params_file" | head -n1)
    case "${ic:-}" in
        ice_slab)        echo "IceSlab" ;;
        enclosed)        echo "EnclosedGrainPair" ;;
        contact_sed)     echo "ContactSed" ;;
        capillary)       echo "CapillaryBridge" ;;
        slab_and_grains) echo "SlabAndGrains" ;;
        ice_cap)         echo "IceCap" ;;
        single_ice)      echo "SingleIceGrain" ;;
        ice_sed_pair)    echo "IceSedPair" ;;
        *)               echo "Other" ;;
    esac
}

# ---------------------------------------------------------------------------
# create_folder
# Creates an output directory under $RESULTS_BASE/<ic_category>/<opts_name>[_tag]_<ts>
# ---------------------------------------------------------------------------
create_folder() {
    echo ""
    echo "--- Creating output folder ---"
    local subfolder ts opts_name tag
    subfolder=$(derive_ic_subfolder)
    opts_name=$(basename "$params_file" .opts)
    ts=$(date +%Y-%m-%d__%H.%M.%S)
    tag="${title:+_${title}}"
    name="${opts_name}${tag}_${ts}"
    folder="$RESULTS_BASE/$subfolder/$name"

    mkdir -p "$folder"
    echo "Output folder: $folder"
}

# ---------------------------------------------------------------------------
# stage_output_folder
# Copies the options file and helper scripts into the output folder before
# the run so the full configuration is captured alongside results
# ---------------------------------------------------------------------------
stage_output_folder() {
    echo ""
    echo "--- Staging output folder ---"

    # Copy the universal opts and the simulation-specific opts used for this run
    [ -f "$UNIVERSAL_OPTS" ] && cp "$UNIVERSAL_OPTS" "$folder/"
    cp "$params_file" "$folder/$(basename "$params_file")"

    # Post-processing scripts — copy the full postprocess/ directory
    local POSTPROCESS="$PROJECT_ROOT/postprocess"
    if [ -d "$POSTPROCESS" ]; then
        cp -r "$POSTPROCESS" "$folder/postprocess"
        echo "  Copied postprocess/ → $folder/postprocess/"
    else
        echo "⚠️  postprocess/ directory not found at $POSTPROCESS — skipping."
    fi

    # Copy this run script itself
    cp "${BASH_SOURCE[0]}" "$folder/run_permafrost.sh"

    echo "✅ Staging complete."
}

# ---------------------------------------------------------------------------
# copy_source_code
# Saves source files from src/ for exact reproducibility
# ---------------------------------------------------------------------------
copy_source_code() {
    echo ""
    echo "--- Copying source code ---"

    local src_dest="$folder/src"
    mkdir -p "$src_dest"

    local found=0

    # Copy all C source and header files from src/
    if [ -d "$SRC_DIR" ]; then
        for ext in "*.c" "*.h" "*.cpp" "*.hpp"; do
            for f in "$SRC_DIR"/$ext; do
                if [ -f "$f" ]; then
                    cp "$f" "$src_dest/"
                    found=1
                fi
            done
        done
    else
        echo "⚠️  Warning: src/ directory not found at $SRC_DIR"
    fi

    # Copy include/ directory if present
    if [ -d "$PROJECT_ROOT/include" ]; then
        cp -r "$PROJECT_ROOT/include" "$src_dest/"
        found=1
    fi

    # Copy the makefile
    for mf in makefile Makefile; do
        if [ -f "$PROJECT_ROOT/$mf" ]; then
            cp "$PROJECT_ROOT/$mf" "$src_dest/"
            found=1
            break
        fi
    done

    if [ "$found" -eq 0 ]; then
        echo "⚠️  Warning: No source files found."
    else
        echo "✅ Source code copied to $src_dest"
    fi
}

# ---------------------------------------------------------------------------
# run_simulation
# Executes the MPI simulation and logs all output
# ---------------------------------------------------------------------------
run_simulation() {
    echo ""
    echo "--- Running simulation ---"
    echo "Executable   : $EXEC"
    echo "Universal    : ${UNIVERSAL_OPTS}"
    echo "Options file : $params_file"
    echo "Output path  : $folder"
    echo "Processes    : $NPROCS"
    echo ""

    export folder

    local universal_arg=""
    [ -f "$UNIVERSAL_OPTS" ] && universal_arg="-options_file $UNIVERSAL_OPTS"

    set +e
    mpiexec -np "$NPROCS" "$EXEC"       \
        $universal_arg                  \
        -options_file "$params_file"    \
        -output_path  "$folder"         \
        | tee "$folder/outp.txt"

    sim_exit=${PIPESTATUS[0]}
    set -e

    if [ "$sim_exit" -ne 0 ]; then
        echo "⚠️  Simulation exited with code $sim_exit (continuing to post-processing)"
    else
        echo "✅ Simulation completed successfully."
    fi
}

# ---------------------------------------------------------------------------
# run_plotting
# Calls the post-processing plotting script if it exists
# ---------------------------------------------------------------------------
run_plotting() {
    # Skip VTK for 1D runs — not meaningful and the file format differs
    local dim
    dim=$(awk '$1 == "-dim" { print $2 }' "$params_file" | head -n1)
    dim=${dim:-2}
    if [[ "$dim" == "1" ]]; then
        return
    fi

    echo ""
    echo "--- Running post-processing (VTK) ---"

    local plot_script="$SCRIPTS_DIR/run_plotpermafrost.sh"

    if [ ! -f "$plot_script" ]; then
        echo "⚠️  Plotting script not found: $plot_script — skipping."
        return
    fi

    set +e
    # Pass the full output folder path so the script works regardless of subfolder depth
    "$plot_script" "$folder"
    local plot_exit=$?
    set -e

    if [ "$plot_exit" -ne 0 ]; then
        echo "⚠️  VTK plotting exited with code $plot_exit"
    else
        echo "✅ VTK post-processing complete."
    fi
}

# ---------------------------------------------------------------------------
# run_1d_plotting
# Detects 1D runs from the opts file and invokes the 1D Python post-processing
# scripts automatically.  Silently skips for 2D/3D runs.
# ---------------------------------------------------------------------------
run_1d_plotting() {
    # Extract the -dim value from the options file (default 2 if absent)
    local dim
    dim=$(awk '$1 == "-dim" { print $2 }' "$params_file" | head -n1)
    dim=${dim:-2}

    if [[ "$dim" != "1" ]]; then
        return   # not a 1D run — nothing to do
    fi

    echo ""
    echo "--- 1D post-processing ---"

    local POSTPROCESS="$PROJECT_ROOT/postprocess"
    local py_exit=0

    if ! command -v python &>/dev/null && ! command -v python3 &>/dev/null; then
        echo "⚠️  python not found — skipping 1D plots."
        return
    fi

    local PYTHON
    PYTHON=$(command -v python3 || command -v python)

    set +e

    # Per-step phase field PNGs
    echo "  Generating phase field images..."
    "$PYTHON" "$POSTPROCESS/plot1D_profiles.py" \
        --dir "$folder" --out-dir "$folder" \
        2>&1 | sed 's/^/    /'
    py_exit=$(( py_exit + $? ))

    # Derived scalar time-series
    echo "  Generating derived quantity plot..."
    "$PYTHON" "$POSTPROCESS/plot1D_profiles.py" \
        --dir "$folder" --derived --save "$folder/derived.png" \
        2>&1 | sed 's/^/    /'
    py_exit=$(( py_exit + $? ))

    # SSA scalar time-series
    if [ -f "$folder/SSA_evo.dat" ]; then
        echo "  Generating SSA scalar plot..."
        "$PYTHON" "$POSTPROCESS/plot_scalars.py" \
            --file "$folder/SSA_evo.dat" --save "$folder/scalars.png" \
            2>&1 | sed 's/^/    /'
        py_exit=$(( py_exit + $? ))
    fi

    set -e

    if [ "$py_exit" -ne 0 ]; then
        echo "⚠️  One or more 1D plots failed (exit sum $py_exit) — check output above."
    else
        echo "✅ 1D post-processing complete."
        echo "   phase_step_*.png  →  $folder/"
        echo "   derived.png       →  $folder/derived.png"
        echo "   scalars.png          →  $folder/scalars.png"
    fi
}

# =============================================================================
# Main workflow
# =============================================================================
echo ""
echo "========================================================================="
echo "  Permafrost simulation workflow"
echo "  Project root : $PROJECT_ROOT"
echo "  Options      : $params_file"
echo "  Title        : $title"
echo "========================================================================="

compute_optimal_nprocs
compile_code
create_folder
stage_output_folder
run_simulation
copy_source_code
run_plotting
run_1d_plotting

PYTHON=$(command -v python3 || command -v python)

# Time step diagnostic — always generated from outp.txt (all dims)
if [ -f "$folder/outp.txt" ]; then
    echo ""
    echo "--- Time step diagnostic ---"
    "$PYTHON" "$PROJECT_ROOT/postprocess/plot_timestep.py" \
        --dir "$folder" --save "$folder/timestep.png" \
        2>&1 | sed 's/^/    /'
fi

# Phase mass plot — always generated when solution snapshots are present
if ls "$folder"/sol_*.dat &>/dev/null; then
    echo ""
    echo "--- Phase mass plot ---"
    "$PYTHON" "$PROJECT_ROOT/postprocess/plot_mass.py" \
        --dir "$folder" --save "$folder/mass.png" \
        2>&1 | sed 's/^/    /'
fi

echo ""
echo "========================================================================="
if [ "${sim_exit:-0}" -ne 0 ]; then
    echo "  ⚠️  Workflow completed with simulation errors (exit code $sim_exit)"
else
    echo "  ✅ Workflow completed successfully"
fi
echo "  Results: $folder"
echo "========================================================================="
echo ""

exit "${sim_exit:-0}"