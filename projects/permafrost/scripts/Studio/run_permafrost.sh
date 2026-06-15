#!/usr/bin/env bash
# =============================================================================
# run_permafrost.sh — Build, run, and post-process a permafrost simulation
#
# Location: $PETIGA_DIR/projects/permafrost/scripts/Studio/run_permafrost.sh
#
# Usage:
#   ./scripts/Studio/run_permafrost.sh <geometry> <experiment> [tag]
#
#   geometry    Name (without .opts) of a file in inputs/geometry/. Sets up
#               the domain/mesh/IC: -dim, -Lx/-Ly/-Lz, -Nx/-Ny/-Nz (or
#               -geom_file for an igakit-generated mesh), -ic_type and its
#               grain/geometry parameters, -eps, -delt_t, etc.
#   experiment  Name (without .opts) of a file in inputs/experiment/. Sets
#               the run conditions: -t_final, -temp, -humidity,
#               -grad_temp0, and (optionally) output cadence (-outp/-n_out).
#   tag         Optional label appended to the run folder name.
#
# Combining the two: the script concatenates solver.opts + the chosen
# geometry .opts + the chosen experiment .opts (later files can override
# earlier settings) and passes the result to `permafrost`.
#
# The output folder is auto-derived from the geometry opts' -ic_type and the
# two opts basenames:
#   $RESULTS_BASE/<ic_type_category>/<geometry>__<experiment>[_<tag>]_<timestamp>/
#
# Examples:
#   ./scripts/Studio/run_permafrost.sh 2D_touching_grains 1day_T-20_h1.00
#   ./scripts/Studio/run_permafrost.sh 1D_ice_slab 1day_T-20_h0.95 sweep_a
#   ./scripts/Studio/run_permafrost.sh 2D_multi_grain_test 2day_T-20_h0.95 multigrain_2day
#
# -----------------------------------------------------------------------
# Adding a new EXPERIMENT (run conditions)
# -----------------------------------------------------------------------
# Copy an existing file, e.g. inputs/experiment/2day_T-20_h0.95.opts, into
# inputs/experiment/<my_experiment>.opts and edit -t_final/-temp/-humidity/
# -grad_temp0 (and -outp/-n_out if you want non-default output cadence —
# solver.opts defaults to -outp 1, i.e. a snapshot every time step). Then:
#   ./scripts/Studio/run_permafrost.sh <geometry> my_experiment [tag]
#
# -----------------------------------------------------------------------
# Adding a new GEOMETRY / initial condition
# -----------------------------------------------------------------------
# Option A — rectangular domain, existing -ic_type:
#   Copy an existing file, e.g. inputs/geometry/2D_two_ice_grains_boundary.opts,
#   into inputs/geometry/<my_geometry>.opts and edit -dim, -Lx/-Ly/-Lz,
#   -Nx/-Ny/-Nz, -ic_type and its parameters (e.g. -ice_grain_cx/cy/R for
#   -ic_type multi_grains), -eps, -delt_t.
#
# Option B — custom (non-rectangular) mesh via igakit:
#   1. Write/adapt a preprocess/build_geometry_*.py script (see
#      preprocess/build_geometry_multi_grain.py for a template) that writes
#      an IGA mesh file under inputs/geometry/<name>.dat. Run it with:
#        source venv_pf311/bin/activate && python3 preprocess/build_geometry_<name>.py
#      (run from the project root — paths inside these scripts are relative
#      to the project root, not preprocess/).
#   2. Create inputs/geometry/<my_geometry>.opts with -geom_file pointing at
#      that .dat file (this overrides -p/-C/-Nx/-Ny/-Nz), plus -ic_type and
#      its parameters, -dim, -Lx/-Ly/-Lz, -eps, -delt_t, -periodic 0.
#
# Then run as usual:
#   ./scripts/Studio/run_permafrost.sh my_geometry <experiment> [tag]
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
SOLVER_OPTS="$INPUTS_DIR/solver.opts"
GEOMETRY_DIR="$INPUTS_DIR/geometry"
EXPERIMENT_DIR="$INPUTS_DIR/experiment"

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    echo ""
    echo "Usage:"
    echo "  ./scripts/Studio/run_permafrost.sh <geometry> <experiment> [tag]"
    echo ""
    echo "Arguments:"
    echo "  geometry    Name (without .opts) of a file in inputs/geometry/"
    echo "              e.g. 2D_touching_grains, 1D_ice_slab"
    echo "  experiment  Name (without .opts) of a file in inputs/experiment/"
    echo "              e.g. 1day_T-20_h1.00"
    echo "  tag         Optional label appended to the run folder name"
    echo ""
    echo "The output folder is:"
    echo "  \$RESULTS_BASE/<ic_type_category>/<geom>__<exp>[_<tag>]_<timestamp>/"
    echo ""
    echo "Example:"
    echo "  ./scripts/Studio/run_permafrost.sh 2D_touching_grains 1day_T-20_h1.00"
    echo "  ./scripts/Studio/run_permafrost.sh 1D_ice_slab 1day_T-20_h0.95 sweep_a"
    echo ""
}

if [ "$#" -lt 2 ]; then
    echo "❌ Error: Missing required arguments."
    usage
    exit 1
fi

geom_name="$1"
exp_name="$2"
title="${3:-}"

GEOM_OPTS="$GEOMETRY_DIR/${geom_name}.opts"
EXP_OPTS="$EXPERIMENT_DIR/${exp_name}.opts"

if [ ! -f "$SOLVER_OPTS" ]; then
    echo "❌ Error: Solver opts not found: $SOLVER_OPTS"
    exit 1
fi
if [ ! -f "$GEOM_OPTS" ]; then
    echo "❌ Error: Geometry opts not found: $GEOM_OPTS"
    echo "   Available geometries:"
    ls "$GEOMETRY_DIR" 2>/dev/null | sed 's/\.opts$//' | sed 's/^/     /'
    exit 1
fi
if [ ! -f "$EXP_OPTS" ]; then
    echo "❌ Error: Experiment opts not found: $EXP_OPTS"
    echo "   Available experiments:"
    ls "$EXPERIMENT_DIR" 2>/dev/null | sed 's/\.opts$//' | sed 's/^/     /'
    exit 1
fi

# For backward-compat with the rest of the script, treat the geometry opts
# as the "params file" — Nx/Ny/Nz and ic_type live there.
params_file="$GEOM_OPTS"

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
        single_sed)      echo "SingleSedGrain" ;;
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
    local subfolder ts tag
    subfolder=$(derive_ic_subfolder)
    ts=$(date +%Y-%m-%d__%H.%M.%S)
    tag="${title:+_${title}}"
    name="${geom_name}__${exp_name}${tag}_${ts}"
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

    # Copy all three opts used by this run
    cp "$SOLVER_OPTS" "$folder/"
    cp "$GEOM_OPTS"   "$folder/"
    cp "$EXP_OPTS"    "$folder/"

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
    echo "Solver       : $SOLVER_OPTS"
    echo "Geometry     : $GEOM_OPTS"
    echo "Experiment   : $EXP_OPTS"
    echo "Output path  : $folder"
    echo "Processes    : $NPROCS"
    echo ""

    export folder

    set +e
    mpiexec -np "$NPROCS" "$EXEC"      \
        -options_file "$SOLVER_OPTS"   \
        -options_file "$GEOM_OPTS"     \
        -options_file "$EXP_OPTS"      \
        -output_path  "$folder"        \
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
echo "  Geometry     : $geom_name  ($GEOM_OPTS)"
echo "  Experiment   : $exp_name  ($EXP_OPTS)"
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