#!/bin/zsh
###############################################################################
# batch_keff_eval.sh
#
# • Batch processing for computing effective k for many samples
# • Keeps PETSc code unchanged — loop handled in bash
###############################################################################

# ----------------------------
# 🔹  Simulation parameters
# ----------------------------
export Nx=1100
export Ny=1100
export Nz=1          # 1 for 2-D

export Lx=2.0e-3
export Ly=2.0e-3
export Lz=2.02e-4    # ignored when dim=2

export FLUX_BOTTOM=-0.1
export TEMP_TOP=$((273.15-30))

export eps=$((9.09629658751972e-07))
export dim=2         # 2 = 2-D, 3 = 3-D

INPUT_DIR="/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/NASAv2_118G_2D_T-20.0_hum0.70_2025-05-31__09.50.53"
OUTPUT_BASE="/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/homog"
timestamp=$(date +%Y-%m-%d__%H.%M.%S)
OUT_FOLDER="ThermalSim_homog_$timestamp"
export OUTPUT_DIR="$OUTPUT_BASE/$OUT_FOLDER"
mkdir -p "$OUTPUT_DIR"

SUMMARY_FILE="$OUTPUT_DIR/keff_summary.csv"
echo "filename,k_eff_xx,k_eff_yy,k_eff_xy,k_eff_yx" > "$SUMMARY_FILE"

export OUTPUT_VTK=0
export OUTPUT_BINARY=1

NUM_PROCS=${NUM_PROCS:-1}

# ----------------------------
# 🔹  Helpers
# ----------------------------
compile_code() {
    echo "🔨 Compiling (release)…"
    make BUILD=release effective_k_ice_homog
}

run_simulation() {
    echo "🏃 Running with $NUM_PROCS MPI proc(s)…"
    mpiexec -np $NUM_PROCS ./effective_k_ice_homog \
        -init_mode "$INIT_MODE" -init_dir "$INIT_DIR"
}

collect_outputs() {
    echo "📂 Moving output files to $1"
    mv *.dat  "$1" 2>/dev/null
    mv *.bin  "$1" 2>/dev/null
    mv *.info "$1" 2>/dev/null
}

# ----------------------------
# 🔹  Main Loop
# ----------------------------
compile_code || { echo "❌ compile failed"; exit 1; }

for SOL in "$INPUT_DIR"/sol_*.dat; do
    BASENAME=$(basename "$SOL" .dat)

    echo "🔁 Processing $BASENAME"
    INIT_MODE="$SOL"
    INIT_DIR="$INPUT_DIR"

    run_simulation || { echo "❌ simulation failed for $BASENAME"; continue; }
    collect_outputs "$OUTPUT_DIR"

    echo "📈 Post-processing…"
    # Assuming your PETSc code already prints k_eff_xx and k_eff_yy to stdout
    # You may want to modify your C code to **also** print k_eff_xy and k_eff_yx

    KLINE=$(python3 postprocess/plot_vector_field.py "$OUT_FOLDER" "$Nx" "$Ny" "$Lx" "$Ly" | tee /dev/tty | grep "Estimated k_eff")
    # Extract k_eff components
    KXX=$(echo "$KLINE" | sed -n 's/.*k_eff_xx = \([^,]*\),.*/\1/p')
    KYY=$(echo "$KLINE" | sed -n 's/.*k_eff_yy = \([^ ]*\).*/\1/p')
    # Placeholders for off-diagonal terms (if you have them)
    KXY="0.0"
    KYX="0.0"

    echo "$BASENAME,$KXX,$KYY,$KXY,$KYX" >> "$SUMMARY_FILE"
done

echo "✅ Batch processing complete. Results saved to $SUMMARY_FILE"