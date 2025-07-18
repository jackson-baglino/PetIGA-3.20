#!/bin/zsh

start_time=$(date +%s)

# ----------------------------
# 🔹  Helper function to read JSON values
# ----------------------------
get_json_value() {
  local key=$1
  jq -r ".$key" "$INIT_DIR/metadata.json"
}

get_nested_json_value() {
  local section=$1
  local key=$2
  jq -r ".$section.$key" "$INIT_DIR/metadata.json"
}

# ----------------------------
# 🔹  Input: INIT_DIR and derived metadata
# ----------------------------
INIT_DIR="/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/scratch/drysnow_2G_Molaro_tight_2D_Tm20.0_hum100_tf28d__2025-06-24__09.31.53"
base_folder=$(basename "$INIT_DIR")

metadata_file="$INIT_DIR/metadata.json"
if [[ ! -f "$metadata_file" ]]; then
  echo "❌ No metadata.json file found in $INIT_DIR"
  exit 1
fi

# Extract from metadata.json
grains="$(get_json_value "number_of_grains")"
dims="$(get_json_value "sim_dimension")"
temp="$(get_json_value "temperature_C")"
hum="$(get_json_value "humidity")"

Lx="$(get_nested_json_value "domain_size_m" "Lx")"
Ly="$(get_nested_json_value "domain_size_m" "Ly")"
Lz="$(get_nested_json_value "domain_size_m" "Lz")"

Nx="$(get_nested_json_value "mesh_resolution" "Nx")"
Ny="$(get_nested_json_value "mesh_resolution" "Ny")"

Nz=1  # still fixed for 2D

eps="$(get_json_value "interface_width_eps")"

echo "Parameters from metadata.json:"
echo "  Grains: $grains"
echo "  Dimensions: $dims"
echo "  Temperature: $temp °C"
echo "  Humidity: $hum"
echo "  Domain size: Lx=$Lx, Ly=$Ly, Lz=$Lz"
echo "  Mesh resolution: Nx=$Nx, Ny=$Ny, Nz=$Nz"
# ----------------------------
# 🔹  Fixed or derived parameters
# ----------------------------
dim=$dims

INIT_MODE="FILE"

timestamp=$(date +%Y-%m-%d__%H.%M.%S)
OUT_FOLDER="ThermalSim_homog_${grains}G_${dims}D_T${temp}_hum${hum}_$timestamp"
OUTPUT_DIR="/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch/Test/$OUT_FOLDER"
mkdir -p "$OUTPUT_DIR"

echo "✅ Using INIT_DIR: $INIT_DIR"
echo "✅ OUT_FOLDER: $OUT_FOLDER"
echo "📏 Domain size: Lx=$Lx, Ly=$Ly"
echo "📊 Mesh size: Nx=$Nx, Ny=$Ny"
echo "🌡️  Temp = $temp °C, Humidity = $hum, Grains = $grains"

# Output flags
export OUTPUT_VTK=1
export OUTPUT_BINARY=1
export SOL_INDEX=-1

export Nx Ny Nz Lx Ly Lz eps dim
export INIT_MODE INIT_DIR OUTPUT_DIR

# MPI ranks
NUM_PROCS=${NUM_PROCS:-1}

# ----------------------------
# 🔹  Helpers
# ----------------------------
compile_code() {
  echo "🔨 Compiling (release)…"
  make BUILD=release effective_k_ice_homog
}

run_simulation() {
  echo "🚀 Running with $NUM_PROCS MPI proc(s)…"
  mpiexec -np $NUM_PROCS ./effective_k_ice_homog \
    -init_mode "$INIT_MODE" \
    -init_dir "$INIT_DIR"
}

collect_outputs() {
  echo "📂 Moving output files to $OUTPUT_DIR"
  mv *.dat  "$OUTPUT_DIR" 2>/dev/null
  mv *.bin  "$OUTPUT_DIR" 2>/dev/null
  mv *.info "$OUTPUT_DIR" 2>/dev/null
  mv *.csv  "$OUTPUT_DIR" 2>/dev/null
}

# ----------------------------
# 🔹  Main workflow
# ----------------------------
compile_code   || { echo "❌ compile failed"; exit 1; }
run_simulation || { echo "❌ simulation failed"; exit 1; }
collect_outputs

echo "📈 Post-processing…"
python3 postprocess/plot_vector_field.py "$OUT_FOLDER" "$Nx" "$Ny" "$Lx" "$Ly"

SSA_FILE="$INIT_DIR/SSA_evo.dat"
cp "$SSA_FILE" "$OUTPUT_DIR/SSA_evo.dat" 2>/dev/null || echo "SSA file not found, skipping copy."
cp "$INIT_DIR"/grainReadFile*.dat "$OUTPUT_DIR/" 2>/dev/null || echo "grainReadFile not found, skipping copy."
cp "$INIT_DIR"/metadata.json "$OUTPUT_DIR/" 2>/dev/null || echo "metadata.json not found, skipping copy."

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))

echo "✅ Simulation completed in $elapsed seconds."
echo "📁 Results saved to: $OUTPUT_DIR"