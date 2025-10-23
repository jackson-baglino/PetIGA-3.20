#!/bin/bash
# batch_generate_grains.sh
#
# This script calls generateCircularGrains.py for multiple (Lx, Ly) pairs
# and multiple seeds. Results are saved inside outputs/ (dat, images, meta).
#
# Example usage:
#   bash batch_generate_grains.sh

# Resolve script directory so we can call the Python generator reliably
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYGEN="$SCRIPT_DIR/generateCircularGrains.py"

# === PARAMETERS YOU CAN EDIT ===
# pairs=(
#   "1.0e-3 1.0e-3"
#   "1.5e-3 1.5e-3"
#   "2.0e-3 2.0e-3"
# )
pairs=("2.0e-3 2.0e-3")

seeds=(7)

porosities=( 0.45)

shape_x=1024
shape_y=1024
mean_r_m=45e-6
sigma_ln=0.5
radius_clip_frac=1.0

# Two-phase controls (SOIL vs ICE)
soil_frac_solid=0.35   # fraction of TOTAL solids assigned to soil (0.0–1.0)
# Optional per-phase overrides (leave blank to inherit)
mean_r_m_ice=""        # e.g., 9.0e-5
sigma_ln_ice=""        # e.g., 0.35
radius_clip_frac_ice=""
mean_r_m_soil=""       # e.g., 6.0e-5
sigma_ln_soil=""       # e.g., 0.30
radius_clip_frac_soil=""

meta_dir="outputs/meta"
mkdir -p "$meta_dir"

# === MAIN LOOP ===
for pair in "${pairs[@]}"; do
  read -r Lx Ly <<< "$pair"

  # Set Nx, Ny based on (Lx, Ly) pair
  case "$Lx $Ly" in
    "1e-3 1e-3")
      shape_x_cur=571; shape_y_cur=571 ;;
    "1.5e-3 1.5e-3"|"1.5e-3 1.5e-4")
      shape_x_cur=857; shape_y_cur=857 ;;
    "2.0e-3 2.0e-3"|"2.e-3 2.0e-3"|"2e-3 2.0e-3"|"2.0e-3 2e-3"|"2.e-3 2.e-3")
      shape_x_cur=1024; shape_y_cur=1024 ;;
    *)
      # Fallback to global defaults if no explicit mapping
      shape_x_cur=$shape_x; shape_y_cur=$shape_y ;;
  esac
  echo "Resolved shape: Nx=${shape_x_cur}, Ny=${shape_y_cur} for Lx=${Lx}, Ly=${Ly}"
  # Note: Nx/Ny only affect preview image; DSM will compute grid internally.

  for seed in "${seeds[@]}"; do
    for porosity in "${porosities[@]}"; do
      echo "Running Lx=${Lx}, Ly=${Ly}, seed=${seed} with porosity=${porosity}"

      # Build argument list
      args=(
        --shape ${shape_x_cur} ${shape_y_cur}
        --Lx ${Lx} --Ly ${Ly}
        --porosity ${porosity}
        --mean_r_m ${mean_r_m} --sigma_ln ${sigma_ln}
        --radius_clip_frac ${radius_clip_frac}
        --soil_frac_solid ${soil_frac_solid}
        --seed ${seed}
        # --require_connectivity
      )
      # Optional per-phase overrides
      [[ -n "$mean_r_m_ice" ]] && args+=(--mean_r_m_ice "$mean_r_m_ice")
      [[ -n "$sigma_ln_ice" ]] && args+=(--sigma_ln_ice "$sigma_ln_ice")
      [[ -n "$radius_clip_frac_ice" ]] && args+=(--radius_clip_frac_ice "$radius_clip_frac_ice")
      [[ -n "$mean_r_m_soil" ]] && args+=(--mean_r_m_soil "$mean_r_m_soil")
      [[ -n "$sigma_ln_soil" ]] && args+=(--sigma_ln_soil "$sigma_ln_soil")
      [[ -n "$radius_clip_frac_soil" ]] && args+=(--radius_clip_frac_soil "$radius_clip_frac_soil")

      # Invoke generator (path-safe)
      python "$PYGEN" "${args[@]}"

      # --- Write informative metadata JSON for this run ---
      # Compute grid spacings and timestamp
      dx=$(awk -v Lx="$Lx" -v Nx="$shape_x_cur" 'BEGIN{printf "%.12g", Lx/Nx}')
      dy=$(awk -v Ly="$Ly" -v Ny="$shape_y_cur" 'BEGIN{printf "%.12g", Ly/Ny}')
      timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

      # Format porosity to 3 sig figs for slug
      slug_por=$(printf "%.3g" "$porosity" | sed 's/[-.]/_/g')
      meta_file="$meta_dir/meta_seed_${seed}_por_${slug_por}.json"

      # Emit JSON (numbers unquoted; booleans true/false)
      cat > "$meta_file" <<JSON
{
  "script": "generateCircularGrains.py",
  "timestamp_utc": "$timestamp",
  "inputs": {
    "Lx": $Lx,
    "Ly": $Ly,
    "dx": $dx,
    "dy": $dy,
    "shape_x": $shape_x_cur,
    "shape_y": $shape_y_cur,
    "porosity": $porosity,
    "soil_frac_solid": $soil_frac_solid,
    "mean_r_m": $mean_r_m,
    "sigma_ln": $sigma_ln,
    "radius_clip_frac": $radius_clip_frac,
    "mean_r_m_ice": "${mean_r_m_ice}",
    "sigma_ln_ice": "${sigma_ln_ice}",
    "radius_clip_frac_ice": "${radius_clip_frac_ice}",
    "mean_r_m_soil": "${mean_r_m_soil}",
    "sigma_ln_soil": "${sigma_ln_soil}",
    "radius_clip_frac_soil": "${radius_clip_frac_soil}",
    "seed": $seed
  },
  "command": "python $PYGEN --shape $shape_x_cur $shape_y_cur --Lx $Lx --Ly $Ly --porosity $porosity --mean_r_m $mean_r_m --sigma_ln $sigma_ln --radius_clip_frac $radius_clip_frac --soil_frac_solid $soil_frac_solid --seed $seed --require_connectivity"
}
JSON

      echo "[meta] Wrote $meta_file"
    done
  done
done