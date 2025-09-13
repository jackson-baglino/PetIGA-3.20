#!/bin/bash
# batch_generate_grains.sh
#
# This script calls generateCircularGrains.py for multiple (Lx, Ly) pairs
# and multiple seeds. Results are saved inside outputs/ (dat, images, meta).
#
# Example usage:
#   bash batch_generate_grains.sh

# === PARAMETERS YOU CAN EDIT ===
# pairs=(
#   "1.0e-3 1.0e-3"
#   "1.5e-3 1.5e-3"
#   "2.0e-3 2.0e-3"
# )
pairs=("3.0e-3 3.0e-3")

seeds=(7 21 22)

porosities=(0.24)

shape_x=1024
shape_y=1024
mean_r_m=45e-6
sigma_ln=0.5
radius_clip_frac=1.0

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

      python preprocess/generateCircularGrains.py \
        --shape ${shape_x_cur} ${shape_y_cur} \
        --Lx ${Lx} --Ly ${Ly} \
        --porosity ${porosity} \
        --mean_r_m ${mean_r_m} --sigma_ln ${sigma_ln} \
        --radius_clip_frac ${radius_clip_frac} \
        --seed ${seed} \
        --require_connectivity

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
    "porosity": $porosity,
    "mean_r_m": $mean_r_m,x
    "sigma_ln": $sigma_ln,
    "radius_clip_frac": $radius_clip_frac,
    "seed": $seed
  },
  "cli": [
    "python", "preprocess/generateCircularGrains.py",
    "--Lx", "$Lx", "--Ly", "$Ly",
    "--porosity", "$porosity",
    "--mean_r_m", "$mean_r_m", "--sigma_ln", "$sigma_ln",
    "--radius_clip_frac", "$radius_clip_frac",
    "--seed", "$seed"
  ]
}
JSON

      echo "[meta] Wrote $meta_file"
    done
  done
done