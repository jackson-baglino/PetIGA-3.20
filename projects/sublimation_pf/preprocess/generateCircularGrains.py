#!/usr/bin/env python3.9
"""
Script: generateCircularGrains.py
Generates non-overlapping, polydisperse circular grain packings and saves output in .dat and .png formats.

Example usage:
python generateCircularGrains.py \
  --shape 1024 1024 \
  --Lx 1e-3 --Ly 1e-3 \
  --porosity 0.38 \
  --mean_r_m 90e-6 --sigma_ln 0.25 \
  --radius_clip_frac 0.15 \
  --seed 7 \
  --show
"""
"""
Note: This script always generates non-overlapping, polydisperse circles using hardcoded packing parameters.
Outputs a boolean image where True == grain (solid) phase.
"""

import argparse
import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import patheffects as pe
from skimage.draw import disk
import sys
try:
    from tqdm import tqdm
except Exception:  # fallback if tqdm not installed
    tqdm = None

# -------------- CLI --------------
parser = argparse.ArgumentParser(description="Generate 2D circular grain packings with PoreSpy")
parser.add_argument("--shape", nargs=2, type=int, metavar=("NX", "NY"), default=(512, 512),
                    help="Image size in pixels (Nx Ny). Default: 512 512")
parser.add_argument("--porosity", type=float, default=0.4, help="Target porosity (void fraction)")
parser.add_argument("--r", type=int, default=8, help="Circle radius for mono-disperse modes")
parser.add_argument("--rmin", type=int, default=4, help="Minimum radius for polydisperse mode")
parser.add_argument("--nbins", type=int, default=6, help="Number of radius bins for polydisperse mode")
parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
parser.add_argument("--save", type=str, default=None, help="Optional path to save the boolean array (.npy)")
parser.add_argument("--png", type=str, default=None, help="Optional path to save a PNG preview")
parser.add_argument("--dat", type=str, default=None, help="Optional path to save centers/radii as .dat (meters): header Lx Ly Lz; then xC yC zC radius per row")
parser.add_argument("--show", action="store_true", help="Show plot interactively")
parser.add_argument("--Lx", type=float, default=1.0e-3, help="Domain length in x (meters)")
parser.add_argument("--Ly", type=float, default=1.0e-3, help="Domain length in y (meters)")
parser.add_argument("--mean_r_m", type=float, default=90e-6, help="Mean grain radius in meters for log-normal distribution")
parser.add_argument("--sigma_ln", type=float, default=0.35, help="Log-normal sigma for grain radii (in log space)")
parser.add_argument("--radius_clip_frac", type=float, default=0.5, help="Limit radii to within (1±frac)*mean_r_m for polydisperse modes (default 0.5 => ±50%)")
parser.add_argument("--soil_frac_solid", type=float, default=0.50,
                    help="Fraction of the *total solid* assigned to soil (rest goes to ice). E.g., 0.3 => 30% soil solids, 70% ice solids. Default: 0.50")
# Phase-specific radius distribution overrides (soil defaults to ice values if not set)
parser.add_argument("--mean_r_m_ice", type=float, default=None,
                    help="Override mean radius for ICE (m); if None, uses --mean_r_m")
parser.add_argument("--sigma_ln_ice", type=float, default=None,
                    help="Override log-normal sigma for ICE; if None, uses --sigma_ln")
parser.add_argument("--radius_clip_frac_ice", type=float, default=None,
                    help="Override clip frac for ICE; if None, uses --radius_clip_frac")
parser.add_argument("--mean_r_m_soil", type=float, default=None,
                    help="Mean radius for SOIL (m); if None, uses ICE's value")
parser.add_argument("--sigma_ln_soil", type=float, default=None,
                    help="Log-normal sigma for SOIL; if None, uses ICE's value")
parser.add_argument("--radius_clip_frac_soil", type=float, default=None,
                    help="Clip frac for SOIL; if None, uses ICE's value")
parser.add_argument("--eps", type=float, default=None,
                    help="Interface width epsilon (meters) to write into generated .env; default uses 0.5*dx")
parser.add_argument("--no_progress", action="store_true", help="Disable fancy progress bars (tqdm)")
parser.add_argument("--require_connectivity", action="store_true",
                    help="Only accept placements that make contact with existing grains; ensures a single connected solid cluster.")
#
# Default roots: resolve relative to the project root (one level above script if in ./preprocess/)
default_project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
parser.add_argument("--inputs-root", type=str, default=os.path.join(default_project_root, "inputs"),
                    help="Root folder to store generated INPUT files (.dat, .env) in run-specific subfolders. Default: ../inputs relative to preprocess/")
parser.add_argument("--metadata-root", type=str, default=os.path.join(default_project_root, "metadata"),
                    help="Root folder to store METADATA for each run (preview.png, metadata.json, checksums.txt, sources/). Default: ../metadata relative to preprocess/")
parser.add_argument("--out-root", type=str, default=None,
                    help="[DEPRECATED] If provided, inputs go to <out-root>/inputs and metadata to <out-root>/metadata.")
args = parser.parse_args()

# Resolve output roots (support deprecated --out-root). All artifacts go under inputs_root.
if args.out_root:
    inputs_root = os.path.join(args.out_root, "inputs")
else:
    inputs_root = args.inputs_root
metadata_root = inputs_root  # unify: save everything under inputs
os.makedirs(inputs_root, exist_ok=True)

# --- Hardcoded packing parameters (previously CLI flags)
MODE = "interleaved"
CONTACT_GAP_M = 1.0e-6
ALLOW_CROSS_BOUNDARY = True
ENFORCE_TANGENCY = False
TOUCH_TOL_M = 0.0
TRIES_PER_CIRCLE = 20000
MAX_ATTEMPTS = 2000000

print("Starting grain generation with parameters:")
print(args)
print({
    'MODE': MODE,
    'CONTACT_GAP_M': CONTACT_GAP_M,
    'ALLOW_CROSS_BOUNDARY': ALLOW_CROSS_BOUNDARY,
    'ENFORCE_TANGENCY': ENFORCE_TANGENCY,
    'TOUCH_TOL_M': TOUCH_TOL_M,
    'TRIES_PER_CIRCLE': TRIES_PER_CIRCLE,
    'MAX_ATTEMPTS': MAX_ATTEMPTS,
})

# --- Structured run directory (timestamp + tags) ---
def _fmt_sig(x, sig=6):
    try:
        return ("%g" % float(f"{float(x):.{sig}g}"))
    except Exception:
        return str(x)

ts = datetime.now().strftime("%Y-%m-%d:%H-%M-%S")
seed_str = str(args.seed) if args.seed is not None else "NA"
# Build concise run name: grains__phi=...__Lxmm=...__Lymm=...__seed=...
phi_str = f"{args.porosity:.3g}"
Lx_mm = args.Lx * 1e3
Ly_mm = args.Ly * 1e3
run_name = (
    f"grains__phi={phi_str}__Lxmm={_fmt_sig(Lx_mm)}__Lymm={_fmt_sig(Ly_mm)}__seed={seed_str}"
)
inputs_run_dir = os.path.join(inputs_root, run_name)
os.makedirs(inputs_run_dir, exist_ok=True)
# Keep a copy of the script for provenance (now stored under inputs_run_dir)
sources_dir = os.path.join(inputs_run_dir, "sources"); os.makedirs(sources_dir, exist_ok=True)
try:
    _self_path = os.path.abspath(__file__)
    if os.path.isfile(_self_path):
        import shutil
        shutil.copy2(_self_path, os.path.join(sources_dir, os.path.basename(_self_path)))
except Exception:
    pass
# Maintain convenient "latest" symlinks for inputs
def _update_latest_link(link_path, target):
    try:
        if os.path.islink(link_path) or os.path.exists(link_path):
            try:
                os.remove(link_path)
            except Exception:
                pass
        os.symlink(target, link_path)
    except Exception:
        pass
_update_latest_link(os.path.join(inputs_root, "latest"), inputs_run_dir)

nx, ny = args.shape
shape = (nx, ny)  # 2D -> circles

# --- Physical mapping (meters ↔ pixels)
dx = args.Lx / nx
dy = args.Ly / ny
if abs(dx - dy) / dx > 1e-6:
    raise ValueError("Non-square pixels detected (Lx/nx != Ly/ny). Use equal aspect for now.")
m2px = lambda r_m: max(1, int(round(r_m / dx)))

def place_nonoverlap_polydisperse(nx, ny, target_solid_frac, rng, mu_ln, sigma_ln,
                                  gap_px, tries_per_circle, max_attempts,
                                  allow_cross_boundary, enforce_tangency, touch_tol_px,
                                  rmin_px, rmax_px,
                                  init_cx=None, init_cy=None, init_r=None, init_solid=None):
    """Random sequential addition of non-overlapping disks with log-normal radii.
    Optionally seed with an existing packing (init_* and init_solid) so new grains
    avoid overlaps across phases. Returns (cx, cy, r_px, solid_mask).
    """
    if init_solid is not None:
        solid = init_solid.copy()
    else:
        solid = np.zeros((nx, ny), dtype=bool)
    cx_list = list([] if init_cx is None else np.asarray(init_cx, dtype=int))
    cy_list = list([] if init_cy is None else np.asarray(init_cy, dtype=int))
    r_list  = list([] if init_r  is None else np.asarray(init_r,  dtype=int))

    area_now = float(solid.sum())
    solid_area_target = target_solid_frac * (nx * ny)

    if getattr(args, 'require_connectivity', False) and len(cx_list) == 0:
        print("[CONNECTIVITY] Enforcing contact: every new grain must touch the existing grains (single connected cluster).")

    attempts = 0
    print(f"Stopping conditions: reach target solids (area≥{solid_area_target:.0f} px) OR attempts≥{max_attempts}.")
    use_pb = (tqdm is not None) and (not args.no_progress)
    pbar_attempts = tqdm(total=max_attempts, desc="Attempts", unit="try", dynamic_ncols=True, leave=False) if use_pb else None
    pbar_solids   = tqdm(total=max(0, int(solid_area_target - area_now)), desc="Solid area", unit="px", dynamic_ncols=True, leave=False) if use_pb else None

    while area_now < solid_area_target and attempts < max_attempts:
        if pbar_attempts:
            pbar_attempts.update(1)
        elif attempts % 10000 == 0 and attempts > 0:
            print(f"Attempt {attempts} / {max_attempts}")
        attempts += 1

        r_m = rng.lognormal(mean=mu_ln, sigma=sigma_ln)
        r_px = int(max(1, round(r_m / dx)))
        if r_px < rmin_px or r_px > rmax_px:
            continue
        if r_px*2 >= min(nx, ny):
            continue

        placed = False
        for _ in range(tries_per_circle):
            if allow_cross_boundary:
                x = rng.integers(0, nx)
                y = rng.integers(0, ny)
            else:
                x = rng.integers(r_px, nx - r_px)
                y = rng.integers(r_px, ny - r_px)

            if len(cx_list) == 0:
                rr, cc = disk((x, y), r_px, shape=(nx, ny))
                if solid[rr, cc].any():
                    continue
                before = solid.sum()
                solid[rr, cc] = True
                added = solid.sum() - before
                placed = True
                cx_list.append(x); cy_list.append(y); r_list.append(r_px)
                area_now += added
                if pbar_solids:
                    pbar_solids.update(int(added))
                    pbar_solids.set_postfix_str(f"grains={len(r_list)} por≈{1.0 - (solid.sum()/(nx*ny)):.4f}")
                    pbar_solids.refresh()
                else:
                    progress = min(1.0, area_now / solid_area_target) if solid_area_target > 0 else 0.0
                    achieved_porosity = 1.0 - (solid.sum() / (nx * ny))
                    print(f"Placed grain #{len(r_list)} at ({x}, {y}) with radius {r_px} | progress: {progress*100:.1f}% | achieved porosity ≈ {achieved_porosity:.4f}")
                break

            dxs = np.array(cx_list) - x
            dys = np.array(cy_list) - y
            d2 = dxs*dxs + dys*dys
            min_allowed = (np.array(r_list) + r_px + gap_px) ** 2
            if np.all(d2 >= min_allowed):
                r_px_final = int(r_px)
                dists = np.sqrt(d2)
                idx_min = int(np.argmin(dists))
                r_needed_raw = dists[idx_min] - (r_list[idx_min] + gap_px)
                r_needed = int(max(1, np.floor(r_needed_raw - touch_tol_px)))
                if not allow_cross_boundary:
                    r_needed = min(r_needed, x, y, nx - 1 - x, ny - 1 - y)
                r_needed = max(rmin_px, min(r_needed, rmax_px))

                if getattr(args, 'require_connectivity', False):
                    if r_needed + touch_tol_px < r_needed_raw:
                        continue
                    r_px_final = r_needed
                elif enforce_tangency:
                    r_px_candidate = max(r_px_final, r_needed)
                    min_allowed_grown = (np.array(r_list) + r_px_candidate + gap_px) ** 2
                    if np.all(d2 >= min_allowed_grown):
                        if not allow_cross_boundary:
                            r_px_candidate = min(r_px_candidate, x, y, nx - 1 - x, ny - 1 - y)
                        r_px_final = int(r_px_candidate)

                rr, cc = disk((x, y), int(r_px_final), shape=(nx, ny))
                if solid[rr, cc].any():
                    continue

                before = solid.sum()
                solid[rr, cc] = True
                added = solid.sum() - before
                placed = True
                cx_list.append(x); cy_list.append(y); r_list.append(int(r_px_final))
                area_now += added
                if pbar_solids:
                    pbar_solids.update(int(added))
                    pbar_solids.set_postfix_str(f"grains={len(r_list)} por≈{1.0 - (solid.sum()/(nx*ny)):.4f}")
                    pbar_solids.refresh()
                else:
                    progress = min(1.0, area_now / solid_area_target) if solid_area_target > 0 else 0.0
                    achieved_porosity = 1.0 - (solid.sum() / (nx * ny))
                    print(f"Placed grain #{len(r_list)} at ({x}, {y}) with radius {int(r_px_final)} | progress: {progress*100:.1f}% | achieved porosity ≈ {achieved_porosity:.4f}")
                break
        if not placed:
            continue

    if area_now >= solid_area_target:
        print(f"Target reached with {len(r_list)} total grains (including seeded)")
    else:
        achieved_porosity = 1.0 - (area_now / (nx * ny))
        print(f"Stopped after max attempts: total grains = {len(r_list)}, porosity = {achieved_porosity:.3f}")

    if pbar_attempts: pbar_attempts.close()
    if pbar_solids:   pbar_solids.close()

    if getattr(args, 'require_connectivity', False):
        try:
            import scipy.ndimage as ndi
            labels, ncomp = ndi.label(solid)
            if ncomp > 1:
                print(f"[CONNECTIVITY][WARN] Final solid has {ncomp} connected components (expected 1). Consider increasing tries_per_circle or max_attempts.")
        except Exception:
            pass

    return np.array(cx_list, dtype=int), np.array(cy_list, dtype=int), np.array(r_list, dtype=int), solid


print(f"Generating {MODE} circular grains with shape={shape}, porosity={args.porosity}, seed={args.seed}")
rng = np.random.default_rng(args.seed)

# ---- Per-phase parameter setup ----
mean_r_ice = args.mean_r_m if args.mean_r_m_ice is None else args.mean_r_m_ice
sigma_ice   = args.sigma_ln if args.sigma_ln_ice is None else args.sigma_ln_ice
clip_ice    = args.radius_clip_frac if args.radius_clip_frac_ice is None else args.radius_clip_frac_ice

mean_r_soil = mean_r_ice if args.mean_r_m_soil is None else args.mean_r_m_soil
sigma_soil   = sigma_ice if args.sigma_ln_soil is None else args.sigma_ln_soil
clip_soil    = clip_ice  if args.radius_clip_frac_soil is None else args.radius_clip_frac_soil

mu_ln_ice  = np.log(mean_r_ice)  - (sigma_ice  ** 2) / 2
mu_ln_soil = np.log(mean_r_soil) - (sigma_soil ** 2) / 2

gap_px = max(0, int(round(CONTACT_GAP_M / dx)))
touch_tol_px = max(0, int(round(TOUCH_TOL_M / dx)))

solid_total = 1.0 - args.porosity
soil_frac   = np.clip(args.soil_frac_solid, 0.0, 1.0)
solid_soil_target = soil_frac * solid_total
solid_ice_target  = (1.0 - soil_frac) * solid_total

rmin_ice_m = max(0.0, mean_r_ice * (1.0 - clip_ice))
rmax_ice_m =            mean_r_ice * (1.0 + clip_ice)
rmin_soil_m = max(0.0, mean_r_soil * (1.0 - clip_soil))
rmax_soil_m =            mean_r_soil * (1.0 + clip_soil)

rmin_ice_px  = max(1, int(round(rmin_ice_m  / dx)))
rmax_ice_px  = max(rmin_ice_px,  int(round(rmax_ice_m  / dx)))
rmin_soil_px = max(1, int(round(rmin_soil_m / dx)))
rmax_soil_px = max(rmin_soil_px, int(round(rmax_soil_m / dx)))

# ---- Two-phase placement ----
cx_ice, cy_ice, r_ice, solid_mask = place_nonoverlap_polydisperse(
    nx, ny, solid_ice_target, rng, mu_ln_ice, sigma_ice, gap_px,
    TRIES_PER_CIRCLE, MAX_ATTEMPTS, ALLOW_CROSS_BOUNDARY,
    ENFORCE_TANGENCY, touch_tol_px, rmin_ice_px, rmax_ice_px,
)

cx_soil, cy_soil, r_soil, solid_mask = place_nonoverlap_polydisperse(
    nx, ny, solid_ice_target + solid_soil_target, rng, mu_ln_soil, sigma_soil, gap_px,
    TRIES_PER_CIRCLE, MAX_ATTEMPTS, ALLOW_CROSS_BOUNDARY,
    ENFORCE_TANGENCY, touch_tol_px, rmin_soil_px, rmax_soil_px,
    init_cx=cx_ice, init_cy=cy_ice, init_r=r_ice, init_solid=solid_mask
)

combined_mask = solid_mask.copy()
solid_ice_mask = np.zeros_like(combined_mask)
solid_soil_mask = np.zeros_like(combined_mask)
for x, y, rpx in zip(cx_ice, cy_ice, r_ice):
    rr, cc = disk((int(x), int(y)), int(rpx), shape=(nx, ny))
    solid_ice_mask[rr, cc] = True
for x, y, rpx in zip(cx_soil[len(cx_ice):], cy_soil[len(cy_ice):], r_soil[len(r_ice):]):
    rr, cc = disk((int(x), int(y)), int(rpx), shape=(nx, ny))
    solid_soil_mask[rr, cc] = True

achieved_por = 1.0 - combined_mask.mean()
print(f"Target porosity={args.porosity:.3f}, achieved≈{achieved_por:.3f}; ice grains={len(cx_ice)}, soil grains={len(cx_soil)-len(cx_ice)}")
sys.stdout.flush()


# --- Save ICE & SOIL .dat files ---
centers_ice_m  = np.column_stack((np.asarray(cx_ice) * dx, np.asarray(cy_ice) * dy))
radii_ice_m    = np.asarray(r_ice) * dx
centers_soil_all = np.column_stack((np.asarray(cx_soil) * dx, np.asarray(cy_soil) * dy))
radii_soil_all   = np.asarray(r_soil) * dx
n_ice = len(cx_ice)
centers_soil_m = centers_soil_all[n_ice:]
radii_soil_m   = radii_soil_all[n_ice:]

Lz = 2.2 * float(max(np.max(radii_ice_m) if radii_ice_m.size else 0.0,
                     np.max(radii_soil_m) if radii_soil_m.size else 0.0))

ice_dat_path  = os.path.join(inputs_run_dir, "grains_ice.dat")
soil_dat_path = os.path.join(inputs_run_dir, "grains_soil.dat")

with open(ice_dat_path, "w") as f:
    f.write(f"{args.Lx} {args.Ly} {Lz}\n")
    for (xC, yC), r in zip(centers_ice_m, radii_ice_m):
        zC = r / 2.0
        f.write(f"{xC} {yC} {zC} {r}\n")
with open(soil_dat_path, "w") as f:
    f.write(f"{args.Lx} {args.Ly} {Lz}\n")
    for (xC, yC), r in zip(centers_soil_m, radii_soil_m):
        zC = r / 2.0
        f.write(f"{xC} {yC} {zC} {r}\n")
print(f"Saved ICE .dat to {ice_dat_path}")
print(f"Saved SOIL .dat to {soil_dat_path}")

# --- Save .env file for downstream model ---
Nx = nx; Ny = ny
Nz = max(1, int(round(Lz / dx)))
eps_val = args.eps if args.eps is not None else 0.5 * dx
env_path = os.path.join(inputs_run_dir, "grains.env")
with open(env_path, "w") as fenv:
    fenv.write(f"Lx={args.Lx:.6e}\n")
    fenv.write(f"Ly={args.Ly:.6e}\n")
    fenv.write(f"Lz={Lz:.6e}\n")
print(f"Saved grains.env to {env_path}")

img_path = os.path.join(inputs_run_dir, "preview.png")
base_name = os.path.splitext(os.path.basename(img_path))[0]

# -------------- Outputs --------------
# Two-phase visualization
viz = np.zeros((nx, ny), dtype=int)
viz[solid_ice_mask] = 1
viz[solid_soil_mask] = 2

img_path = os.path.join(inputs_run_dir, "preview.png")

fig = plt.figure(figsize=(6.8, 6.8), dpi=150)
ax = fig.add_subplot(111)
from matplotlib.colors import ListedColormap
ax.imshow(viz, cmap=ListedColormap(["white", "lightgray", "black"]), interpolation="nearest", vmin=0, vmax=2)
ax.set_aspect("equal")
ax.set_axis_off()

n_grains_ice  = int(n_ice)
n_grains_soil = int(len(cx_soil) - n_ice)
ax.set_title(rf"$N_{{ice}}={n_grains_ice}$   $N_{{soil}}={n_grains_soil}$   $\phi={achieved_por:.3f}$",
             fontsize=12, pad=8)

short_side_m = min(args.Lx, args.Ly)
ideal_len_m = 0.25 * short_side_m
choices_um = np.array([25, 50, 100, 200, 500, 1000, 2000])
length_um = choices_um[choices_um <= ideal_len_m * 1e6]
length_um = int(length_um[-1]) if length_um.size else int(max(25, round(ideal_len_m * 1e6)))
len_px = max(1, int(round((length_um * 1e-6) / dx)))
bar_h_px = max(2, int(round(0.015 * ny)))
margin_x = int(round(0.06 * nx)); margin_y = int(round(0.06 * ny))
rect = patches.Rectangle((margin_x, ny - margin_y - bar_h_px), len_px, bar_h_px,
                         linewidth=1.5, edgecolor="black", facecolor="white", alpha=1.0)
ax.add_patch(rect)
txt = ax.text(margin_x + len_px / 2,
              ny - margin_y - bar_h_px - max(4, int(0.01 * ny)),
              f"{length_um} µm",
              fontsize=10, ha="center", va="bottom", color="black")
txt.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])

fig.tight_layout(pad=0.2)
plt.savefig(img_path, dpi=300)
print(f"Saved image to {img_path}")
plt.close(fig)

# --- Save metadata & checksums ---

import json, hashlib
metadata = {
    "timestamp": datetime.now().astimezone().isoformat(),
    "generator": {
        "name": "porespy",
        "script": os.path.basename(__file__),
        "mode": MODE,
        "seed": None if args.seed is None else int(args.seed)
    },
    "domain": {"Lx": float(args.Lx), "Ly": float(args.Ly), "Lz": float(Lz) if 'Lz' in locals() else None},
    "structure": {"porosity_target": float(args.porosity),
                   "porosity_achieved": float(achieved_por if 'achieved_por' in locals() else (1.0 - combined_mask.mean())),
                   "mean_r_m_ice": float(mean_r_ice), "sigma_ln_ice": float(sigma_ice),
                   "mean_r_m_soil": float(mean_r_soil), "sigma_ln_soil": float(sigma_soil),
                   "soil_frac_solid": float(soil_frac)},
    "artifacts": {
        "grains_ice_dat": os.path.basename(ice_dat_path),
        "grains_soil_dat": os.path.basename(soil_dat_path),
        "env_file": os.path.basename(env_path) if 'env_path' in locals() else None,
        "preview_png": os.path.basename(img_path)
    }
}
# Derived and informative fields
px_scale_m = float(dx)
rmin_ice_m = float(rmin_ice_m)
rmax_ice_m = float(rmax_ice_m)
rmin_soil_m = float(rmin_soil_m)
rmax_soil_m = float(rmax_soil_m)
metadata.update({
    "pixel_scale_m": px_scale_m,
    "radii_ice": {
        "distribution": "lognormal",
        "mean_r_m": float(mean_r_ice),
        "sigma_ln": float(sigma_ice),
        "rmin_m": rmin_ice_m,
        "rmax_m": rmax_ice_m
    },
    "radii_soil": {
        "distribution": "lognormal",
        "mean_r_m": float(mean_r_soil),
        "sigma_ln": float(sigma_soil),
        "rmin_m": rmin_soil_m,
        "rmax_m": rmax_soil_m
    },
    "packing": {
        "mode": MODE,
        "contact_gap_m": float(CONTACT_GAP_M),
        "allow_cross_boundary": bool(ALLOW_CROSS_BOUNDARY),
        "enforce_tangency": bool(ENFORCE_TANGENCY),
        "touch_tol_m": float(TOUCH_TOL_M),
        "tries_per_circle": int(TRIES_PER_CIRCLE),
        "max_attempts": int(MAX_ATTEMPTS)
    },
    "counts": {
        "num_grains_ice": int(len(cx_ice)),
        "num_grains_soil": int(len(cx_soil) - len(cx_ice))
    },
    "artifacts_abs": {
        "inputs_dir": os.path.abspath(inputs_run_dir),
        "metadata_dir": os.path.abspath(inputs_run_dir),
        "grains_ice_dat": os.path.abspath(ice_dat_path),
        "grains_soil_dat": os.path.abspath(soil_dat_path),
        "env_file": os.path.abspath(env_path) if 'env_path' in locals() else None,
        "preview_png": os.path.abspath(img_path)
    }
})
meta_json_path = os.path.join(inputs_run_dir, "metadata.json")
with open(meta_json_path, "w") as jf:
    json.dump(metadata, jf, indent=2)

# Checksums
chk_paths = [p for p in [meta_json_path, img_path, ice_dat_path, soil_dat_path, locals().get('env_path')] if p]
chk_text = []
for p in chk_paths:
    try:
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            h.update(fh.read())
        chk_text.append(f"{h.hexdigest()}  {os.path.basename(p)}")
    except Exception:
        pass
with open(os.path.join(inputs_run_dir, "checksums.txt"), "w") as cf:
    cf.write("\n".join(chk_text) + "\n")

print(f"Saved metadata to {meta_json_path}")

# Manifest (optional): manifest.csv under inputs_root
try:
    man_path = os.path.join(inputs_root, "manifest.csv")
    hdr = (
        "timestamp,inputs_dir,metadata_dir,grains_ice_dat,grains_soil_dat,env_file,preview_png,"
        "mode,Lx,Ly,Lz,phi_target,phi_achieved,soil_frac_solid,seed\n"
    )
    row = (
        f"{metadata['timestamp']},{os.path.abspath(inputs_run_dir)},{os.path.abspath(inputs_run_dir)},"
        f"{os.path.abspath(ice_dat_path)},{os.path.abspath(soil_dat_path)},{os.path.abspath(env_path)},{os.path.abspath(img_path)},"
        f"{MODE},{args.Lx},{args.Ly},{Lz if 'Lz' in locals() else ''},"
        f"{args.porosity},{metadata['structure']['porosity_achieved']},{soil_frac},{seed_str}\n"
    )
    need_hdr = not os.path.exists(man_path)
    with open(man_path, "a") as mf:
        if need_hdr:
            mf.write(hdr)
        mf.write(row)
except Exception:
    pass

if args.show:
  plt.show()