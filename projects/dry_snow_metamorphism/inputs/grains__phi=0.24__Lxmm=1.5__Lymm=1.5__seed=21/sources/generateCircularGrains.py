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
parser.add_argument("--eps", type=float, default=None,
                    help="Interface width epsilon (meters) to write into generated .env; default uses 0.5*dx")
parser.add_argument("--no_progress", action="store_true", help="Disable fancy progress bars (tqdm)")
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
MODE = "nonoverlap_poly"
CONTACT_GAP_M = 0.0
ALLOW_CROSS_BOUNDARY = True
ENFORCE_TANGENCY = True
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

def place_nonoverlap_polydisperse(nx, ny, target_solid_frac, rng, mu_ln, sigma_ln, gap_px, tries_per_circle, max_attempts, allow_cross_boundary, enforce_tangency, touch_tol_px, rmin_px, rmax_px):
    """Random sequential addition of non-overlapping disks with log-normal radii.
    Returns arrays (cx, cy, r_px) and a boolean mask 'solid'.
    """
    solid = np.zeros((nx, ny), dtype=bool)
    cx_list, cy_list, r_list = [], [], []
    solid_area_target = target_solid_frac * (nx * ny)
    area_now = 0.0
    attempts = 0
    print(f"Stopping conditions: reach target solids (area≥{solid_area_target:.0f} px) OR attempts≥{max_attempts}.")
    use_pb = (tqdm is not None) and (not args.no_progress)
    pbar_attempts = tqdm(total=max_attempts, desc="Attempts", unit="try", dynamic_ncols=True, leave=False) if use_pb else None
    pbar_solids = tqdm(total=int(solid_area_target), desc="Solid area", unit="px", dynamic_ncols=True, leave=False) if use_pb else None
    # Pre-allocate mask at the end for speed; use vector checks for overlap during placement
    while area_now < solid_area_target and attempts < max_attempts:
        if pbar_attempts:
            pbar_attempts.update(1)
        elif attempts % 10000 == 0 and attempts > 0:
            print(f"Attempt {attempts} / {max_attempts}")
        attempts += 1
        # sample radius (meters), convert to pixels
        r_m = rng.lognormal(mean=mu_ln, sigma=sigma_ln)
        r_px = int(max(1, round(r_m / dx)))
        if r_px < rmin_px or r_px > rmax_px:
            continue
        # quick skip if too large to fit
        if r_px*2 >= min(nx, ny):
            continue
        placed = False
        for _ in range(tries_per_circle):
            if allow_cross_boundary:
                # Center anywhere inside the image; disk may extend beyond bounds and will be clipped on rasterization
                x = rng.integers(0, nx)
                y = rng.integers(0, ny)
            else:
                # Constrain center so the entire disk fits inside the image
                x = rng.integers(r_px, nx - r_px)
                y = rng.integers(r_px, ny - r_px)
            if not cx_list:
                rr, cc = disk((x, y), r_px, shape=(nx, ny))
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
                    print(f"Placed grain #{len(r_list)} at ({x}, {y}) with radius {r_px} | progress: {progress*100:.1f}% of target solids | achieved porosity ≈ {achieved_porosity:.4f}")
                break
            # vectorized distance check
            dxs = np.array(cx_list) - x
            dys = np.array(cy_list) - y
            d2 = dxs*dxs + dys*dys
            min_allowed = (np.array(r_list) + r_px + gap_px)**2
            if np.all(d2 >= min_allowed):
                # Optionally expand radius to touch nearest neighbor (without overlap)
                if enforce_tangency and cx_list:
                    dists = np.sqrt(d2)
                    r_allowed = np.min(dists - (np.array(r_list) + gap_px))
                    r_target = int(max(1, round(r_allowed - touch_tol_px)))
                    if not allow_cross_boundary:
                        r_target = min(r_target, x, y, nx - 1 - x, ny - 1 - y)
                    # Clamp to requested radius bounds
                    r_target = max(rmin_px, min(r_target, rmax_px))
                    r_px = max(r_px, r_target)
                # Rasterize this grain and increment area using clipped area (respects boundary)
                rr, cc = disk((x, y), int(r_px), shape=(nx, ny))
                before = solid.sum()
                solid[rr, cc] = True
                added = solid.sum() - before
                placed = True
                cx_list.append(x); cy_list.append(y); r_list.append(int(r_px))
                area_now += added
                if pbar_solids:
                    pbar_solids.update(int(added))
                    pbar_solids.set_postfix_str(f"grains={len(r_list)} por≈{1.0 - (solid.sum()/(nx*ny)):.4f}")
                    pbar_solids.refresh()
                else:
                    progress = min(1.0, area_now / solid_area_target) if solid_area_target > 0 else 0.0
                    achieved_porosity = 1.0 - (solid.sum() / (nx * ny))
                    print(f"Placed grain #{len(r_list)} at ({x}, {y}) with radius {int(r_px)} | progress: {progress*100:.1f}% of target solids | achieved porosity ≈ {achieved_porosity:.4f}")
                break
        if not placed:
            continue
    if area_now >= solid_area_target:
        print(f"Target porosity reached with {len(r_list)} grains placed")
    else:
        achieved_porosity = 1.0 - (area_now / (nx * ny))
        print(f"Stopped after max attempts: placed {len(r_list)} grains, porosity = {achieved_porosity:.3f}")
    if pbar_attempts:
        pbar_attempts.close()
    if pbar_solids:
        pbar_solids.close()
    return np.array(cx_list), np.array(cy_list), np.array(r_list), solid

print(f"Generating {MODE} circular grains with shape={shape}, porosity={args.porosity}, seed={args.seed}")
# Only one mode supported: nonoverlap_poly (non-overlapping, polydisperse)
rng = np.random.default_rng(args.seed)
# For log-normal: mean = exp(mu + sigma^2/2) => mu = ln(mean) - sigma^2/2
mu_ln = np.log(args.mean_r_m) - (args.sigma_ln ** 2) / 2
gap_px = max(0, int(round(CONTACT_GAP_M / dx)))
target_solid = 1.0 - args.porosity
touch_tol_px = max(0, int(round(TOUCH_TOL_M / dx)))
rmin_m = max(0.0, args.mean_r_m * (1.0 - args.radius_clip_frac))
rmax_m = args.mean_r_m * (1.0 + args.radius_clip_frac)
rmin_px = max(1, int(round(rmin_m / dx)))
rmax_px = max(rmin_px, int(round(rmax_m / dx)))
cx_arr, cy_arr, r_arr, solid_mask = place_nonoverlap_polydisperse(
    nx, ny, target_solid, rng, mu_ln, args.sigma_ln, gap_px,
    TRIES_PER_CIRCLE, MAX_ATTEMPTS, ALLOW_CROSS_BOUNDARY,
    ENFORCE_TANGENCY, touch_tol_px, rmin_px, rmax_px
)
achieved_por = 1.0 - solid_mask.mean()
print(f"Target porosity={args.porosity:.3f}, achieved≈{achieved_por:.3f}; grains={len(r_arr)}")
sys.stdout.flush()

# --- Save .dat file ---
centers = np.column_stack((cx_arr * dx, cy_arr * dy))
radii = r_arr * dx
Lz = 2.2 * float(np.max(radii))  # meters (DSM convention)
poro_str = f"p{int(round(args.porosity*1000)):03d}"  # e.g., 0.38 -> p038
base_name = f"grainset_N{len(centers)}_seed{seed_str}_{poro_str}"

dat_path = os.path.join(inputs_run_dir, "grains.dat")
with open(dat_path, "w") as f:
    f.write(f"{args.Lx} {args.Ly} {Lz}\n")
    for (xC, yC), r in zip(centers, radii):
        zC = r / 2.0
        f.write(f"{xC} {yC} {zC} {r}\n")
print(f"Saved .dat to {dat_path}")

# --- Save .env file for downstream model ---
Nx = nx; Ny = ny
Nz = max(1, int(round(Lz / dx)))
eps_val = args.eps if args.eps is not None else 0.5 * dx
env_path = os.path.join(inputs_run_dir, "grains.env")
with open(env_path, "w") as fenv:
    fenv.write(f"Lx={args.Lx:.6e}\n")
    fenv.write(f"Ly={args.Ly:.6e}\n")
    fenv.write(f"Lz={Lz:.6e}\n")
print(f"Saved .env to {env_path}")

# Construct void mask for consistency with earlier logic
im = ~solid_mask

# im is boolean with True == void in many PoreSpy generators. For grains, invert if needed.
# According to PoreSpy docs, generators return True == pore (void). We want True == grains, so invert:
solid = ~im

img_path = os.path.join(inputs_run_dir, "preview.png")
base_name = os.path.splitext(os.path.basename(img_path))[0]

# -------------- Outputs --------------
# Only save .dat file already saved above for nonoverlap_poly mode
# No saving of .npy or other files

 # Save presentation-quality PNG
 # Compute numbers for concise title
n_grains = len(r_arr)
achieved_porosity = achieved_por

fig = plt.figure(figsize=(6.8, 6.8), dpi=150)
ax = fig.add_subplot(111)

# Use a clean binary colormap, square pixels, no axes
ax.imshow(solid, cmap="binary", interpolation="nearest", vmin=0, vmax=1)
ax.set_aspect("equal")
ax.set_axis_off()

# Short, mathy title: N, phi, R_ave
ax.set_title(rf"$N={n_grains}$   $\phi={achieved_porosity:.3f}$   $R_{{ave}}={args.mean_r_m*1e6:.0f}\,\mu m$",
             fontsize=12, pad=8)

# --- Add a simple scale bar (nice round length close to 20–30% of the short side) ---
short_side_m = min(args.Lx, args.Ly)
ideal_len_m = 0.25 * short_side_m
# choose a round micrometer length from common steps
choices_um = np.array([25, 50, 100, 200, 500, 1000, 2000])
length_um = choices_um[choices_um <= ideal_len_m * 1e6]
length_um = int(length_um[-1]) if length_um.size else int(max(25, round(ideal_len_m * 1e6)))
len_px = max(1, int(round((length_um * 1e-6) / dx)))
bar_h_px = max(2, int(round(0.015 * ny)))

# bottom-left margin
margin_x = int(round(0.06 * nx))
margin_y = int(round(0.06 * ny))
# High-contrast scale bar: white bar with black edge, outlined text for visibility
rect = patches.Rectangle((margin_x, ny - margin_y - bar_h_px), len_px, bar_h_px,
                         linewidth=1.5, edgecolor="black", facecolor="white", alpha=1.0)
ax.add_patch(rect)
txt = ax.text(margin_x + len_px / 2,
              ny - margin_y - bar_h_px - max(4, int(0.01 * ny)),
              f"{length_um} µm",
              fontsize=10, ha="center", va="bottom", color="black")
# Add a white outline so text stays visible over dark grains
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
                   "porosity_achieved": float(achieved_por if 'achieved_por' in locals() else (1.0 - solid.mean())),
                   "mean_r_m": float(args.mean_r_m), "sigma_ln": float(args.sigma_ln)},
    "artifacts": {"grains_dat": os.path.basename(dat_path) if 'dat_path' in locals() else None,
                   "env_file": os.path.basename(env_path) if 'env_path' in locals() else None,
                   "preview_png": os.path.basename(img_path)}
}
# Derived and informative fields
px_scale_m = float(dx)
rmin_m = float(max(0.0, args.mean_r_m * (1.0 - args.radius_clip_frac)))
rmax_m = float(args.mean_r_m * (1.0 + args.radius_clip_frac))
metadata.update({
    "pixel_scale_m": px_scale_m,
    # "shape_px": {"Nx": int(nx), "Ny": int(ny)},  # REMOVED
    "radii": {
        "distribution": "lognormal",
        "mean_r_m": float(args.mean_r_m),
        "sigma_ln": float(args.sigma_ln),
        "rmin_m": rmin_m,
        "rmax_m": rmax_m
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
        "num_grains": int(len(r_arr))
    },
    "artifacts_abs": {
        "inputs_dir": os.path.abspath(inputs_run_dir),
        "metadata_dir": os.path.abspath(inputs_run_dir),
        "grains_dat": os.path.abspath(dat_path) if 'dat_path' in locals() else None,
        "env_file": os.path.abspath(env_path) if 'env_path' in locals() else None,
        "preview_png": os.path.abspath(img_path)
    }
})
meta_json_path = os.path.join(inputs_run_dir, "metadata.json")
with open(meta_json_path, "w") as jf:
    json.dump(metadata, jf, indent=2)

# Checksums
chk_paths = [p for p in [meta_json_path, img_path, locals().get('dat_path'), locals().get('env_path')] if p]
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
        "timestamp,inputs_dir,metadata_dir,grains_dat,env_file,preview_png,"
        "mode,Lx,Ly,Lz,phi_target,phi_achieved,seed\n"
    )
    row = (
        f"{metadata['timestamp']},{os.path.abspath(inputs_run_dir)},{os.path.abspath(inputs_run_dir)},"
        f"{os.path.abspath(dat_path)},{os.path.abspath(env_path)},{os.path.abspath(img_path)},"
        f"{MODE},{args.Lx},{args.Ly},{Lz if 'Lz' in locals() else ''},"
        f"{args.porosity},{metadata['structure']['porosity_achieved']},{seed_str}\n"
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
else:
  plt.close(fig)