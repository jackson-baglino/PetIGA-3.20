#!/usr/bin/env python3
"""
convert_env_to_opts.py — one-time migration of .env files to PETSc .opts format.

Scans configs/porespy/ and configs/porespy2/ for .env files and writes
corresponding .opts files under inputs/geometry/.

Usage:
    python3 preprocess/convert_env_to_opts.py

Delete this script after migration is complete.
"""
import os
import glob

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS_GEOM = os.path.join(BASE_DIR, "inputs", "geometry")
os.makedirs(INPUTS_GEOM, exist_ok=True)

ENV_DIRS = [
    os.path.join(BASE_DIR, "configs", "porespy"),
    os.path.join(BASE_DIR, "configs", "porespy2"),
    os.path.join(BASE_DIR, "configs", "porespy2", "to_run"),
]

def parse_env(path):
    """Return a dict of key=value from a shell .env file."""
    params = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                params[k.strip()] = v.strip()
    return params

def infer_dim(params):
    """Infer dimension from Nz: Nz=1 → 2D, else 3D."""
    nz = int(params.get("Nz", "1"))
    return 2 if nz <= 1 else 3

converted = 0
for env_dir in ENV_DIRS:
    if not os.path.isdir(env_dir):
        continue
    for env_path in glob.glob(os.path.join(env_dir, "*.env")):
        stem = os.path.splitext(os.path.basename(env_path))[0]
        params = parse_env(env_path)
        dim = infer_dim(params)

        # Build the expected grains.dat path (heuristic; user may need to adjust)
        grain_dir = os.path.join(BASE_DIR, "inputs", stem)
        grains_file = os.path.join(grain_dir, "grains.dat")

        opts_path = os.path.join(INPUTS_GEOM, f"{stem}.opts")
        with open(opts_path, "w") as f:
            f.write(f"# Geometry config converted from {os.path.relpath(env_path, BASE_DIR)}\n")
            f.write(f"# Use with: -options_file inputs/solver.opts \\\n")
            f.write(f"#           -options_file inputs/geometry/{stem}.opts \\\n")
            f.write(f"#           -options_file inputs/experiment/<name>.opts\n\n")
            f.write(f"-dim {dim}\n")
            if "Nx" in params:
                f.write(f"-Nx {params['Nx']}\n")
            if "Ny" in params:
                f.write(f"-Ny {params['Ny']}\n")
            if "Nz" in params:
                f.write(f"-Nz {params['Nz']}\n")
            if "Lx" in params:
                f.write(f"-Lx {params['Lx']}\n")
            if "Ly" in params:
                f.write(f"-Ly {params['Ly']}\n")
            if "Lz" in params:
                f.write(f"-Lz {params['Lz']}\n")
            if "eps" in params:
                f.write(f"-eps {params['eps']}\n")
            f.write(f"-readFlag 1\n")
            f.write(f"-grains_file {grains_file}\n")

        print(f"  {os.path.relpath(env_path, BASE_DIR)} → inputs/geometry/{stem}.opts")
        converted += 1

print(f"\nConverted {converted} .env files to inputs/geometry/*.opts")
print("Review -grains_file paths above — adjust any that don't match your inputs/ layout.")
