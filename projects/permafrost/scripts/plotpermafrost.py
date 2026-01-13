#!/usr/bin/env python3
import os
import glob
from igakit.io import PetIGA, VTK
import re
import numpy as np

# ========================
# User settings (edit me)
# ========================
# If True, sample the NURBS on a regular grid before writing VTK.
# We choose the number of samples per direction so that the *physical* spacing is ~square (dx≈dy).
RESAMPLE = True

# Target number of samples along the longest physical dimension (for 2D/3D).
# The other directions are chosen to keep dx≈dy (and dz≈dx when 3D).
NPTS_LONG = 400

# Output directory for VTK files
VTK_DIR = "vtkOut"


def extract_trailing_number(stem: str) -> str:
    """Extract trailing digits from a filename stem; return '' if none."""
    m = re.search(r"(\d+)$", stem)
    return m.group(1) if m else ""


def make_square_sampler(nrb, npts_long: int):
    """Return a sampler callable and (nu,nv,nw) chosen for ~square physical spacing.

    Igakit's VTK writer accepts `sampler=<callable>` where callable(U) returns a 1D array of
    parameter values for that parametric axis.

    We choose counts based on the physical extents estimated from the control points.
    """
    # Estimate physical extents from control points
    P = np.asarray(nrb.points)
    xyz_min = P[..., :3].reshape(-1, 3).min(axis=0)
    xyz_max = P[..., :3].reshape(-1, 3).max(axis=0)
    L = np.maximum(xyz_max - xyz_min, 1e-30)

    # Determine how many parametric directions we have (2D surface vs 3D volume)
    dim_param = len(nrb.knots)

    # Choose counts proportional to physical lengths to make spacing ~equal
    Lmax = float(L[:dim_param].max())
    if Lmax <= 0.0:
        Lmax = 1.0

    counts = []
    for d in range(dim_param):
        nd = int(round((L[d] / Lmax) * (npts_long - 1))) + 1
        counts.append(max(2, nd))

    # Build a sampler that returns uniform parametric samples per axis
    def sampler(U, _counts=counts):
        # U is the knot vector for this axis; pick count based on which knot vector it is
        # Igakit calls sampler separately for each axis, so we infer axis index by matching object id.
        # Fallback: use npts_long.
        try:
            axis = nrb.knots.index(U)
            n = _counts[axis]
        except Exception:
            n = npts_long
        return np.linspace(U[0], U[-1], n)

    # Pad counts for consistent reporting
    while len(counts) < 3:
        counts.append(1)

    return sampler, tuple(counts[:3])


def convert_solution_files(nrb_path, pattern, scalar_map, vtk_prefix, vtk_dir="vtkOut", resample=False, npts_long=400):
    """
    Convert solution .dat files to VTK format if not already present.

    Parameters:
        nrb_path (str): Path to the IGA geometry file.
        pattern (str): Glob pattern to match input solution files.
        scalar_map (dict): Mapping of scalar field names to their indices.
        vtk_prefix (str): Prefix for VTK output filenames.
        vtk_dir (str): Directory to store output VTK files.
        resample (bool): If True, sample onto a regular grid with ~square physical spacing.
        npts_long (int): Target samples along the longest physical dimension.
    """
    print(f"Processing files matching pattern: {pattern}")
    
    # Read IGA geometry
    try:
        nrb = PetIGA().read(nrb_path)
    except Exception as e:
        print(f"❌ Failed to read NURBS from {nrb_path}: {e}")
        return

    # Ensure output directory exists
    os.makedirs(vtk_dir, exist_ok=True)

    sampler = None
    if resample:
        sampler, (nu, nv, nw) = make_square_sampler(nrb, int(npts_long))
        if nw > 1:
            print(f"Resampling to ~square grid: {nu} x {nv} x {nw}")
        else:
            print(f"Resampling to ~square grid: {nu} x {nv}")

    for infile in glob.glob(pattern):
        name = os.path.splitext(os.path.basename(infile))[0]
        number = extract_trailing_number(name)
        # For files like "soil.dat" there may be no trailing number; that's OK.

        outfile = os.path.join(vtk_dir, f"{vtk_prefix}{number}.vtk")
        if not os.path.isfile(outfile):
            try:
                sol = PetIGA().read_vec(infile, nrb)
                if sampler is not None:
                    VTK().write(outfile, nrb, fields=sol, scalars=scalar_map, sampler=sampler)
                else:
                    VTK().write(outfile, nrb, fields=sol, scalars=scalar_map)
                print(f"✅ Wrote: {outfile}")
            except Exception as e:
                print(f"❌ Error processing {infile}: {e}")
        else:
            print(f"⏭️ Skipping existing file: {outfile}")

def main():
    convert_solution_files(
        nrb_path="igasol.dat",
        pattern="sol*.dat",
        scalar_map={'IcePhase': 0, 'Temperature': 1, 'VaporDensity': 2},
        vtk_prefix="solV_",
        vtk_dir=VTK_DIR,
        resample=RESAMPLE,
        npts_long=NPTS_LONG,
    )

    convert_solution_files(
        nrb_path="igasoil.dat",
        pattern="soil.dat",
        scalar_map={'Sediment': 0},
        vtk_prefix="soil",
        vtk_dir=VTK_DIR,
        resample=RESAMPLE,
        npts_long=NPTS_LONG,
    )

if __name__ == "__main__":
    main()