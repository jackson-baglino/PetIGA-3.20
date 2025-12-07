import numpy as np
from petsc4py import PETSc
from scipy.ndimage import distance_transform_edt, gaussian_filter, zoom
import matplotlib.pyplot as plt


def load_petsc_vec(filename, comm=PETSc.COMM_SELF):
    """Load a PETSc Vec from a binary file using petsc4py."""
    viewer = PETSc.Viewer().createBinary(filename, 'r', comm=comm)
    vec = PETSc.Vec().create(comm=comm)
    vec.load(viewer)
    viewer.destroy()
    return vec


def save_petsc_vec(vec, filename, comm=PETSc.COMM_SELF):
    """Save a PETSc Vec to a binary file."""
    viewer = PETSc.Viewer().createBinary(filename, 'w', comm=comm)
    vec.view(viewer)
    viewer.destroy()


def vec_to_array_2d3(vec, Nx, Ny, dof=3):
    """Convert a PETSc Vec to a NumPy array of shape (Ny, Nx, dof).

    Assumes standard DMDA / IGA lexicographic ordering with x fastest,
    then y, and 'dof' components per node.

    total_size = Nx * Ny * dof  (this is checked).
    """
    Nglobal = vec.getSize()
    expected = Nx * Ny * dof
    if Nglobal != expected:
        raise ValueError(
            f"Global Vec size {Nglobal} != Nx*Ny*dof = {expected} "
            f"(Nx={Nx}, Ny={Ny}, dof={dof})"
        )

    # Get a copy of the global array
    arr_1d = vec.getArray().copy()
    # Reshape to (Ny, Nx, dof); x is fastest, then y
    arr = arr_1d.reshape((Ny, Nx, dof))
    return arr


def array_2d3_to_vec(arr, comm=PETSc.COMM_SELF):
    """Convert a NumPy array of shape (Ny, Nx, dof) to a PETSc Vec
    with the same global layout (x fastest, then y, then dof)."""
    Ny, Nx, dof = arr.shape
    Nglobal = Nx * Ny * dof

    vec = PETSc.Vec().create(comm=comm)
    vec.setSizes(Nglobal)
    vec.setFromOptions()

    flat = arr.reshape(-1)
    vec.setArray(flat.copy())
    return vec


def compute_signed_distance(binary_mask, dx=1.0, dy=1.0):
    """Compute a signed distance field from a binary mask.

    binary_mask: bool array where True = "inside" (e.g., ice)
    dx, dy:     grid spacing in x and y (for physical scaling)

    Returns:
        signed_dist (same shape as mask): negative inside, positive outside.
    """
    inside = binary_mask
    outside = ~binary_mask

    # distance_transform_edt computes distance in pixel units to the nearest zero.
    # Standard SDF construction:
    #   dist_outside = distance to closest True (inside)
    #   dist_inside  = distance to closest False (outside)
    # Then signed_dist = dist_outside - dist_inside
    # so that it is negative inside, positive outside.
    dist_outside = distance_transform_edt(outside, sampling=(dy, dx))
    dist_inside = distance_transform_edt(inside, sampling=(dy, dx))

    signed_dist = dist_outside - dist_inside
    return signed_dist


def reinitialize_phasefield_from_levelset(signed_dist, eps_new):
    """Given a signed distance field 'signed_dist' and a desired interface width
    eps_new, reconstruct a smooth phase-field using a tanh profile.

    We use the usual Cahn-Hilliard-style profile:
        phi(x) = 0.5 * (1 - tanh( s / (2*eps_new) ))
    where s < 0 is "inside" (phi ~ 1), s > 0 is "outside" (phi ~ 0).
    """
    return 0.5 * (1.0 - np.tanh(signed_dist / (2.0 * eps_new)))


def save_vtk_image(filename, ice, sediment, temp, Lx, Ly, upsample=1):
    Ny, Nx = ice.shape

    # Optional upsampling for smoother visualization in ParaView only (not for downsampling)
    if upsample > 1:
        from scipy.ndimage import zoom
        zoom_factor = float(upsample)
        ice_vtk = zoom(ice, zoom_factor, order=1)
        sediment_vtk = zoom(sediment, zoom_factor, order=1)
        temp_vtk = zoom(temp, zoom_factor, order=1)
        Ny_vtk, Nx_vtk = ice_vtk.shape
    else:
        ice_vtk = ice
        sediment_vtk = sediment
        temp_vtk = temp
        Ny_vtk, Nx_vtk = Ny, Nx

    # Use (Nx-1), (Ny-1) so that the physical size is Lx, Ly
    if Nx_vtk > 1:
        dx = Lx / (Nx_vtk - 1)
    else:
        dx = 1.0
    if Ny_vtk > 1:
        dy = Ly / (Ny_vtk - 1)
    else:
        dy = 1.0

    with open(filename, "w", encoding="ascii", newline="\n") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Reinitialized Phase Field Output\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {Nx_vtk} {Ny_vtk} 1\n")
        f.write("ORIGIN 0 0 0\n")
        f.write(f"SPACING {dx} {dy} 1\n")
        f.write(f"POINT_DATA {Nx_vtk * Ny_vtk}\n")

        # Ice field
        f.write("SCALARS IcePhase double\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(Ny_vtk):
            row = ice_vtk[j, :]
            f.write(" ".join(f"{val:.6e}" for val in row) + "\n")
        f.write("\n")

        # Sediment field
        f.write("SCALARS SedPhase double\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(Ny_vtk):
            row = sediment_vtk[j, :]
            f.write(" ".join(f"{val:.6e}" for val in row) + "\n")
        f.write("\n")

        # Temperature field
        f.write("SCALARS Temperature double\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(Ny_vtk):
            row = temp_vtk[j, :]
            f.write(" ".join(f"{val:.6e}" for val in row) + "\n")

    print(f"VTK file written to: {filename}")


def check_vtk_scalar_counts(path):
    """Quick sanity check that each SCALARS block has exactly POINT_DATA values."""
    with open(path, "r", encoding="ascii", errors="strict") as f:
        lines = f.readlines()

    # Find DIMENSIONS and POINT_DATA
    dims_line = next(l for l in lines if l.startswith("DIMENSIONS"))
    _, nx, ny, nz = dims_line.split()
    nx, ny, nz = int(nx), int(ny), int(nz)

    pd_line = next(l for l in lines if l.startswith("POINT_DATA"))
    _, npts = pd_line.split()
    npts = int(npts)

    print(f"DIMENSIONS: {nx} x {ny} x {nz}, POINT_DATA: {npts}")
    assert npts == nx * ny * nz, "POINT_DATA != DIMENSIONS product!"

    # Count how many float tokens each SCALARS block has
    i = lines.index(pd_line) + 1
    while i < len(lines):
        if lines[i].startswith("SCALARS"):
            name = lines[i].split()[1]
            i += 1  # skip LOOKUP_TABLE
            assert lines[i].startswith("LOOKUP_TABLE")
            i += 1
            values = []
            while i < len(lines) and not lines[i].startswith("SCALARS"):
                values.extend(lines[i].split())
                i += 1
            print(f"{name}: found {len(values)} values")
            assert len(values) == npts, f"{name} has wrong number of values"
        else:
            i += 1


def reinitialize_ice_phase(
    infile,
    outfile,
    Nx_in,
    Ny_in,
    Lx,
    Ly,
    eps_new,
    threshold=0.5,
    dof=3,
    ice_idx=0,
    sed_idx=1,
    temp_idx=2,
    comm=PETSc.COMM_SELF,
    plot=False,
    smooth_before_threshold=True,
    sigma_threshold=0.75,
    smooth_signed_dist=True,
    sigma_sdf=0.0,
    Nx_out=None,
    Ny_out=None,
):
    """High-level driver.

    - infile:    existing PETSc Vec file (e.g., [ph1, ph2, temp])
    - outfile:   output PETSc Vec file with reinitialized ice only
    - Nx_in, Ny_in:    input (fine) grid dimensions
    - Lx, Ly:    physical domain size (for dx, dy)
    - eps_new:   desired interface width for reconstructed phase-field
    - threshold: binarization threshold for current ice field
    - dof:       number of dofs (default 3: ice, sediment, temp)
    - ice_idx:   component index for ice in the Vec (default 0)
    - sed_idx:   component index for sediment in the Vec (default 1)
    - temp_idx:  component index for temperature in the Vec (default 2)
    - Nx_out, Ny_out: optional output (downsampled) grid size; if None, use input size

    Only the ice component is reinitialized; sediment and temperature
    are carried through unchanged.

    Note: VTK output is always generated at the output grid size.
    """
    print(f"Loading Vec from {infile} ...")
    vec_in = load_petsc_vec(infile, comm=comm)

    print("Converting Vec to array (Ny_in, Nx_in, dof)...")
    arr = vec_to_array_2d3(vec_in, Nx=Nx_in, Ny=Ny_in, dof=dof)

    # Extract existing fields using explicit component indices
    ice = arr[:, :, ice_idx]
    sediment = arr[:, :, sed_idx]
    temp = arr[:, :, temp_idx]

    # Grid spacing for signed distance computation
    dx = Lx / max(Nx_in - 1, 1)
    dy = Ly / max(Ny_in - 1, 1)

    # Optional smoothing before binarization to reduce pixel noise
    if smooth_before_threshold:
        print(f"Smoothing ice field before thresholding (sigma={sigma_threshold})...")
        ice_smooth = gaussian_filter(ice, sigma=sigma_threshold, mode="nearest")
    else:
        ice_smooth = ice

    print("Binarizing ice phase...")
    ice_binary = ice_smooth >= threshold

    print("Computing signed distance field...")
    signed_dist = compute_signed_distance(ice_binary, dx=dx, dy=dy)

    if smooth_signed_dist and sigma_sdf > 0.0:
        print(f"Smoothing signed distance field (sigma={sigma_sdf})...")
        signed_dist = gaussian_filter(signed_dist, sigma=sigma_sdf, mode="nearest")

    print(f"Reinitializing ice phase with eps_new = {eps_new:.3e} ...")
    ice_new = reinitialize_phasefield_from_levelset(signed_dist, eps_new=eps_new)

    # Clamp to [0,1] just for safety
    ice_new = np.clip(ice_new, 0.0, 1.0)

    # Decide output grid size (for PETSc Vec and VTK)
    if Nx_out is None:
        Nx_out = Nx_in
    if Ny_out is None:
        Ny_out = Ny_in

    # If requested, downsample fields from (Ny_in, Nx_in) to (Ny_out, Nx_out)
    if (Nx_out != Nx_in) or (Ny_out != Ny_in):
        print(f"Downsampling from ({Nx_in}, {Ny_in}) to ({Nx_out}, {Ny_out}) ...")

        # Compute zoom factors to map input grid -> output grid by coordinates
        # We use (Nx_out - 1)/(Nx_in - 1) so that the physical domain still
        # spans [0, Lx] x [0, Ly] after resampling.
        zoom_y = (Ny_out - 1) / max(Ny_in - 1, 1)
        zoom_x = (Nx_out - 1) / max(Nx_in - 1, 1)

        # Use linear interpolation (order=1) to avoid aliasing small features
        ice_use = zoom(ice_new, (zoom_y, zoom_x), order=1, mode="nearest")
        sed_use = zoom(sediment, (zoom_y, zoom_x), order=1, mode="nearest")
        temp_use = zoom(temp, (zoom_y, zoom_x), order=1, mode="nearest")

        # Make sure the shapes are exactly (Ny_out, Nx_out) after rounding
        Ny_out_eff, Nx_out_eff = ice_use.shape
        if (Ny_out_eff != Ny_out) or (Nx_out_eff != Nx_out):
            # Crop or pad by one cell if necessary due to rounding
            ice_use = ice_use[:Ny_out, :Nx_out]
            sed_use = sed_use[:Ny_out, :Nx_out]
            temp_use = temp_use[:Ny_out, :Nx_out]
    else:
        ice_use = ice_new
        sed_use = sediment
        temp_use = temp

    arr_new = np.zeros((Ny_out, Nx_out, dof), dtype=arr.dtype)

    # Safely handle minor shape mismatches after zoom (due to rounding)
    ny_i, nx_i = ice_use.shape
    ny_s, nx_s = sed_use.shape
    ny_t, nx_t = temp_use.shape

    ny_min = min(Ny_out, ny_i, ny_s, ny_t)
    nx_min = min(Nx_out, nx_i, nx_s, nx_t)

    arr_new[:ny_min, :nx_min, ice_idx] = ice_use[:ny_min, :nx_min]
    arr_new[:ny_min, :nx_min, sed_idx] = sed_use[:ny_min, :nx_min]
    arr_new[:ny_min, :nx_min, temp_idx] = temp_use[:ny_min, :nx_min]

    print(f"Ice (new)     range: [{ice_use.min():.3e}, {ice_use.max():.3e}]")
    print(f"Sediment      range: [{sed_use.min():.3e}, {sed_use.max():.3e}]")
    print(f"Temperature   range: [{temp_use.min():.3e}, {temp_use.max():.3e}]")
    print("Converting array back to Vec ...")
    vec_out = array_2d3_to_vec(arr_new, comm=comm)

    print(f"Saving reinitialized Vec to {outfile} ...")
    save_petsc_vec(vec_out, outfile, comm=comm)

    # Also save a VTK visualization file
    vtkfile = outfile.replace(".dat", ".vtk")
    print(f"Saving VTK file to {vtkfile} ...")
    save_vtk_image(
        vtkfile,
        arr_new[:, :, ice_idx],
        arr_new[:, :, sed_idx],
        arr_new[:, :, temp_idx],
        Lx,
        Ly,
        upsample=1,
    )

    # Clean up PETSc Vecs
    vec_in.destroy()
    vec_out.destroy()

    print("Done.")


if __name__ == "__main__":
    # Example usage – change these to match your case
    infile = "/Users/jacksonbaglino/PetIGA-3.20/projects/permafrost/inputs/tests/SthavisthaInputs/AssemblePhase.dat"
    outfile = "/Users/jacksonbaglino/PetIGA-3.20/projects/permafrost/inputs/tests/SthavisthaInputs/AssemblePhase_reinit.dat"

    # Input (fine) grid and physical domain – must match the original AssemblePhase.dat file
    Nx_in = 868 + 1
    Ny_in = 1736 + 1
    Lx = 4.04e-4
    Ly = 8.08e-4

    # Target (coarse) grid for the new Vec / VTK
    Nx_out = (434+1)
    Ny_out = (868+1)

    eps_new = 9.3e-7   # target epsilon
    threshold = 0.5

    reinitialize_ice_phase(
        infile=infile,
        outfile=outfile,
        Nx_in=Nx_in,
        Ny_in=Ny_in,
        Lx=Lx,
        Ly=Ly,
        eps_new=eps_new,
        threshold=threshold,
        dof=3,
        ice_idx=0,
        sed_idx=1,
        temp_idx=2,
        comm=PETSc.COMM_SELF,
        plot=True,
        smooth_before_threshold=True,
        sigma_threshold=1.0,
        smooth_signed_dist=True,
        sigma_sdf=2.0,
        Nx_out=Nx_out,
        Ny_out=Ny_out,
    )

    vtkfile = outfile.replace(".dat", ".vtk")
    check_vtk_scalar_counts(path=vtkfile)