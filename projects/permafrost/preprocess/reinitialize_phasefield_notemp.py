import numpy as np
from petsc4py import PETSc
from scipy.ndimage import distance_transform_edt, gaussian_filter, zoom
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Basic PETSc Vec I/O helpers
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# Vec <-> NumPy array conversion (2D grid, Ndof per point)
# -----------------------------------------------------------------------------

def vec_to_array_2d(vec, Nx, Ny, dof):
    """Convert a PETSc Vec to a NumPy array of shape (Ny, Nx, dof).

    Assumes standard lexicographic ordering: x fastest, then y, then dof.

    total_size = Nx * Ny * dof  (this is checked).
    """
    Nglobal = vec.getSize()
    expected = Nx * Ny * dof
    if Nglobal != expected:
        raise ValueError(
            f"Global Vec size {Nglobal} != Nx*Ny*dof = {expected} "
            f"(Nx={Nx}, Ny={Ny}, dof={dof})"
        )

    arr_1d = vec.getArray().copy()
    arr = arr_1d.reshape((Ny, Nx, dof))
    return arr


def array_2d_to_vec(arr, comm=PETSc.COMM_SELF):
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


# -----------------------------------------------------------------------------
# Geometry helpers: signed distance + phase-field reconstruction
# -----------------------------------------------------------------------------

def compute_signed_distance(binary_mask, dx=1.0, dy=1.0):
    """Compute a signed distance field from a binary mask.

    binary_mask: bool array where True = "inside" (e.g., ice)
    dx, dy:     grid spacing in x and y (for physical scaling)

    Returns:
        signed_dist (same shape as mask): negative inside, positive outside.
    """
    inside = binary_mask
    outside = ~binary_mask

    dist_outside = distance_transform_edt(outside, sampling=(dy, dx))
    dist_inside = distance_transform_edt(inside, sampling=(dy, dx))

    signed_dist = dist_outside - dist_inside
    return signed_dist


def reinitialize_phasefield_from_levelset(signed_dist, eps_new):
    """Given a signed distance field 'signed_dist' and a desired interface width
    eps_new, reconstruct a smooth phase-field using a tanh profile:

        phi(x) = 0.5 * (1 - tanh( s / (2*eps_new) ))

    where s < 0 is "inside" (phi ~ 1), s > 0 is "outside" (phi ~ 0).
    """
    return 0.5 * (1.0 - np.tanh(signed_dist / (2.0 * eps_new)))


# -----------------------------------------------------------------------------
# VTK writer + checker
# -----------------------------------------------------------------------------

def save_vtk_image(filename, ice, sediment, temp, Lx, Ly, upsample=1):
    """Write a STRUCTURED_POINTS VTK file with 3 scalar fields:
       IcePhase, SedPhase, Temperature.

    ice, sediment, temp: 2D arrays (Ny, Nx)
    Lx, Ly: physical domain sizes
    upsample: optional upsampling factor for nicer visualization (default 1)
    """
    Ny, Nx = ice.shape

    if upsample > 1:
        zoom_factor = float(upsample)
        ice_vtk = zoom(ice, zoom_factor, order=1, mode="nearest")
        sediment_vtk = zoom(sediment, zoom_factor, order=1, mode="nearest")
        temp_vtk = zoom(temp, zoom_factor, order=1, mode="nearest")
        Ny_vtk, Nx_vtk = ice_vtk.shape
    else:
        ice_vtk = ice
        sediment_vtk = sediment
        temp_vtk = temp
        Ny_vtk, Nx_vtk = Ny, Nx

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
        f.write("Reinitialized Phase Field Output (no temp in input)\n")
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

    dims_line = next(l for l in lines if l.startswith("DIMENSIONS"))
    _, nx, ny, nz = dims_line.split()
    nx, ny, nz = int(nx), int(ny), int(nz)

    pd_line = next(l for l in lines if l.startswith("POINT_DATA"))
    _, npts = pd_line.split()
    npts = int(npts)

    print(f"DIMENSIONS: {nx} x {ny} x {nz}, POINT_DATA: {npts}")
    assert npts == nx * ny * nz, "POINT_DATA != DIMENSIONS product!"

    i = lines.index(pd_line) + 1
    while i < len(lines):
        if lines[i].startswith("SCALARS"):
            name = lines[i].split()[1]
            i += 1
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


# -----------------------------------------------------------------------------
# High-level driver: input has NO temp, output DOES
# -----------------------------------------------------------------------------

def reinitialize_ice_phase_notemp(
    infile,
    outfile,
    Nx_in,
    Ny_in,
    Lx,
    Ly,
    eps_new,
    threshold=0.5,
    dof_in=2,
    ice_idx_in=0,
    sed_idx_in=1,
    dof_out=3,
    ice_idx_out=0,
    sed_idx_out=1,
    temp_idx_out=2,
    T0_out=0.0,
    comm=PETSc.COMM_SELF,
    plot=False,
    smooth_before_threshold=False,
    sigma_threshold=0.75,
    smooth_signed_dist=True,
    sigma_sdf=2.0,
    Nx_out=None,
    Ny_out=None,
):
    """
    High-level driver for the case where the *input* Vec has only:
        ph1 (ice), ph2 (sediment)
    and the *output* Vec should have:
        IcePhase (reinitialized), SedPhase (copied), Temperature (constant T0_out).

    Parameters
    ----------
    infile : str
        Path to input PETSc Vec file (binary).
    outfile : str
        Path to output PETSc Vec file (binary).
    Nx_in, Ny_in : int
        Input grid dimensions (must match infile).
    Lx, Ly : float
        Physical domain size in x and y.
    eps_new : float
        Desired interface width for reconstructed phase-field.
    threshold : float
        Threshold for binarizing the ice phase (default 0.5).
    dof_in : int
        Number of DOFs in the input Vec (default 2).
    dof_out : int
        Number of DOFs in the output Vec (default 3).
    T0_out : float
        Constant temperature value for the output Temperature field.
    Nx_out, Ny_out : int or None
        Optional output grid size; if None, use input size.
    """
    print(f"Loading Vec from {infile} ...")
    vec_in = load_petsc_vec(infile, comm=comm)

    print("Converting Vec to array (Ny_in, Nx_in, dof_in)...")
    arr_in = vec_to_array_2d(vec_in, Nx=Nx_in, Ny=Ny_in, dof=dof_in)

    # Extract existing fields
    ice = arr_in[:, :, ice_idx_in]
    sediment = arr_in[:, :, sed_idx_in]

    # Grid spacing for signed distance computation
    dx = Lx / max(Nx_in - 1, 1)
    dy = Ly / max(Ny_in - 1, 1)

    # Optional smoothing before threshold
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
    ice_new = np.clip(ice_new, 0.0, 1.0)

    # Decide output grid size
    if Nx_out is None:
        Nx_out = Nx_in
    if Ny_out is None:
        Ny_out = Ny_in

    # Downsample or keep as is
    if (Nx_out != Nx_in) or (Ny_out != Ny_in):
        print(f"Downsampling from ({Nx_in}, {Ny_in}) to ({Nx_out}, {Ny_out}) ...")

        zoom_y = (Ny_out - 1) / max(Ny_in - 1, 1)
        zoom_x = (Nx_out - 1) / max(Nx_in - 1, 1)

        ice_use = zoom(ice_new, (zoom_y, zoom_x), order=1, mode="nearest")
        sed_use = zoom(sediment, (zoom_y, zoom_x), order=1, mode="nearest")
    else:
        ice_use = ice_new
        sed_use = sediment

    # Build output array with 3 DOFs
    arr_out = np.zeros((Ny_out, Nx_out, dof_out), dtype=arr_in.dtype)

    ny_i, nx_i = ice_use.shape
    ny_s, nx_s = sed_use.shape

    ny_min = min(Ny_out, ny_i, ny_s)
    nx_min = min(Nx_out, nx_i, nx_s)

    # Place ice and sediment
    arr_out[:ny_min, :nx_min, ice_idx_out] = ice_use[:ny_min, :nx_min]
    arr_out[:ny_min, :nx_min, sed_idx_out] = sed_use[:ny_min, :nx_min]

    # Temperature field: uniform T0_out everywhere
    arr_out[:, :, temp_idx_out] = T0_out

    print(f"Ice (new)     range: [{ice_use.min():.3e}, {ice_use.max():.3e}]")
    print(f"Sediment      range: [{sed_use.min():.3e}, {sed_use.max():.3e}]")
    print(f"Temperature   uniform: {T0_out:.3e}")
    print("Converting array back to Vec ...")
    vec_out = array_2d_to_vec(arr_out, comm=comm)

    print(f"Saving reinitialized Vec to {outfile} ...")
    save_petsc_vec(vec_out, outfile, comm=comm)

    vtkfile = outfile.replace(".dat", ".vtk")
    print(f"Saving VTK file to {vtkfile} ...")
    save_vtk_image(
        vtkfile,
        arr_out[:, :, ice_idx_out],
        arr_out[:, :, sed_idx_out],
        arr_out[:, :, temp_idx_out],
        Lx,
        Ly,
        upsample=1,
    )

    # Optional quick diagnostics plot
    if plot:
        fig, axs = plt.subplots(1, 3, figsize=(12, 4))
        im0 = axs[0].imshow(ice, origin="lower")
        axs[0].set_title("Ice (input)")
        plt.colorbar(im0, ax=axs[0], shrink=0.8)

        im1 = axs[1].imshow(ice_new, origin="lower")
        axs[1].set_title("Ice (reinitialized, fine grid)")
        plt.colorbar(im1, ax=axs[1], shrink=0.8)

        im2 = axs[2].imshow(arr_out[:, :, ice_idx_out], origin="lower")
        axs[2].set_title("Ice (output grid)")
        plt.colorbar(im2, ax=axs[2], shrink=0.8)

        plt.tight_layout()
        plt.show()

    vec_in.destroy()
    vec_out.destroy()
    print("Done.")


# -----------------------------------------------------------------------------
# Example main (edit paths + parameters as needed)
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    infile = (
        "/Users/jacksonbaglino/PetIGA-3.20/projects/permafrost/"
        "inputs/tests/SthavisthaInputs/AssemblePhase__first.dat"
    )
    outfile = (
        "/Users/jacksonbaglino/PetIGA-3.20/projects/permafrost/"
        "inputs/tests/SthavisthaInputs/AssemblePhase__first_reinit.dat"
    )

    # These must match the original input Vec
    Nx_in = 217 + 1
    Ny_in = 434 + 1
    Lx = 4.04e-4
    Ly = 8.08e-4

    # Output grid (can be the same as input or coarser)
    Nx_out = 434 + 1
    Ny_out = 868 + 1

    eps_new = 9.3e-7
    threshold = 0.5
    T0_out = 0.0  # or e.g. -20.0 if you want a physical temperature

    reinitialize_ice_phase_notemp(
        infile=infile,
        outfile=outfile,
        Nx_in=Nx_in,
        Ny_in=Ny_in,
        Lx=Lx,
        Ly=Ly,
        eps_new=eps_new,
        threshold=threshold,
        dof_in=2,
        ice_idx_in=0,
        sed_idx_in=1,
        dof_out=3,
        ice_idx_out=0,
        sed_idx_out=1,
        temp_idx_out=2,
        T0_out=T0_out,
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