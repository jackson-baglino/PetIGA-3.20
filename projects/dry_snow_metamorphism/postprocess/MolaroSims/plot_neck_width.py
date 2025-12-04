import matplotlib.pyplot as plt
from pathlib import Path
import re
from typing import Tuple, Union

import numpy as np
import pyvista as pv
from skimage.morphology import binary_opening, binary_closing, disk
from skimage.measure import label, regionprops, find_contours
from scipy.ndimage import zoom


# ============================================================
# USER PARAMETERS / FLAGS
# ============================================================

# Interpolation factor for mid-plane slice resolution
INTERP_FACTOR = 4.0

# Threshold for binarizing the IcePhase field
ICE_THRESHOLD = 0.5

# Minimum connected-component area to be considered a grain
MIN_REGION_AREA = 50

# Require fused grains to measure a neck
REQUIRE_FUSION = True

# Toggle diagnostic slice plotting during neck-width evaluation
DISPLAY_SLICE_PLOTS = True

# Toggle saving midplane cross-section plots
SAVE_SLICE_PLOTS = True


STYLE_PATH = Path(__file__).with_name("plot_style.mplstyle")
plt.style.use(STYLE_PATH)

def read_fields_from_vtk(
    filename: Union[str, Path],
    phi_name: str = "IcePhase",
    T_name: str = "Temperature",
    rhov_name: str = "VaporDensity",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read three 3D fields (phi, T, rhov) from a VTK file written by the solver.

    Parameters
    ----------
    filename : str or Path
        Path to the VTK file, e.g. 'vtkOut/solV_00010.vtk'.
    phi_name : str
        Name of the ice phase field in the VTK point_data.
    T_name : str
        Name of the temperature field in the VTK point_data.
    rhov_name : str
        Name of the vapor density field in the VTK point_data.

    Returns
    -------
    phi : ndarray, shape (nx, ny, nz)
        Ice phase field.
    T : ndarray, shape (nx, ny, nz)
        Temperature field.
    rhov : ndarray, shape (nx, ny, nz)
        Vapor density field.

    Notes
    -----
    The returned arrays are shaped according to the VTK grid dimensions
    (nx, ny, nz), where nx, ny, nz are the number of points in the x, y,
    and z directions respectively. If you prefer a different axis order,
    you can transpose them after calling this function.
    """
    filename = Path(filename)
    grid = pv.read(filename)

    # VTK UniformGrid / ImageData: dimensions is (nx, ny, nz)
    nx, ny, nz = grid.dimensions

    def _get_field(name: str) -> np.ndarray:
        if name not in grid.point_data:
            raise KeyError(
                f"Field '{name}' not found in {filename}. "
                f"Available fields: {list(grid.point_data.keys())}"
            )
        data = grid.point_data[name]
        if data.size != nx * ny * nz:
            raise ValueError(
                f"Field '{name}' size {data.size} does not match nx*ny*nz = {nx*ny*nz}"
            )
        # VTK uses x-fastest (Fortran-like) ordering for point data
        arr = data.reshape((nx, ny, nz), order="F")
        return arr

    phi = _get_field(phi_name)
    T = _get_field(T_name)
    rhov = _get_field(rhov_name)

    return phi, T, rhov

def compute_neck_width_from_slice(
    phi_slice: np.ndarray,
    pixel2length: float,
    threshold: float = 0.5,
    min_region_area: int = 50,
    require_fusion: bool = True,
    disp: bool = False,
) -> float:
    """
    Compute neck width from a 2D mid-plane slice of the ice phase field,
    using the phi = threshold contour (continuous interface).
    The neck width is defined as the minimum horizontal span (width) of the contour within a central vertical band where the grains are actually in contact. We estimate this band by taking a neighborhood around the vertical center-of-mass of the fused ice mask (so the band follows small vertical drift through time) and then, for each horizontal cross-section (y band) in that region, finding the leftmost and rightmost contour intersections. The minimum such span in that central band is taken as the neck width. This avoids measuring narrow features at the very top/bottom of the grains and focuses on the physical neck near the center.

    Parameters
    ----------
    phi_slice : ndarray, shape (nx, ny)
        Ice phase field on a 2D slice.
    pixel2length : float
        Conversion factor from pixels to physical length (meters per pixel).
    threshold : float
        Contour level for the ice interface, e.g. 0.5.
    min_region_area : int
        Minimum number of pixels to consider a 'big' grain component when
        checking whether grains have fused.
    require_fusion : bool
        If True, returns 0 until the two main grains have merged into
        one connected component (based on a binary threshold).
    disp : bool
        If True, show a diagnostic plot of the slice and neck measurement.

    Returns
    -------
    neck_length : float
        Neck width in physical units (meters). Returns 0 if not yet fused
        (when require_fusion=True) or if the geometry is degenerate.
    """
    phi_img = phi_slice.T  # shape (ny, nx), where ny is vertical (rows), nx is horizontal (cols)

    # 1) Fusion check using a binary mask
    mask = phi_img >= threshold
    mask = binary_opening(mask, disk(1))
    mask = binary_closing(mask, disk(1))

    labels_img = label(mask, connectivity=2)
    props = regionprops(labels_img)

    # Keep only "large" components, ignore tiny noise
    large_props = [p for p in props if p.area >= min_region_area]

    if not large_props:
        return 0.0

    if require_fusion:
        # Sort by area, descending
        big_props = sorted(large_props, key=lambda p: p.area, reverse=True)

        if len(big_props) >= 2:
            # Two separate grains: treat as no neck yet
            if disp:
                print("Two large components detected; treating neck width as 0 and plotting slice.")
                plt.figure()
                plt.imshow(phi_img, origin="lower", aspect="equal")
                cbar = plt.colorbar()
                cbar.set_label("IcePhase", fontweight="bold")
                plt.contour(mask, levels=[0.5], colors="w", linewidths=0.8)
                plt.title("Mid-plane slice (grains not yet fused)", fontweight="bold")
                plt.xlabel("x (pixels)")
                plt.ylabel("y (pixels)")
                plt.tight_layout()
            return 0.0

    # 2) Extract phi = threshold contours (continuous interface)
    contours = find_contours(phi_img, level=threshold)
    if len(contours) == 0:
        # No interface found
        return 0.0

    # Concatenate all contour segments into one array of points
    # Each contour is an array of shape (N_i, 2) with (row, col) = (y, x)
    all_points = np.vstack(contours)

    # Scan across y to find the narrowest horizontal span of the contour,
    # but only within a central vertical band where the grains are expected
    # to be in contact (near the middle of the image). The band is centered
    # at the vertical center-of-mass of the fused ice mask and follows small
    # vertical drift through time.
    ny, nx = phi_img.shape

    ys = all_points[:, 0]
    xs = all_points[:, 1]

    # Vertical extent of the contour as a whole
    y_all_min, y_all_max = ys.min(), ys.max()
    if y_all_max <= y_all_min:
        return 0.0

    # Compute vertical center-of-mass of the fused ice mask
    mask_y, mask_x = np.nonzero(mask)
    if mask_y.size == 0:
        return 0.0
    com_y = mask_y.mean()

    # Define a central band around com_y (e.g., +/- 25% of image height)
    band_half_height = 0.25 * ny
    band_min = max(y_all_min, com_y - band_half_height)
    band_max = min(y_all_max, com_y + band_half_height)
    if band_max <= band_min:
        # Fallback: use full contour extent if something goes wrong
        band_min, band_max = y_all_min, y_all_max

    # Use a set of horizontal bands spanning only this central vertical region.
    # Choose the number of bands proportional to the band height, but not more
    # than ny to keep a ~1-pixel vertical resolution.
    band_height = band_max - band_min
    n_bins = int(min(ny, max(10, np.ceil(band_height))))
    if n_bins <= 0:
        return 0.0

    bins = np.linspace(band_min, band_max, n_bins + 1)

    best_width = None
    best_x_left = None
    best_x_right = None
    best_y = None

    for i_bin in range(n_bins):
        y0 = 0.5 * (bins[i_bin] + bins[i_bin + 1])
        mask_band = (ys >= bins[i_bin]) & (ys < bins[i_bin + 1])
        if mask_band.sum() < 2:
            continue

        x_vals = xs[mask_band]
        width = x_vals.max() - x_vals.min()
        if best_width is None or width < best_width:
            best_width = width
            best_x_left = x_vals.min()
            best_x_right = x_vals.max()
            best_y = y0

    # If we could not find any valid band, return 0
    if best_width is None:
        return 0.0

    neck_pixels = float(best_width)
    neck_length = neck_pixels * pixel2length

    # For plotting, define two points at the narrowest cross-section
    x1 = float(best_x_left)
    x2 = float(best_x_right)
    p1_y = float(best_y)
    p2_y = float(best_y)

    if disp:
        # Convert neck length to microns for annotation
        neck_um = neck_length * 1e6

        plt.figure()
        # phi_img is indexed in the standard (y, x) = (row, col) convention.
        # imshow(phi_img, origin="lower") uses x = column index, y = row index,
        # so we can plot contour coordinates (x, y) directly.
        plt.imshow(phi_img, origin="lower", aspect="equal")
        cbar = plt.colorbar()
        cbar.set_label("IcePhase", fontweight="bold")

        x1_plot, y1_plot = x1, p1_y
        x2_plot, y2_plot = x2, p2_y

        # Plot the two neck points and the segment between them
        plt.plot([x1_plot, x2_plot], [y1_plot, y2_plot],
                 "ro", markersize=5, label="Neck points")
        plt.plot([x1_plot, x2_plot], [y1_plot, y2_plot],
                 "r-", linewidth=2, label=f"Neck width ≈ {neck_um:.2f} μm")

        # Also overlay the phi=threshold contour(s) for context
        for c in contours:
            # c is (y, x); plot as (x, y) directly
            plt.plot(c[:, 1], c[:, 0], "w-", linewidth=0.8)

        plt.title("Mid-plane slice with contour-based neck measurement", fontweight="bold")
        plt.xlabel("x (pixels)")
        plt.ylabel("y (pixels)")
        plt.legend(loc="best", fontsize=8)
        plt.tight_layout()

    return neck_length


def compute_neck_width_from_midplane(
    phi3d: np.ndarray,
    pixel2length: float,
    threshold: float = 0.5,
    min_region_area: int = 50,
    require_fusion: bool = True,
    disp: bool = False,
    interp_factor: float = 1.0,
) -> float:
    """
    Convenience wrapper: take the mid-plane slice in z and compute neck width.

    Assumes phi3d has shape (nx, ny, nz), with the third axis being z.

    Parameters
    ----------
    phi3d : ndarray
        3D ice phase field with shape (nx, ny, nz).
    pixel2length : float
        Conversion factor from pixels to physical length (meters per pixel)
        on the original grid.
    threshold : float
        Threshold for ice vs non-ice, e.g. 0.5.
    min_region_area : int
        Minimum number of pixels to consider a 'big' grain component.
    require_fusion : bool
        If True, returns 0 until the two main grains have merged into
        one connected component.
    disp : bool
        If True, show a diagnostic plot of the slice, skeleton, and neck.
    interp_factor : float
        Optional interpolation factor to increase in-plane resolution.
        For example, interp_factor=4.0 will upsample the mid-plane slice
        by a factor of 4 in x and y using cubic interpolation, and adjust
        the pixel2length accordingly.
    """
    nx, ny, nz = phi3d.shape
    k_mid = nz // 2
    phi_slice = phi3d[:, :, k_mid]

    # Optionally upsample the mid-plane slice to improve geometric resolution.
    # This mirrors the MATLAB approach where an interp_factor is used.
    if interp_factor != 1.0 and interp_factor > 0.0:
        phi_slice_hi = zoom(phi_slice, zoom=interp_factor, order=3)
        pixel2length_hi = pixel2length / interp_factor
    else:
        phi_slice_hi = phi_slice
        pixel2length_hi = pixel2length

    return compute_neck_width_from_slice(
        phi_slice=phi_slice_hi,
        pixel2length=pixel2length_hi,
        threshold=threshold,
        min_region_area=min_region_area,
        require_fusion=require_fusion,
        disp=disp,
    )

def main():
    parent_dir = "/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/archived_results/" \
        "BestParams/Group1/NASAv2-Molaro2D_T-20.0_hum_2024-11-26__09.41.38"

    Lx = 0.0002424      
    Ly = 0.0003884
    Lz = 0.0002424

    Nx = 134
    Ny = 214
    Nz = 134

    # Directory containing VTK files
    vtk_dir = Path(parent_dir) / "vtkOut"

    # Collect all VTK files that match the solver naming pattern
    vtk_files = sorted(vtk_dir.glob("solV_*.vtk"))
    if not vtk_files:
        raise FileNotFoundError(f"No VTK files matching 'solV_*.vtk' found in {vtk_dir}")

    print(f"Found {len(vtk_files)} VTK files in {vtk_dir}")

    # Compute pixel-to-length conversion (match MATLAB convention)
    pixel2length = np.mean([Lx / Nx, Ly / Ny, Lz / Nz])
    print(f"pixel2length = {pixel2length:.6e} m/pixel")

    # Arrays to store results
    steps = []
    neck_lengths = []

    # Create folder for cross‑section plots
    cross_dir = Path(parent_dir) / "cross_sections"
    cross_dir.mkdir(exist_ok=True)

    # Loop over all VTK files and compute neck widths
    for i, f in enumerate(vtk_files, start=1):
        # Extract the numeric step from the filename, e.g. solV_00010.vtk -> 10
        m = re.search(r"\d+", f.stem)
        if m is None:
            print(f"Warning: could not parse step index from filename {f.name}; skipping.")
            continue
        step_idx = int(m.group())

        # Read fields from VTK; we only need phi for the neck computation
        phi, T, rhov = read_fields_from_vtk(f)

        # For debugging, print shape/min/max for the first file
        if i == 1:
            print(f"[DEBUG] {f.name}: phi shape = {phi.shape}, min = {phi.min():.3e}, max = {phi.max():.3e}")

        # Compute neck width from the mid-plane slice (no plotting in loop)
        neck_length = compute_neck_width_from_midplane(
            phi,
            pixel2length=pixel2length,
            threshold=ICE_THRESHOLD,
            min_region_area=MIN_REGION_AREA,
            require_fusion=REQUIRE_FUSION,
            disp=DISPLAY_SLICE_PLOTS,
            interp_factor=INTERP_FACTOR,
        )

        steps.append(step_idx)
        neck_lengths.append(neck_length)

        print(f"{f.name}: step {step_idx:6d}, neck width = {neck_length:.6e} m")

        if SAVE_SLICE_PLOTS:
            _ = compute_neck_width_from_midplane(
                phi,
                pixel2length=pixel2length,
                threshold=ICE_THRESHOLD,
                min_region_area=MIN_REGION_AREA,
                require_fusion=REQUIRE_FUSION,
                disp=DISPLAY_SLICE_PLOTS,
                interp_factor=INTERP_FACTOR,
            )
            fig_path = cross_dir / f"cross_section_{step_idx:06d}.png"
            plt.savefig(fig_path, dpi=200, bbox_inches="tight")
            plt.close()

    # Convert to numpy arrays and sort by step index
    steps = np.array(steps, dtype=int)
    neck_lengths = np.array(neck_lengths, dtype=float)

    sort_idx = np.argsort(steps)
    steps_sorted = steps[sort_idx]
    neck_sorted = neck_lengths[sort_idx]

    # Save results to CSV in the parent directory
    output_file = Path(parent_dir) / "neck_widths_vs_step.csv"
    data_to_save = np.column_stack([steps_sorted, neck_sorted])
    np.savetxt(
        output_file,
        data_to_save,
        delimiter=",",
        header="step,neck_width_m",
        comments="",
    )


    print(f"\nSaved neck widths to {output_file}")

    # ------------------------------------------------------------------
    # Load SSA_evo.dat to map steps to physical time and build time series
    # ------------------------------------------------------------------
    ssa_file = Path(parent_dir) / "SSA_evo.dat"
    if not ssa_file.exists():
        raise FileNotFoundError(f"SSA_evo.dat not found in {parent_dir}")

    # SSA_evo.dat has at least four columns; we use:
    #  - column 3 (index 2) as time
    #  - column 4 (index 3) as time step index (matching VTK step)
    ssa_data = np.loadtxt(ssa_file)
    if ssa_data.ndim == 1:
        ssa_data = ssa_data[None, :]  # ensure 2D

    time_col = ssa_data[:, 2]
    step_col = ssa_data[:, 3].astype(int)

    # Map each measured step to its corresponding time
    times_for_steps = np.full_like(steps_sorted, fill_value=np.nan, dtype=float)
    for i_step, step in enumerate(steps_sorted):
        mask = step_col == step
        if not np.any(mask):
            print(f"Warning: no matching time found in SSA_evo.dat for step {step}")
            continue
        # If multiple rows match, take the first
        times_for_steps[i_step] = time_col[mask][0]

    # Save time vs neck width to a new CSV
    output_time_file = Path(parent_dir) / "neck_widths_vs_time.csv"
    data_time_to_save = np.column_stack([times_for_steps, steps_sorted, neck_sorted])
    np.savetxt(
        output_time_file,
        data_time_to_save,
        delimiter=",",
        header="time_s,step,neck_width_m",
        comments="",
    )
    print(f"Saved neck widths vs time to {output_time_file}")

    # ------------------------------------------------------------------
    # Plot neck width vs time (in minutes) and save the figure
    # Only times corresponding to actual VTK/neck measurements are used.
    # ------------------------------------------------------------------
    valid_mask = ~np.isnan(times_for_steps)
    if np.any(valid_mask):
        t_sec = times_for_steps[valid_mask]
        w = neck_sorted[valid_mask]

        # Convert to minutes for plotting
        t_min = t_sec / 60.0

        # Sort by time in case step/time ordering differs
        sort_t_idx = np.argsort(t_min)
        t_min_sorted = t_min[sort_t_idx]
        w_sorted = w[sort_t_idx]

        fig, ax = plt.subplots()
        ax.plot(
            t_min_sorted,
            w_sorted,
            linestyle="-",
        )
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Neck width (m)")
        ax.set_title("Neck width evolution (VTK outputs only)")
        ax.grid(True)

        fig_path_time = Path(parent_dir) / "neck_width_vs_time.png"
        fig.savefig(fig_path_time, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved neck width vs time plot to {fig_path_time}")
    else:
        print("No valid time data found to plot neck width vs time.")

if __name__ == "__main__":
    main()