from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
#!/usr/bin/env python3
"""
Plot vector‑field components stored in binary file.

The binary file contains 2·Nx·Ny float64 values (big‑endian):
    t_x(0,0), t_y(0,0), t_x(1,0), t_y(1,0), …   (x varies fastest)

Usage
-----
    python plot_vector_field.py  <output_folder>  [Nx] [Ny] [Lx] [Ly]

If Nx,Ny are omitted the script assumes N=32 → Nx=Ny=32.
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmocean
from scipy.interpolate import griddata

# --- material conductivities (W/m·K) --------------------------
K_ICE = 2.29
K_AIR = 0.02

# ---------- helper -------------------------------------------------
def imsave(field, fname, Lx=1.0, Ly=1.0, vmin=None, vmax=None):
    """Save heat‑map; if vmin/vmax are None, use data limits."""
    ny, nx = field.shape
    fig, ax = plt.subplots(figsize=(6,5))
    if vmin is None: vmin = field.min()
    if vmax is None: vmax = field.max()
    im = ax.imshow(field,
                   origin="lower",
                   cmap=cmocean.cm.ice,
                   extent=[0, Lx, 0, Ly],
                   aspect="equal",
                   norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    plt.colorbar(im, ax=ax).set_label("Value")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.tight_layout(); plt.savefig(fname); plt.close()
    print("saved", fname)

# plot_ice_field helper
def plot_ice_field(datapath, nx, ny, save_path, Lx=1.0, Ly=1.0):
    """Read x y ice  (format: x y z ice) and plot heat‑map."""
    data = np.loadtxt(datapath)
    x, y, ice = data[:,0], data[:,1], data[:,-1]

    grid_x = np.linspace(0, Lx, nx, endpoint=True)
    grid_y = np.linspace(0, Ly, ny, endpoint=True)
    grid_X, grid_Y = np.meshgrid(grid_x, grid_y)

    ice_grid = griddata((x, y), ice, (grid_X, grid_Y),
                        method='linear', fill_value=np.nan)
    
    print("The range of values in the ice field is", ice.min(), ice.max())
    imsave(ice_grid, save_path, Lx, Ly, vmin=0.0, vmax=1.0)

def compute_keff(tx, ty, ice, Lx=1.0, Ly=1.0, order=4):
    """
    Estimate full k_eff tensor with global 2‑D Gauss‑Legendre quadrature.

    * tx, ty, ice are Ny×Nx arrays on a uniform grid.
    * order = number of Gauss points per direction (default 4).
    """
    # Fill NaNs in ice
    if np.isnan(ice).any():
        from scipy import ndimage
        mask = np.isnan(ice)
        ice = ice.copy()
        ice[mask] = ndimage.distance_transform_edt(mask,
                        return_distances=False,
                        return_indices=True)[0][mask]

    ny, nx = tx.shape
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)

    # Finite-difference gradients on grid
    dtx_dx = np.gradient(tx, dx, axis=1)  # ∂t_x/∂x
    dtx_dy = np.gradient(tx, dy, axis=0)  # ∂t_x/∂y
    dty_dx = np.gradient(ty, dx, axis=1)  # ∂t_y/∂x
    dty_dy = np.gradient(ty, dy, axis=0)  # ∂t_y/∂y

    # Material conductivity on grid
    k_grid = K_ICE * ice + K_AIR * (1.0 - ice)

    # Integrand values on grid
    f_xx_grid = k_grid * (dtx_dx + 1.0)
    f_yy_grid = k_grid * (dty_dy + 1.0)
    f_xy_grid = k_grid * dtx_dy
    f_yx_grid = k_grid * dty_dx

    # Build bilinear interpolators
    x_nodes = np.linspace(0, Lx, nx)
    y_nodes = np.linspace(0, Ly, ny)
    fx_interp = RegularGridInterpolator((y_nodes, x_nodes), f_xx_grid)
    fy_interp = RegularGridInterpolator((y_nodes, x_nodes), f_yy_grid)
    fxy_interp = RegularGridInterpolator((y_nodes, x_nodes), f_xy_grid)
    fyx_interp = RegularGridInterpolator((y_nodes, x_nodes), f_yx_grid)

    # Gauss-Legendre nodes and weights
    gxi, gwi = leggauss(order)
    x_gl = 0.5 * Lx * (gxi + 1.0)
    y_gl = 0.5 * Ly * (gxi + 1.0)
    wx = 0.5 * Lx * gwi
    wy = 0.5 * Ly * gwi

    # Tensor-product quadrature
    I_xx = 0.0
    I_yy = 0.0
    I_xy = 0.0
    I_yx = 0.0
    for j, y in enumerate(y_gl):
        for i, x in enumerate(x_gl):
            w = wx[i] * wy[j]
            pt = np.array([y, x])
            I_xx += w * fx_interp(pt).item()
            I_yy += w * fy_interp(pt).item()
            I_xy += w * fxy_interp(pt).item()
            I_yx += w * fyx_interp(pt).item()

    area = Lx * Ly
    keff_xx = I_xx / area
    keff_yy = I_yy / area
    keff_xy = I_xy / area
    keff_yx = I_yx / area

    return keff_xx, keff_yy, keff_xy, keff_yx

# ---------- main ---------------------------------------------------
def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python plot_vector_field.py <folder> [Nx] [Ny] [Lx] [Ly]")

    folder = sys.argv[1]
    Nx = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    Ny = int(sys.argv[3]) if len(sys.argv) > 3 else Nx
    Lx = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    Ly = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0

    base = f"./outputs/homog/{folder}"
    binfile = os.path.join(base, "t_vec.dat")
    if not os.path.isfile(binfile):
        sys.exit(f"{binfile} not found")

    # read raw big‑endian float64
    raw = np.fromfile(binfile, dtype=">f8")

    # Remove the first data point
    raw = raw[1:]
    expected = 2 * Nx * Ny
    if raw.size != expected:
        sys.exit(f"File length {raw.size} != 2*{Nx}*{Ny}={expected}")

    # de‑interleave
    t_x = raw[0::2].reshape(Ny, Nx)
    t_y = raw[1::2].reshape(Ny, Nx)

    # Print the mean values
    print("The mean values of t_x and t_y are", t_x.mean(), t_y.mean())

    # save images
    os.makedirs(os.path.join(base, "plots"), exist_ok=True)
    imsave(t_x, os.path.join(base, "plots", "t_x.svg"), Lx, Ly)
    imsave(t_y, os.path.join(base, "plots", "t_y.svg"), Lx, Ly)

    # --- optional ice field -------------------------------------------------
    icefile = os.path.join(base, "ice_data.dat")
    if os.path.isfile(icefile):
        ice_svg = os.path.join(base, "plots", "ice_field.svg")
        plot_ice_field(icefile, Nx, Ny, ice_svg, Lx, Ly)

    if os.path.isfile(icefile):
        ice_grid = np.loadtxt(icefile)[:,-1]    # quick 1‑D read
        # Re‑interpolate to grid matching t_x
        x_pts = np.loadtxt(icefile)[:,0]
        y_pts = np.loadtxt(icefile)[:,1]
        grid_x = np.linspace(0, Lx, Nx, endpoint=True)
        grid_y = np.linspace(0, Ly, Ny, endpoint=True)
        grid_X, grid_Y = np.meshgrid(grid_x, grid_y)
        ice_grid = griddata((x_pts, y_pts), ice_grid, (grid_X, grid_Y),
                            method='linear', fill_value=0.0)
        kxx, kyy, kxy, kyx = compute_keff(t_x, t_y, ice_grid, Lx, Ly)
        print(f"Estimated k_eff_xx = {float(kxx):.6f}, k_eff_yy = {float(kyy):.6f}, "
              f"k_eff_xy = {float(kxy):.6f}, k_eff_yx = {float(kyx):.6f}")
    else:
        print("No ice_data.dat — k_eff not computed.")

if __name__ == "__main__":
    main()
