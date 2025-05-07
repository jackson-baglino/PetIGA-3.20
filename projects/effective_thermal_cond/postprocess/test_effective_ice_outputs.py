# ================================================
# 🔥 Thermal Simulation Post-Processing Script
# Usage: python test_effective_ice_outputs.py <folder_name>
# Goals: (1) Plot 2D temperature field, (2) Best midline fit, (3) Report midline error
# ================================================

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import Polynomial

# ---------------------------
# ⚙️ Global Configurations
# ---------------------------
# Lx, Ly = 0.5e-3, 0.5e-3
Lx, Ly = 1.0, 1.0
N = 128
Nx, Ny = N + 1, N + 1
k = 2.29         # Thermal conductivity (W/m/K)
q_out = 2 * 15.0   # Outgoing flux (W/m²) (linear case)
# q_out = 0.0      # Outgoing flux (W/m²) (quadratic case)
T_in = 240.15    # Inlet temperature (K) 
s = 0.0          # Internal source term (W/m³) (linear case)
# s = 100.0 / Lx   # Internal source term (W/m³) (quadratic case)

plot_profile_flag = True

# Configure matplotlib styles globally
plt.style.use("seaborn-v0_8-whitegrid")
import matplotlib as mpl
mpl.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "lines.linewidth": 2,
    "figure.dpi": 150
})

# ---------------------------
# 📂 Utility Functions
# ---------------------------

def load_temperature(file_path, Nx, Ny):
    """Load and reshape temperature data."""
    temperature = np.fromfile(file_path, dtype=">f8")
    while temperature.size > Nx * Ny:
        temperature = temperature[1:]
    return temperature.reshape((Ny, Nx))

import matplotlib.colors as mcolors

def plot_temperature_field(temperature, save_dir):
    """Plot and save the full 2D temperature field with absolute colorbar values."""
    temp_min = np.min(temperature)
    temp_max = np.max(temperature)

    fig, ax = plt.subplots(figsize=(7, 6))

    norm = mcolors.Normalize(vmin=temp_min, vmax=temp_max)
    im = ax.imshow(temperature, cmap="magma", origin="lower", interpolation="nearest",
                   extent=[0, Lx*1e3, 0, Ly*1e3], aspect="equal", norm=norm)  # pass norm explicitly!

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Temperature (K)")
    cbar.ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.4f'))  # Format tick labels normally, no scientific notation

    ax.set(xlabel="X (mm)", ylabel="Y (mm)", title="Temperature Field")
    output_path = os.path.join(save_dir, "temperature_field_TS.png")
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved 2D temperature plot: {output_path}")

def fit_temperature_profile(y_coords, temp_profile, threshold=1e-1):
    """Fit linear and quadratic models and select the best fit intelligently."""
    # Fit both
    linear_coefs = Polynomial.fit(y_coords, temp_profile, deg=1).convert().coef
    quad_coefs = Polynomial.fit(y_coords, temp_profile, deg=2).convert().coef

    linear_fit = linear_coefs[1] * y_coords + linear_coefs[0]
    quad_fit = quad_coefs[2] * y_coords**2 + quad_coefs[1] * y_coords + quad_coefs[0]

    rss_linear = np.sum((temp_profile - linear_fit)**2)
    rss_quad = np.sum((temp_profile - quad_fit)**2)

    # Calculate relative importance of quadratic term
    quad_strength = abs(quad_coefs[2]) / (abs(linear_coefs[1]) + 1e-12)  # add tiny number to avoid div by zero

    # Decide
    if (rss_quad < rss_linear) and (quad_strength > threshold):
        # Quadratic fit is both better and significant
        return quad_fit, "Quadratic Fit", "crimson", quad_coefs
    else:
        # Prefer linear fit
        return linear_fit, "Linear Fit", "darkorange", linear_coefs

def analytical_profile(y_coords):
    """Compute analytical temperature profile."""
    a = s / 2 / k
    b = q_out / k
    c = T_in - b * Ly - a * Ly**2
    return a * y_coords**2 + b * y_coords + c, (a, b, c)

def plot_centerline_profile(y_coords, temp_profile, best_fit, best_label, best_color,
                             analytical, save_dir, fit_coefs, ana_coefs):
    """Plot the temperature centerline profile and fits."""
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(y_coords*1e3, temp_profile, label="Raw Profile", color="black")
    ax.plot(y_coords*1e3, best_fit, label=best_label, linestyle="--", color=best_color)
    if s == 0.0:
        ax.plot(y_coords*1e3, analytical, label="Analytical (linear)", linestyle=":", color="blue")
    else:
        ax.plot(y_coords*1e3, analytical, label="Analytical (quadratic)", linestyle=":", color="blue")
    ax.set(xlabel="Vertical Position (mm)", ylabel="Temperature (K)", 
           title="Centerline Temperature Profile")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="best", frameon=False, fontsize=12)

    profile_path = os.path.join(save_dir, "temperature_profile_centerline.png")
    plt.savefig(profile_path, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved centerline profile plot: {profile_path}")

    # Report fits
    if len(fit_coefs) == 3:
        fit_eqn = f"Fit: y = {fit_coefs[2]:.2e}·x² + {fit_coefs[1]:.2e}·x + {fit_coefs[0]:.2f}"
    else:
        fit_eqn = f"Fit: y = {fit_coefs[1]:.2e}·x + {fit_coefs[0]:.2f}"
    ana_eqn = f"Analytical: y = {ana_coefs[0]:.2e}·x² + {ana_coefs[1]:.2e}·x + {ana_coefs[2]:.2f}"

    print(fit_eqn)
    print(ana_eqn)

def plot_error_profile(y_coords, temp_profile, analytical, fit,save_dir):
    """Plot error between raw data and analytical solution."""
    error_analytical = (temp_profile - analytical) / (np.max(temp_profile) - np.min(temp_profile)) * 100
    error_fit = (temp_profile - fit) / (np.max(temp_profile) - np.min(temp_profile)) * 100
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(error_analytical, y_coords*1e3, color="darkgreen", label="Error (Raw - Analytical)")
    ax.plot(error_fit, y_coords*1e3, color="crimson", label="Error (Raw - Fit)")
    ax.set(xlabel="Temperature Difference (%)", ylabel="Vertical Position (mm)",
           title="Centerline Error Profile")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="best", frameon=False, fontsize=12)

    error_path = os.path.join(save_dir, "temperature_profile_error.png")
    plt.savefig(error_path, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved error plot: {error_path}")

# ---------------------------
# 🚀 Main Execution
# ---------------------------
def main():
    if len(sys.argv) != 2:
        print("Usage: python test_effective_ice_outputs.py <folder_name>")
        sys.exit(1)

    folder = sys.argv[1]
    base_dir = f"/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/{folder}/"
    save_dir = os.path.join(base_dir, "plots")
    os.makedirs(save_dir, exist_ok=True)

    # Load data
    temp_file = os.path.join(base_dir, "temperature.bin")
    temperature = load_temperature(temp_file, Nx, Ny)

    # Plot temperature field
    plot_temperature_field(temperature, save_dir)

    # Midline profile analysis
    if plot_profile_flag:
        mid_x = Nx // 2
        y_coords = np.linspace(0, Ly, Ny)
        temp_profile = temperature[:, mid_x]

        best_fit, best_label, best_color, fit_coefs = fit_temperature_profile(y_coords, temp_profile)
        analytical, ana_coefs = analytical_profile(y_coords)

        plot_centerline_profile(y_coords, temp_profile, best_fit, best_label, best_color,
                                analytical, save_dir, fit_coefs, ana_coefs)
        plot_error_profile(y_coords, temp_profile, analytical, best_fit,save_dir)

if __name__ == "__main__":
    main()