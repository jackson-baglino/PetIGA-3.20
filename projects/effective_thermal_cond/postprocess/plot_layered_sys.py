#!/usr/bin/env python3
# ================================================
# Thermal Simulation Post-Processing Script (SVG + Custom Analytical)
# ================================================

import cmocean
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata

# -------------------------------------
# Domain and resolution settings
# -------------------------------------
Lx, Ly = 1.0, 1.0
# N = 2**8
# N = 60
N = 550
Nx, Ny = N + 1, N + 1

# -------------------------------------
# Plot control flags
# -------------------------------------
plot_profile_flag     = 1
plot_raw_curve        = 1
plot_analytical_curve = 1
plot_fitted_curve     = 0   # enable piecewise linear fit overlay

legend_outside        = 0  # Set to False to keep legend inside plot

# -------------------------------------
# Analytical slopes and intercepts
# -------------------------------------
k_i = 2.29
k_a = 0.02
q_out = -0.1
# q_out = 1.0
T_in = 273.15 - 30

interface_y = Ly / 2.0
# interface_y = 0

#
# Slopes:
m_ice = -q_out / k_i
m_air = -q_out / k_a

# Corrected intercepts to satisfy T_top = T_in and continuity at interface:
T_top = T_in
T_interface = T_top - m_air * (Ly - interface_y)

b_air = T_top - m_air * Ly
b_ice = T_interface - m_ice * interface_y

print("Using analytical piecewise solution:")
print(f"Ice region (bottom):     T(y) = {m_ice:.4f} * y + {b_ice:.4f}")
print(f"Air region (top):        T(y) = {m_air:.4f} * y + {b_air:.4f}")

# -------------------------------------
# Matplotlib styling
# -------------------------------------
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

# -------------------------------------
# Utility Functions
# -------------------------------------

def load_temperature(file_path, Nx, Ny):
    temperature = np.fromfile(file_path, dtype=">f8")
    while temperature.size > Nx * Ny:
        temperature = temperature[1:]
    return temperature.reshape((Ny, Nx))


def load_ice_field(file_path):
    data = np.loadtxt(file_path)
    x, y, _, ice = data[:,0], data[:,1], data[:,2], data[:,3]
    return x, y, ice


def plot_temperature_field(temperature, save_dir):
    fig, ax = plt.subplots(figsize=(7,6))
    norm = mcolors.Normalize(vmin=np.min(temperature), vmax=np.max(temperature))
    im = ax.imshow(
        temperature,
        cmap=cmocean.cm.thermal,
        origin="lower",
        extent=[
            -0.5 * Lx * 1e3 / (Nx - 1), Lx * 1e3 + 0.5 * Lx * 1e3 / (Nx - 1),
            -0.5 * Ly * 1e3 / (Ny - 1), Ly * 1e3 + 0.5 * Ly * 1e3 / (Ny - 1)
        ],
        aspect="equal",
        norm=norm,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Temperature (K)")
    ax.set(xlabel="X (mm)", ylabel="Y (mm)", title="Temperature Field")
    ax.grid(False)
    plt.savefig(os.path.join(save_dir, "temperature_field_TS.svg"), bbox_inches="tight")
    plt.close()
    print("Saved: temperature_field_TS.svg")


def plot_ice_field(x, y, ice, save_dir):
    grid_x = np.linspace(0, Lx, Nx, endpoint=True)
    grid_y = np.linspace(0, Ly, Ny, endpoint=True)
    grid_X, grid_Y = np.meshgrid(grid_x, grid_y)
    ice_grid = griddata((x, y), ice, (grid_X, grid_Y), method='linear', fill_value=np.nan, rescale=True)

    fig, ax = plt.subplots(figsize=(7,6))
    norm = mcolors.Normalize(vmin=np.nanmin(ice_grid), vmax=np.nanmax(ice_grid))
    im = ax.imshow(
        ice_grid,
        cmap=cmocean.cm.ice,
        origin="lower",
        extent=[
            -0.5 * Lx * 1e3 / (Nx - 1), Lx * 1e3 + 0.5 * Lx * 1e3 / (Nx - 1),
            -0.5 * Ly * 1e3 / (Ny - 1), Ly * 1e3 + 0.5 * Ly * 1e3 / (Ny - 1)
        ],
        aspect="equal",
        norm=norm,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Ice Phase Field")
    ax.set(xlabel="X (mm)", ylabel="Y (mm)", title="Ice Phase Field")
    ax.grid(False)
    plt.savefig(os.path.join(save_dir, "ice_field_TS.svg"), bbox_inches="tight")
    plt.close()
    print("Saved: ice_field_TS.svg")


def analytical_piecewise(y_coords, interface_y):
    analytical = np.zeros_like(y_coords)
    for i, y in enumerate(y_coords):
        if y >= interface_y:
            analytical[i] = m_air * y + b_air
        else:
            analytical[i] = m_ice * y + b_ice
    return analytical


def plot_centerline_profile(y_coords, temp_profile, analytical, save_dir):
    fig_profile, ax_profile = plt.subplots(figsize=(6,5))
    if plot_raw_curve:
        ax_profile.plot(y_coords * 1e3, temp_profile, label="Simulation", color="tab:purple")
    if plot_analytical_curve:
        ax_profile.plot(y_coords * 1e3, analytical, linestyle=":", color="tab:blue", label="Analytical")
    ax_profile.set(xlabel="Vertical Position (mm)", ylabel="Temperature (K)", title="Centerline Temperature Profile")
    ax_profile.grid(True, linestyle="--", alpha=0.4)

    if legend_outside:
        ax_profile.legend(loc="best", bbox_to_anchor=(1.0, 0.5), frameon=False)
        fig_profile.tight_layout(rect=[0, 0, 0.8, 1])
    else:
        ax_profile.legend(loc="best", frameon=False)

    plt.savefig(os.path.join(save_dir, "temperature_profile_only.svg"), bbox_inches="tight")
    plt.close()
    print("Saved: temperature_profile_only.svg")

    fit = None
    # --- piecewise linear fit across interface_y ---
    if plot_fitted_curve:
        mask_bot = y_coords < interface_y
        mask_top = ~mask_bot
        fit = np.empty_like(temp_profile)
        if mask_bot.any():
            m1, b1 = np.polyfit(y_coords[mask_bot], temp_profile[mask_bot], 1)
            fit[mask_bot] = m1 * y_coords[mask_bot] + b1
            print(f"Fitted bottom: T(y) = {m1:.4f} * y + {b1:.4f}")
        if mask_top.any():
            m2, b2 = np.polyfit(y_coords[mask_top], temp_profile[mask_top], 1)
            fit[mask_top] = m2 * y_coords[mask_top] + b2
            print(f"Fitted top:    T(y) = {m2:.4f} * y + {b2:.4f}")

    fig, ax1 = plt.subplots()
    fig.set_size_inches(6, 6)
    error = np.abs(analytical - temp_profile)

    ax1.plot(y_coords * 1e3, error, linestyle=":", color="tab:red", label="Error")
    ax1.set_xlabel("Vertical Position (mm)")
    ax1.set_ylabel("Error (K)", color="tab:red")
    ax1.tick_params(axis='y', labelcolor="tab:red")
    ax1.grid(True, linestyle="--", alpha=0.4)
    ax1.set_title("Centerline Profile and Error")

    ax2 = ax1.twinx()
    if plot_raw_curve:
        ax2.plot(y_coords * 1e3, temp_profile, label="Simulation", color="tab:purple")
    if plot_analytical_curve:
        ax2.plot(y_coords * 1e3, analytical, linestyle=":", color="tab:blue", label="Analytical")
    if plot_fitted_curve and fit is not None:
        ax2.plot(y_coords * 1e3, fit, linestyle="--", color="green", label="Fitted")
    ax2.set_ylabel("Temperature (K)", color="tab:purple")
    ax2.tick_params(axis='y', labelcolor="tab:purple")

    if legend_outside:
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, loc="best", bbox_to_anchor=(1.0, 0.5), frameon=False)
        fig.tight_layout(rect=[0, 0, 0.8, 1])
    else:
        ax1.legend(loc="best", frameon=False)

    plt.savefig(os.path.join(save_dir, "temperature_profile_centerline.svg"), bbox_inches="tight")
    plt.close()
    print("Saved: temperature_profile_centerline.svg")

# -------------------------------------
# Main Execution
# -------------------------------------

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <folder_name>")
        sys.exit(1)

    folder = sys.argv[1]
    base_dir = f"/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/{folder}/"
    save_dir = os.path.join(base_dir, "plots")
    os.makedirs(save_dir, exist_ok=True)

    temp_file = os.path.join(base_dir, "temperature.bin")
    ice_file = os.path.join(base_dir, "ice_data.dat")

    temperature = load_temperature(temp_file, Nx, Ny)
    x, y, ice = load_ice_field(ice_file)

    plot_temperature_field(temperature, save_dir)
    plot_ice_field(x, y, ice, save_dir)

    max_analytical = np.max(analytical_piecewise(np.array([0, Ly]), interface_y)) - 273.15
    min_analytical = np.min(analytical_piecewise(np.array([0, Ly]), interface_y)) - 273.15

    max_simulation = np.max(temperature) - 273.15
    min_simulation = np.min(temperature) - 273.15

    print(f"The temperature range from the analytical solution is: [{min_analytical:.4f}, {max_analytical:.4f}]")
    print(f"The temperature range from the simulation is: [{min_simulation:.4f}, {max_simulation:.4f}]")

    if plot_profile_flag:
        mid_x = Nx // 2
        # Use vertical coordinate vector spanning 0 to Ly inclusive, matching grid points
        y_coords = np.linspace(0, Ly, Ny)
        temp_profile = temperature[:, mid_x]
        analytical = analytical_piecewise(y_coords, interface_y)
        plot_centerline_profile(y_coords, temp_profile, analytical, save_dir)

if __name__ == "__main__":
    main()
