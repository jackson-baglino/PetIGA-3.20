import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
Plot effective thermal conductivity (k_eff) vs. time and SSA.

This script pairs rows from `k_eff.csv` with SSA/time entries in `SSA_evo.dat`.
It does not scan `sol_*.dat`; instead it uses the number of rows in k_eff.csv
and matches by timestep index to SSA_evo.dat (with a small offset heuristic).

Outputs:
  - k_eff_vs_time.png
  - k_eff_vs_ssa.png
  - ssa_vs_time.png

Adjust the parameters below (folder_path, plotting options) as needed.
"""

# --- Parameters ---
folder_path = '/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch/grains_2D_Tm-15_hum95_tf28d__2025-09-16__09.34.10'
ssa_file = os.path.join(folder_path, 'SSA_evo.dat')
k_eff_file = os.path.join(folder_path, 'k_eff.csv')  # CSV file with header

# --- Plotting Options ---
use_scatter_for_ssa = True    # Set True for scatter plots, False for line plots
scatter_step = 1             # Plot every Nth point if using scatter
save_figures = True           # Save figures to disk
dpi = 600                     # DPI for saved figures
normalize_values = False       # If True, normalize SSA and k_eff by initial values

# --- Output Directory Setup ---
parent_folder_name = os.path.basename(os.path.normpath(folder_path))
output_dir = f'./outputs/homog/plots/{parent_folder_name}'
os.makedirs(output_dir, exist_ok=True)

# --- Step 1: Determine counts from outputs (no need to scan sol_*.dat) ---
# We'll use k_eff.csv row count for the number of k_eff outputs, and SSA_evo.dat for SSA steps.

# --- Step 2: Load SSA_evo.dat ---
ssa_data = np.loadtxt(ssa_file)
ssa_time = ssa_data[:, 2]        # time column
ssa_time_steps = ssa_data[:, 3]  # time step column
ssa_values = ssa_data[:, 0]      # SSA column (1st column)

# --- Step 3: Load k_eff.csv ---
k_eff_df = pd.read_csv(k_eff_file)
k_xx = k_eff_df['k_00'].values
k_yy = k_eff_df['k_11'].values

# Determine how many k_eff outputs we have
n_keff = len(k_xx)

# --- Step 4: Match by index using counts from k_eff and SSA_evo ---
# Normalize SSA time steps to integers in case they were stored as floats
ssa_steps_int = np.rint(ssa_time_steps).astype(int)
step_to_time = {int(s): t for s, t in zip(ssa_steps_int, ssa_time)}
step_to_ssa  = {int(s): v for s, v in zip(ssa_steps_int, ssa_values)}

# Try a small set of offsets to account for logging that started at 0 or 1
def try_with_offset(off):
    mt, ms, missing = [], [], []
    for ts in range(n_keff):
        key = ts + off
        if key in step_to_time:
            mt.append(step_to_time[key])
            ms.append(step_to_ssa[key])
        else:
            missing.append(ts)
    return np.array(mt), np.array(ms), missing

best = (np.array([]), np.array([]), [])
best_off = 0
for off in (0, -1, 1, -2, 2):
    mt, ms, miss = try_with_offset(off)
    if mt.size > best[0].size:
        best = (mt, ms, miss)
        best_off = off

matched_times, matched_ssa, missing_steps = best

# --- Step 4b: Optional normalization by initial values ---
if normalize_values:
    if matched_ssa.size > 0:
        matched_ssa = matched_ssa / matched_ssa[0]
    if k_xx.size > 0:
        k_xx = k_xx / k_xx[0]
    if k_yy.size > 0:
        k_yy = k_yy / k_yy[0]

# --- Step 4c: Auto-select time units and scale times for plotting ---
scaled_times = matched_times
time_unit_label = 'seconds'
if matched_times.size > 0:
    tmax = float(np.nanmax(matched_times))
    if tmax >= 2 * 24 * 3600:
        scaled_times = matched_times / 86400.0
        time_unit_label = 'days'
    elif tmax >= 2 * 3600:
        scaled_times = matched_times / 3600.0
        time_unit_label = 'hours'
    elif tmax >= 120:
        scaled_times = matched_times / 60.0
        time_unit_label = 'minutes'
    else:
        scaled_times = matched_times
        time_unit_label = 'seconds'

# --- Plot Settings ---
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.linewidth': 1.2,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'figure.dpi': 100
})

def beautify_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, linestyle='--', alpha=0.6)
    return ax

# --- Step 5: Plot k_eff vs Time ---
fig, ax = plt.subplots(figsize=(8, 6))
if matched_times.size == 0:
    print("WARNING: No matched times; skipping k_eff vs time plot.")
else:
    ax.plot(scaled_times, k_xx, label=r'$k_{xx}$')
    ax.plot(scaled_times, k_yy, label=r'$k_{yy}$')
    ax.set_xlabel(f'Time ({time_unit_label})')
    if normalize_values:
        ax.set_ylabel(r'Normalized Effective Thermal Conductivity, $k/k_0$')
        ax.set_title('Normalized Effective Thermal Conductivity vs Time')
    else:
        ax.set_ylabel(r'Effective Thermal Conductivity')
        ax.set_title('Effective Thermal Conductivity vs Time')
    ax.legend()
    beautify_axes(ax)
    fig.tight_layout()
    if save_figures:
        fig.savefig(os.path.join(output_dir, 'k_eff_vs_time.png'), dpi=dpi, bbox_inches='tight', transparent=True)
plt.show()

# --- Step 6: SSA vs k_eff, scatter or curve plot ---
fig, ax = plt.subplots(figsize=(8, 6))
if matched_ssa.size == 0:
    print("WARNING: No matched SSA; skipping k_eff vs SSA plot.")
else:
    if use_scatter_for_ssa:
        scatter_indices = np.arange(0, len(matched_ssa), scatter_step)
        scatter = ax.scatter(matched_ssa[scatter_indices], k_xx[scatter_indices], 
                             c=scaled_times[scatter_indices], cmap='plasma', 
                             label=r'$k_{xx}$', edgecolor='k', s=50, alpha=0.8)
        ax.scatter(matched_ssa[scatter_indices], k_yy[scatter_indices], 
                   c=scaled_times[scatter_indices], cmap='plasma', 
                   label=r'$k_{yy}$', marker='s', edgecolor='k', s=50, alpha=0.8)
        cbar = fig.colorbar(scatter)
        cbar.set_label(f'Time ({time_unit_label})')
    else:
        ax.plot(matched_ssa, k_xx, 'o-', label=r'$k_{xx}$')
        ax.plot(matched_ssa, k_yy, 's-', label=r'$k_{yy}$')

    if normalize_values:
        ax.set_xlabel(r'Normalized Specific Surface Area, $\mathrm{SSA}/\mathrm{SSA}_0$')
        ax.set_ylabel(r'Normalized Effective Thermal Conductivity, $k/k_0$')
        ax.set_title('Normalized Effective Thermal Conductivity vs Normalized SSA')
    else:
        ax.set_xlabel('Specific Surface Area')
        ax.set_ylabel('Effective Thermal Conductivity')
        ax.set_title('Effective Thermal Conductivity vs Specific Surface Area')
    ax.legend()
    beautify_axes(ax)
    fig.tight_layout()
    if save_figures:
        fig.savefig(os.path.join(output_dir, 'k_eff_vs_ssa.png'), dpi=dpi, bbox_inches='tight', transparent=True)
plt.show()

# --- Step 7: SSA vs Time, scatter or curve plot ---
fig, ax = plt.subplots(figsize=(8, 6))
if matched_times.size == 0:
    print("WARNING: No matched times; skipping SSA vs time plot.")
else:
    ax.plot(scaled_times, matched_ssa, '-', color='tab:green')
    ax.set_xlabel(f'Time ({time_unit_label})')
    if normalize_values:
        ax.set_ylabel(r'Normalized Specific Surface Area, $\mathrm{SSA}/\mathrm{SSA}_0$')
        ax.set_title('Normalized SSA vs Time')
    else:
        ax.set_ylabel('Specific Surface Area')
        ax.set_title('SSA vs Time')
    beautify_axes(ax)
    fig.tight_layout()
    if save_figures:
        fig.savefig(os.path.join(output_dir, 'ssa_vs_time.png'), dpi=dpi, bbox_inches='tight', transparent=True)
plt.show()