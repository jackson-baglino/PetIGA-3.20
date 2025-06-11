import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Parameters ---
folder_path = '/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/NASAv2_30G-2D_T-80.0_hum0.70_2025-06-09__16.04.46'
ssa_file = os.path.join(folder_path, 'SSA_evo.dat')
k_eff_file = os.path.join(folder_path, 'k_eff.csv')  # CSV file with header

# --- Plotting Options ---
use_scatter_for_ssa = True    # Set True for scatter plots, False for line plots
scatter_step = 10             # Plot every Nth point if using scatter
save_figures = True           # Save figures to disk
dpi = 300                     # DPI for saved figures

# --- Output Directory Setup ---
parent_folder_name = os.path.basename(os.path.normpath(folder_path))
output_dir = f'./outputs/homog/plots/{parent_folder_name}'
os.makedirs(output_dir, exist_ok=True)

# --- Step 1: Find all sol_#.dat files and extract time steps ---
sol_files = [f for f in os.listdir(folder_path) if re.match(r'sol_\d+\.dat', f)]
time_steps = sorted([int(re.findall(r'\d+', f)[0]) for f in sol_files])

# --- Step 2: Load SSA_evo.dat ---
ssa_data = np.loadtxt(ssa_file)
ssa_time = ssa_data[:, 2]        # time column
ssa_time_steps = ssa_data[:, 3]  # time step column
ssa_values = ssa_data[:, 0]      # SSA column (1st column)

# --- Step 3: Load k_eff.csv ---
k_eff_df = pd.read_csv(k_eff_file)
k_xx = k_eff_df['k_00'].values
k_yy = k_eff_df['k_11'].values

# --- Step 4: Match time and SSA for each time step in sol_#.dat ---
time_step_to_time = {int(ts): t for ts, t in zip(ssa_time_steps, ssa_time)}
time_step_to_ssa = {int(ts): ssa for ts, ssa in zip(ssa_time_steps, ssa_values)}

matched_times = []
matched_ssa = []
for ts in time_steps[:len(k_xx)]:
    if ts in time_step_to_time:
        matched_times.append(time_step_to_time[ts])
        matched_ssa.append(time_step_to_ssa[ts])
    else:
        print(f"Warning: time step {ts} not found in SSA_evo.dat")

matched_times = np.array(matched_times)
matched_ssa = np.array(matched_ssa)

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
ax.plot(matched_times, k_xx, label=r'$k_{xx}$')
ax.plot(matched_times, k_yy, label=r'$k_{yy}$')
ax.set_xlabel('Time')
ax.set_ylabel('Effective Thermal Conductivity')
ax.set_title('Effective Thermal Conductivity vs Time')
ax.legend()
beautify_axes(ax)
fig.tight_layout()
if save_figures:
    fig.savefig(os.path.join(output_dir, 'k_eff_vs_time.png'), dpi=dpi, bbox_inches='tight', transparent=True)
plt.show()

# --- Step 6: SSA vs k_eff, scatter or curve plot ---
fig, ax = plt.subplots(figsize=(8, 6))

if use_scatter_for_ssa:
    scatter_indices = np.arange(0, len(matched_ssa), scatter_step)
    scatter = ax.scatter(matched_ssa[scatter_indices], k_xx[scatter_indices], 
                         c=matched_times[scatter_indices], cmap='plasma', 
                         label=r'$k_{xx}$', edgecolor='k', s=50, alpha=0.8)
    ax.scatter(matched_ssa[scatter_indices], k_yy[scatter_indices], 
               c=matched_times[scatter_indices], cmap='plasma', 
               label=r'$k_{yy}$', marker='s', edgecolor='k', s=50, alpha=0.8)
    cbar = fig.colorbar(scatter)
    cbar.set_label('Time')
else:
    ax.plot(matched_ssa, k_xx, 'o-', label=r'$k_{xx}$')
    ax.plot(matched_ssa, k_yy, 's-', label=r'$k_{yy}$')

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
ax.plot(matched_times, matched_ssa, '-', color='tab:green')
ax.set_xlabel('Time')
ax.set_ylabel('Specific Surface Area')
ax.set_title('Specific Surface Area vs Time')
beautify_axes(ax)
fig.tight_layout()
if save_figures:
    fig.savefig(os.path.join(output_dir, 'ssa_vs_time.png'), dpi=dpi, bbox_inches='tight', transparent=True)
plt.show()