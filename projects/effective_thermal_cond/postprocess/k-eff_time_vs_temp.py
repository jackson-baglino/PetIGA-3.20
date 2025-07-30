import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# --- Parameters ---
parent_dir = '/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/'
folder_pattern = r'NASAv2_2G-Molaro-2D_T([-\d.]+)_hum[\d.]+_'
save_figures = True  # If True, save plots to disk; if False, show plots interactively
scatter_step = 10
dpi = 300

# --- Storage for global data ---
collected_data = []  # Will hold dicts with 'T', 'time', 'k_xx'

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

# --- Loop through matching folders ---
for subfolder in sorted(os.listdir(parent_dir)):
    match = re.match(folder_pattern, subfolder)
    if not match:
        continue  # skip non-matching folders

    temperature = float(match.group(1))
    folder_path = os.path.join(parent_dir, subfolder)
    print(f'Processing {subfolder} at T = {temperature} K')

    ssa_file = os.path.join(folder_path, 'SSA_evo.dat')
    k_eff_file = os.path.join(folder_path, 'k_eff.csv')

    if not os.path.exists(ssa_file) or not os.path.exists(k_eff_file):
        print(f"Missing SSA or k_eff file in {subfolder}, skipping.")
        continue

    # Step 1: Load SSA
    ssa_data = np.loadtxt(ssa_file)
    ssa_time = ssa_data[:, 2]
    ssa_time_steps = ssa_data[:, 3]
    ssa_values = ssa_data[:, 0]

    # Step 2: Load k_eff (now includes a "step" column!)
    k_eff_df = pd.read_csv(k_eff_file)
    k_xx = k_eff_df['k_00'].values
    steps = k_eff_df['sol_index'].values  # This column should match SSA steps

    # Step 3: Build a lookup from SSA time step → time
    time_step_to_time = {int(ts): t for ts, t in zip(ssa_time_steps, ssa_time)}

    # Step 4: Match times using step column from k_eff.csv
    matched_times = []
    matched_kxx = []

    for i, step in enumerate(steps):
        if step in time_step_to_time:
            matched_times.append(time_step_to_time[step])
            matched_kxx.append(k_xx[i])
        else:
            print(f"Step {step} not found in SSA data for folder {subfolder}")

    matched_times = np.array(matched_times)
    matched_kxx = np.array(matched_kxx)

    # Step 5: Store for colormap later
    collected_data.append({
        'temperature': temperature,
        'time': matched_times,
        'k_xx': matched_kxx
    })

    # Step 6: Plot and save/show individual plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(matched_times, matched_kxx, marker='o', lw=1.5)
    ax.set_xlabel('Time')
    ax.set_ylabel(r'$k_{xx}$ (W/m·K)')
    ax.set_title(f'T = {temperature} K')
    beautify_axes(ax)
    fig.tight_layout()

    if save_figures:
        fig_path = os.path.join(folder_path, 'kxx_vs_time.png')
        fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)
    else:
        plt.show()

# --- Summary Output ---
print(f"\n✅ Done collecting data from {len(collected_data)} simulations.")
for entry in collected_data:
    print(f"T={entry['temperature']} K | timesteps={len(entry['time'])} | k_xx range=({entry['k_xx'].min():.3g}, {entry['k_xx'].max():.3g})")

# --- Step 7: Build masked heatmap with per-simulation extents ---

# Sort by temperature
collected_data.sort(key=lambda d: d['temperature'])
temperatures = np.array([d['temperature'] for d in collected_data])

# Create common time grid (but we won't force all data to fill it)
final_times = np.array([d['time'][-1] for d in collected_data])
start_times = np.array([d['time'][0] for d in collected_data])
t_min = max(start_times)
t_max = min(final_times)
Nt = 200
common_time = np.linspace(t_min, t_max, Nt)

# Interpolate to common time with masking
k_map = []
mask = []

for d in collected_data:
    f_interp = interp1d(d['time'], d['k_xx'], bounds_error=False, fill_value=np.nan)
    k_interp = f_interp(common_time)
    k_map.append(k_interp)

    valid_mask = np.isfinite(f_interp(common_time))
    mask.append(valid_mask)

k_map = np.ma.array(k_map, mask=~np.array(mask))  # Mask invalid entries
k_map = np.ma.masked_invalid(k_map)

# --- Plot masked heatmap ---
fig, ax = plt.subplots(figsize=(10, 6))
c = ax.pcolormesh(temperatures, common_time, k_map.T, shading='auto', cmap='plasma')
cb = fig.colorbar(c, ax=ax, label=r'$k_{xx}$ (W/m·K)')

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Time')
ax.set_title('Effective Thermal Conductivity $k_{xx}$ vs Temperature and Time')
beautify_axes(ax)
fig.tight_layout()
fig.savefig('./outputs/homog/kxx_colormap_masked.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# --- Plot all k_xx vs time curves together ---
fig, ax = plt.subplots(figsize=(10, 6))

for entry in collected_data:
    ax.plot(entry['time'], entry['k_xx'], label=f"T = {entry['temperature']} K")

ax.set_xlabel('Time')
ax.set_ylabel(r'$k_{xx}$ (W/m·K)')
ax.set_title(r'Effective Thermal Conductivity $k_{xx}$ for All Simulations')
beautify_axes(ax)
ax.legend(title='Temperature', loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()

# Save to the parent_dir (one level above the subfolders)
multi_plot_path = os.path.join(parent_dir, 'kxx_vs_time_all.png')
fig.savefig(multi_plot_path, dpi=dpi, bbox_inches='tight', transparent=True)

if not save_figures:
    plt.show()
else:
    plt.close(fig)