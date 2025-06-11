import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Parameters ---
parent_dir = '/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/'
folder_pattern = r'NASAv2_30G-2D_T([-\d.]+)_hum[\d.]+_'
save_figures = False  # Skip plotting for now
scatter_step = 10
dpi = 300

# --- Storage for global data ---
collected_data = []  # Will hold dicts with 'T', 'time', 'k_xx'

# --- Plot Settings (still useful if you re-enable plots) ---
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

    # Step 2: Load k_eff
    k_eff_df = pd.read_csv(k_eff_file)
    k_xx = k_eff_df['k_00'].values

    # Step 3: Get time steps from sol_#.dat files
    sol_files = [f for f in os.listdir(folder_path) if re.match(r'sol_\d+\.dat', f)]
    time_steps = sorted([int(re.findall(r'\d+', f)[0]) for f in sol_files])

    # Step 4: Match time and k_eff
    time_step_to_time = {int(ts): t for ts, t in zip(ssa_time_steps, ssa_time)}
    matched_times = []
    matched_kxx = []
    for ts in time_steps[:len(k_xx)]:
        if ts in time_step_to_time:
            matched_times.append(time_step_to_time[ts])
            matched_kxx.append(k_xx[len(matched_times)-1])
        else:
            print(f"Time step {ts} not found in SSA_evo.dat for {subfolder}")

    matched_times = np.array(matched_times)
    matched_kxx = np.array(matched_kxx)

    # Step 5: Store for colormap later
    collected_data.append({
        'temperature': temperature,
        'time': matched_times,
        'k_xx': matched_kxx
    })

print(f"\n✅ Done collecting data from {len(collected_data)} simulations.")

# Optional: show structure of collected_data
for entry in collected_data:
    print(f"T={entry['temperature']} K | timesteps={len(entry['time'])} | k_xx range=({entry['k_xx'].min():.3g}, {entry['k_xx'].max():.3g})")

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# --- Step 1: Sort data by temperature ---
collected_data.sort(key=lambda d: d['temperature'])
temperatures = np.array([d['temperature'] for d in collected_data])

# --- Step 2: Create a common time grid ---
# Use union of all time points, or uniform linspace over global min/max time
all_times = np.concatenate([d['time'] for d in collected_data])
t_min, t_max = all_times.min(), all_times.max()
Nt = 200  # resolution in time
common_time = np.linspace(t_min, t_max, Nt)

# --- Step 3: Interpolate k_xx to common time for each simulation ---
k_map = []  # rows: different temperatures, columns: interpolated k_xx over common_time
for d in collected_data:
    f_interp = interp1d(d['time'], d['k_xx'], bounds_error=False, fill_value=np.nan)
    k_interp = f_interp(common_time)
    k_map.append(k_interp)

k_map = np.array(k_map)  # shape: (n_temperatures, Nt)

# --- Step 4: Plot colormap ---
fig, ax = plt.subplots(figsize=(10, 6))
c = ax.pcolormesh(temperatures, common_time, k_map.T, shading='auto', cmap='plasma')
cb = fig.colorbar(c, ax=ax, label=r'$k_{xx}$ (W/m·K)')

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Time')
ax.set_title('Effective Thermal Conductivity $k_{xx}$ vs Temperature and Time')
beautify_axes(ax)
fig.tight_layout()
fig.savefig('./outputs/homog/kxx_colormap.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()