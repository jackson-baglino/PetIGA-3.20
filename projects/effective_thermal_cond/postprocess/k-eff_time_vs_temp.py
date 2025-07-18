import os
import re
import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# --- Parameters ---
parent_dir = '/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch/Test'
save_figures = False
dpi = 1200

# --- Matplotlib settings ---
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

def determine_time_unit(all_times):
    max_time = max(t[-1] for t in all_times)
    return ('days', 86400) if max_time > 2 * 86400 else ('hours', 3600)

# --- Data Collection ---
collected_data = []
all_times_for_unit = []

for subfolder in sorted(os.listdir(parent_dir)):
    folder_path = os.path.join(parent_dir, subfolder)
    if not os.path.isdir(folder_path):
        continue

    temp_match = re.search(r'Tm?(-?\d+\.?\d*)', subfolder)
    if not temp_match:
        print(f"⚠️ Skipping '{subfolder}': No temperature found.")
        continue
    temperature = float(temp_match.group(1))

    ssa_file = os.path.join(folder_path, 'SSA_evo.dat')
    k_eff_file = os.path.join(folder_path, 'k_eff.csv')
    metadata_file = os.path.join(folder_path, 'metadata.json')

    if not os.path.exists(ssa_file) or not os.path.exists(k_eff_file):
        print(f"⚠️ Skipping '{subfolder}': Missing SSA or k_eff data.")
        continue

    # Read metadata.json
    if os.path.exists(metadata_file):
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        domain_size = metadata.get("domain_size_m", {})
        Lx = domain_size.get("Lx", None)
        Ly = domain_size.get("Ly", None)
        Lz = domain_size.get("Lz", None)
    else:
        print(f"⚠️ No metadata.json found in {folder_path}")
        Lx = Ly = Lz = None

    # Read SSA and k_eff
    ssa_data = np.loadtxt(ssa_file)
    ssa_time = ssa_data[:, 2]
    ssa_steps = ssa_data[:, 3]
    time_lookup = {int(step): t for step, t in zip(ssa_steps, ssa_time)}

    k_df = pd.read_csv(k_eff_file)
    k_xx = k_df['k_00'].values
    steps = k_df['sol_index'].values

    matched_times = [time_lookup[s] for s in steps if s in time_lookup]
    matched_kxx = [k for s, k in zip(steps, k_xx) if s in time_lookup]

    if len(matched_times) == 0:
        print(f"⚠️ Skipping '{subfolder}': No matched steps.")
        continue

    matched_porosity = [
        ssa_data[np.where(ssa_steps == s)[0][0], 1]
        for s in steps if s in time_lookup
    ]

    collected_data.append({
        'temperature': temperature,
        'time_sec': np.array(matched_times),
        'k_xx': np.array(matched_kxx),
        'porosity': np.array(matched_porosity),
        'folder': folder_path,
        'Lx': Lx,
        'Ly': Ly,
        'Lz': Lz
    })

    all_times_for_unit.append(np.array(matched_times))