import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# --- Parameters ---
parent_dir = '/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch'
save_figures = True
dpi = 300

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

# --- Determine time unit globally based on max time ---
def determine_time_unit(all_times):
    max_time = max(t[-1] for t in all_times)
    if max_time > 2 * 86400:  # more than 2 days
        return 'days', 86400
    else:
        return 'hours', 3600

# --- Collect simulation data ---
collected_data = []
all_times_for_unit = []

for subfolder in sorted(os.listdir(parent_dir)):
    folder_path = os.path.join(parent_dir, subfolder)
    if not os.path.isdir(folder_path):
        continue

    # Extract temperature from folder name (e.g., Tm80.0 → -80.0)
    temp_match = re.search(r'Tm?(-?\d+\.?\d*)', subfolder)
    if not temp_match:
        print(f"⚠️ Skipping '{subfolder}': No temperature found.")
        continue
    temperature = float(temp_match.group(1))

    ssa_file = os.path.join(folder_path, 'SSA_evo.dat')
    k_eff_file = os.path.join(folder_path, 'k_eff.csv')

    if not os.path.exists(ssa_file) or not os.path.exists(k_eff_file):
        print(f"⚠️ Skipping '{subfolder}': Missing SSA or k_eff data.")
        continue

    # Load SSA data
    ssa_data = np.loadtxt(ssa_file)
    ssa_time = ssa_data[:, 2]  # in seconds
    ssa_steps = ssa_data[:, 3]
    time_lookup = {int(step): t for step, t in zip(ssa_steps, ssa_time)}

    # Load k_eff data
    k_df = pd.read_csv(k_eff_file)
    k_xx = k_df['k_00'].values
    steps = k_df['sol_index'].values

    matched_times = [time_lookup[s] for s in steps if s in time_lookup]
    matched_kxx = [k for s, k in zip(steps, k_xx) if s in time_lookup]

    if len(matched_times) == 0:
        print(f"⚠️ Skipping '{subfolder}': No matched steps.")
        continue

    matched_times = np.array(matched_times)
    matched_kxx = np.array(matched_kxx)
    matched_porosity = [ssa_data[np.where(ssa_steps == s)[0][0], 1] for s in steps if s in time_lookup]

    collected_data.append({
        'temperature': temperature,
        'time_sec': matched_times,
        'k_xx': matched_kxx,
        'porosity': np.array(matched_porosity),
        'folder': folder_path
    })

    all_times_for_unit.append(matched_times)

# --- Decide time unit ---
time_unit, time_scale = determine_time_unit(all_times_for_unit)
print(f"\n📏 Using time unit: {time_unit}")

# --- Plot individual curves ---
for entry in collected_data:
  times = entry['time_sec'] / time_scale
  kxx = entry['k_xx']
  porosity = entry['porosity']
  temp = entry['temperature']

  # --- Plot kxx ---
  fig1, ax1 = plt.subplots(figsize=(6, 4))
  ax1.plot(times, kxx, marker='o', lw=1.5)
  ax1.set_xlabel(f'Time ({time_unit})')
  ax1.set_ylabel(r'$k_{xx}$ (W/m·K)')
  ax1.set_title(f'T = {temp:.1f} °C — $k_{{xx}}$')
  beautify_axes(ax1)
  fig1.tight_layout()
  if save_figures:
    fig1.savefig(os.path.join(entry['folder'], 'kxx_vs_time.png'),
                  dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig1)
  else:
    plt.show()

  # --- Normalize porosity ---
  porosity = entry['porosity']
  norm_porosity = porosity / np.max(porosity)

  # --- Plot normalized porosity ---
  fig2, ax2 = plt.subplots(figsize=(6, 4))
  ax2.plot(times, norm_porosity, marker='s', color='teal', lw=1.5)
  ax2.set_xlabel(f'Time ({time_unit})')
  ax2.set_ylabel('Normalized Porosity')
  ax2.set_title(f'T = {temp:.1f} °C — Normalized Porosity')
  ax2.set_ylim(0, 1.05)
  beautify_axes(ax2)
  fig2.tight_layout()
  if save_figures:
    fig2.savefig(os.path.join(entry['folder'], 'porosity_vs_time_normalized.png'),
                  dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig2)
  else:
    plt.show()

# --- Interpolate onto common time grid ---
collected_data.sort(key=lambda d: d['temperature'])
temps = np.array([d['temperature'] for d in collected_data])
common_t = np.linspace(
    max(d['time_sec'][0] for d in collected_data) / time_scale,
    min(d['time_sec'][-1] for d in collected_data) / time_scale,
    200
)

k_map, mask = [], []
for d in collected_data:
    times = d['time_sec'] / time_scale
    f = interp1d(times, d['k_xx'], bounds_error=False, fill_value=np.nan)
    values = f(common_t)
    k_map.append(values)
    mask.append(np.isfinite(values))

k_map = np.ma.array(k_map, mask=~np.array(mask))

# --- Heatmap ---
# Sort temperatures and k_map
sorted_indices = np.argsort(temps)
sorted_kmap = k_map[sorted_indices]
sorted_temps = temps[sorted_indices]

# Calculate uniform temperature spacing (dT)
if len(sorted_temps) >= 2:
    temp_deltas = np.diff(sorted_temps)
    print(f"\n📏 Temperature deltas: {temp_deltas}")
    dT = np.mean(temp_deltas)
else:
    raise ValueError("At least two temperature entries are required.")

# Compute padded temperature bin edges (each temp centered in a bin)
T_min = sorted_temps[0] - dT / 2
T_max = sorted_temps[-1] + dT / 2
temp_edges = np.linspace(T_min, T_max, len(sorted_temps) + 1)

# Debug prints
print("\n📊 Heatmap Debug Info:")
print(f"  Temperatures: {sorted_temps}")
print(f"  dT: {dT}")
print(f"  Temp edges: {temp_edges}")
print(f"  k_map shape: {sorted_kmap.shape} (should be Nt x Ntemps = {common_t.size} x {len(sorted_temps)})")

# Compute time bin edges
dt = common_t[1] - common_t[0]
time_edges = np.concatenate(([common_t[0] - dt / 2], common_t + dt / 2))

# Plot the heatmap
fig, ax = plt.subplots(figsize=(10, 6))
mesh = ax.pcolormesh(
    temp_edges,
    time_edges,
    sorted_kmap.T,
    shading='auto',
    cmap='plasma'
)

cb = fig.colorbar(mesh, ax=ax, label=r'$k_{xx}$ (W/m·K)')
ax.set_xlabel('Temperature (°C)')
ax.set_ylabel(f'Time ({time_unit})')
ax.set_title(r'$k_{xx}$ vs Temperature and Time')
beautify_axes(ax)
fig.tight_layout()
fig.savefig(os.path.join(parent_dir, 'kxx_colormap.png'),
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()

# --- Combined line plot ---
fig, ax = plt.subplots(figsize=(10, 6))
for d in collected_data:
    t = d['time_sec'] / time_scale
    ax.plot(t, d['k_xx'], label=f"{d['temperature']} °C")

ax.set_xlabel(f'Time ({time_unit})')
ax.set_ylabel(r'$k_{xx}$ (W/m·K)')
ax.set_title(r'$k_{xx}$ for All Simulations')
beautify_axes(ax)
ax.legend(title='Temperature', loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()

fig.savefig(os.path.join(parent_dir, 'kxx_vs_time_all.png'),
            dpi=dpi, bbox_inches='tight', transparent=True)

if not save_figures:
    plt.show()
else:
    plt.close(fig)