import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# --- Parameters ---
parent_dir = os.path.expanduser('/Users/jacksonbaglino/SimulationResults/effective_thermal_cond/scratch')
folder_pattern = r'^grains_2D_Tm([+-]?\d+(?:\.\d+)?)_hum\d+_tf[^_]+__\d{4}-\d{2}-\d{2}__\d{2}\.\d{2}(?:\.\d{2})?$'
save_figures = True   # If True, save plots to disk; if False, show plots interactively
dpi = 600             # DPI for saved figures
normalize_values = True   # Plot relative quantities: normalize k_xx, k_yy, and SSA by their initial values per simulation

# --- Storage for global data ---
collected_data = []  # Will hold dicts with 'T', 'time', 'k_yy'

# --- Plot Settings ---
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 20,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'legend.fontsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'axes.linewidth': 1.2,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'figure.dpi': 100
})

# --- High-contrast palette and helper functions ---
import matplotlib as mpl
# Build a high-contrast palette that scales beyond 10 series
PALETTE = list(mpl.colors.TABLEAU_COLORS.values()) \
          + [mpl.cm.Dark2(i) for i in range(mpl.cm.Dark2.N)] \
          + [mpl.cm.Set1(i) for i in range(mpl.cm.Set1.N)]

def color_for(i: int):
    return PALETTE[i % len(PALETTE)]

def sparse_markevery(x, target=12):
    """Return a stride so only ~`target` markers are drawn across series length."""
    n = len(x)
    if n <= 0:
        return 1
    if n <= target:
        return 1
    return max(1, n // target)

def beautify_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Remove dotted grid lines
    ax.grid(False)
    return ax

# --- Legend/SSA helpers ---
from typing import Sequence, Tuple

def place_legend_outside(ax, handles: Sequence, labels: Sequence[str], loc='center left', anchor=(0.98, 0.5), right_pad=0.80):
    """Place a *figure-level* legend outside and reserve right margin so it never overlaps.
    Returns the recommended rect for tight_layout.
    """
    fig = ax.get_figure()
    if handles and labels:
        fig.legend(handles, labels, loc=loc, bbox_to_anchor=anchor)
        # Reserve right margin for the legend
        fig.subplots_adjust(right=right_pad)
    # Recommend a tight_layout rect matching the reserved margin
    return [0.0, 0.0, right_pad, 1.0]


def add_secondary_ssa_axis(ax, time: np.ndarray, ssa: np.ndarray) -> Tuple[bool, object]:
    """Attach a top x-axis (SSA) that's functionally linked to time via interpolation.
    Returns (ok, axis). If SSA is not monotonic, falls back to labeling a few ticks.
    """
    # Clean
    valid = np.isfinite(time) & np.isfinite(ssa)
    t = np.asarray(time)[valid]
    s = np.asarray(ssa)[valid]
    if t.size < 2 or s.size < 2:
        return False, None

    # Ensure strictly increasing time for interpolation
    order = np.argsort(t)
    t = t[order]
    s = s[order]

    # Check monotonic SSA to define inverse
    inc = np.all(np.diff(s) > 0)
    dec = np.all(np.diff(s) < 0)
    if inc or dec:
        f_ts = interp1d(t, s, kind='linear', bounds_error=False, fill_value='extrapolate')
        f_st = interp1d(s, t, kind='linear', bounds_error=False, fill_value='extrapolate')
        secax = ax.secondary_xaxis('top', functions=(lambda tt: f_ts(tt), lambda ss: f_st(ss)))
        secax.set_xlabel(r'SSA/SSA$_0$' if normalize_values else 'Specific Surface Area')
        return True, secax
    else:
        # Fallback: create a twin axis and label a few corresponding ticks without strict mapping
        secax = ax.twiny()
        # Choose ~5 ticks evenly in time and label with SSA values at those times
        ticks_t = np.linspace(t.min(), t.max(), 5)
        ticks_s = np.interp(ticks_t, t, s)
        secax.set_xlim(ax.get_xlim())
        secax.set_xticks(ticks_t)
        secax.set_xticklabels([f"{val:.2g}" for val in ticks_s])
        secax.set_xlabel(r'SSA/SSA$_0$' if normalize_values else 'Specific Surface Area')
        return False, secax

# --- Helpers ---
def autoscale_time(arr):
    unit = 'seconds'
    scaled = arr
    if arr.size > 0:
        tmax = float(np.nanmax(arr))
        if tmax >= 2 * 24 * 3600:
            scaled, unit = arr / 86400.0, 'days'
        elif tmax >= 2 * 3600:
            scaled, unit = arr / 3600.0, 'hours'
        elif tmax >= 120:
            scaled, unit = arr / 60.0, 'minutes'
    return scaled, unit

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
    k_yy = k_eff_df['k_11'].values
    steps = k_eff_df['sol_index'].astype(int).values  # ensure int keys to match SSA steps

    # Step 3: Build a lookup from SSA time step → time
    time_step_to_time = {int(ts): t for ts, t in zip(ssa_time_steps, ssa_time)}
    time_step_to_ssa  = {int(ts): v for ts, v in zip(ssa_time_steps, ssa_values)}  # NEW: step → SSA

    # Step 4: Match times using step column from k_eff.csv
    matched_times = []
    matched_kxx = []
    matched_kyy = []
    matched_ssa = []  # SSA aligned to k

    for i, step in enumerate(steps):
        if step in time_step_to_time:
            matched_times.append(time_step_to_time[step])
            matched_kxx.append(k_xx[i])
            matched_kyy.append(k_yy[i])
            if step in time_step_to_ssa:
                matched_ssa.append(time_step_to_ssa[step])
            else:
                matched_ssa.append(np.nan)
        else:
            print(f"Step {step} not found in SSA data for folder {subfolder}")

    matched_times = np.array(matched_times)
    matched_kxx = np.array(matched_kxx)
    matched_kyy = np.array(matched_kyy)
    matched_ssa = np.array(matched_ssa)

    # Optional normalization by the first value in these series
    if normalize_values:
        if matched_kxx.size > 0:
            matched_kxx = matched_kxx / matched_kxx[0]
        if matched_kyy.size > 0:
            matched_kyy = matched_kyy / matched_kyy[0]
        if matched_ssa.size > 0:
            matched_ssa = matched_ssa / matched_ssa[0]

    # Auto-scale time units for plotting
    scaled_times, time_unit_label = autoscale_time(matched_times)

    # Store for later summary/plots
    collected_data.append({
        'temperature': temperature,
        'time': matched_times,
        'time_scaled': scaled_times,
        'time_unit': time_unit_label,
        'k_xx': matched_kxx,
        'k_yy': matched_kyy,
        'ssa': matched_ssa,
    })

    # Plot and save/show individual plot
    fig, ax = plt.subplots(figsize=(11.28, 4.2))
    ax.plot(scaled_times, matched_kyy, marker='o', lw=1.5)
    ax.set_xlabel(f'Time ({time_unit_label})')
    if normalize_values:
        ax.set_ylabel(r'$k_{yy}/k_{yy,0}$')
        ax.set_title(rf'T = {temperature} K')
    else:
        ax.set_ylabel(r'$k_{yy}$ (W/m·K)')
        ax.set_title(rf'T = {temperature} K')
    beautify_axes(ax)
    fig.tight_layout()

    if save_figures:
        fig_path = os.path.join(folder_path, 'kyy_vs_time.png')
        fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)
    else:
        plt.show()

    # Plot and save/show k_xx vs Time for this folder
    fig, ax = plt.subplots(figsize=(10.0, 3.0))
    ax.plot(scaled_times, matched_kxx, marker='o', lw=1.5)
    ax.set_xlabel(f'Time ({time_unit_label})')
    if normalize_values:
        ax.set_ylabel(r'$k_{xx}/k_{xx,0}$')
        ax.set_title(rf'T = {temperature} K')
    else:
        ax.set_ylabel(r'$k_{xx}$ (W/m·K)')
        ax.set_title(rf'T = {temperature} K')
    beautify_axes(ax)
    fig.tight_layout()
    if save_figures:
        fig_path = os.path.join(folder_path, 'kxx_vs_time.png')
        fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)
    else:
        plt.show()

    # Combined k_xx and k_yy vs Time (dual y-axes)
    fig, ax1 = plt.subplots(figsize=(11.28, 4.0))
    ln1 = ax1.plot(scaled_times, matched_kxx, marker='o', lw=1.5, label=r'$k_{xx}$')
    ax1.set_xlabel(f'Time ({time_unit_label})')
    if normalize_values:
        ax1.set_ylabel(r'$k_{xx}/k_{xx,0}$')
        ax1.set_title(rf'T = {temperature} K')
    else:
        ax1.set_ylabel(r'$k_{xx}$ (W/m·K)')
        ax1.set_title(rf'T = {temperature} K')
    beautify_axes(ax1)

    ax2 = ax1.twinx()
    ln2 = ax2.plot(scaled_times, matched_kyy, marker='s', lw=1.5, label=r'$k_{yy}$')
    if normalize_values:
        ax2.set_ylabel(r'$k_{yy}/k_{yy,0}$')
    else:
        ax2.set_ylabel(r'$k_{yy}$ (W/m·K)')

    # shared legend combining both axes
    lines = ln1 + ln2
    labels = [l.get_label() for l in lines]
    rect = place_legend_outside(ax1, lines, labels)

    fig.tight_layout(rect=rect)
    if save_figures:
        fig_path = os.path.join(folder_path, 'kxx_kyy_vs_time_dualy.png')
        fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)
    else:
        plt.show()

    # Plot and save/show k_eff vs SSA for this folder (aligned pairs)
    fig, ax = plt.subplots(figsize=(11.28, 3.8))
    valid = np.isfinite(matched_ssa) & np.isfinite(matched_kyy)
    if np.count_nonzero(valid) == 0:
        print(f"No valid SSA↔k pairs to plot in {subfolder}.")
    else:
        ax.plot(matched_ssa[valid], matched_kyy[valid], 'o-', lw=1.5)
        if normalize_values:
            ax.set_xlabel(r'SSA/SSA$_0$')
            ax.set_ylabel(r'$k_{yy}/k_{yy,0}$')
            ax.set_title(rf'T = {temperature} K')
        else:
            ax.set_xlabel('Specific Surface Area')
            ax.set_ylabel(r'$k_{yy}$ (W/m·K)')
            ax.set_title(rf'T = {temperature} K')
        beautify_axes(ax)
        fig.tight_layout()
        if save_figures:
            fig_path = os.path.join(folder_path, 'kyy_vs_ssa.png')
            fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
            plt.close(fig)
        else:
            plt.show()

    # Plot and save/show k_xx vs SSA for this folder (aligned pairs)
    fig, ax = plt.subplots(figsize=(11.28, 3.8))
    valid = np.isfinite(matched_ssa) & np.isfinite(matched_kxx)
    if np.count_nonzero(valid) == 0:
        print(f"No valid SSA↔k_xx pairs to plot in {subfolder}.")
    else:
        ax.plot(matched_ssa[valid], matched_kxx[valid], 'o-', lw=1.5)
        if normalize_values:
            ax.set_xlabel(r'SSA/SSA$_0$')
            ax.set_ylabel(r'$k_{xx}/k_{xx,0}$')
            ax.set_title(rf'T = {temperature} K')
        else:
            ax.set_xlabel('Specific Surface Area')
            ax.set_ylabel(r'$k_{xx}$ (W/m·K)')
            ax.set_title(rf'T = {temperature} K')
        beautify_axes(ax)
        fig.tight_layout()
        if save_figures:
            fig_path = os.path.join(folder_path, 'kxx_vs_ssa.png')
            fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
            plt.close(fig)
        else:
            plt.show()

    # Dual X-Axis combined plot for k_yy (Time bottom, SSA top)
    fig, ax = plt.subplots(figsize=(11.28, 4.0))
    me = sparse_markevery(scaled_times)
    line, = ax.plot(scaled_times, matched_kyy, lw=2.0, marker='o', markevery=me, ms=4, label=r'$k_{yy}$')
    ax.set_xlabel(f'Time ({time_unit_label})')
    ax.set_ylabel(r'$k_{yy}/k_{yy,0}$' if normalize_values else r'$k_{yy}$ (W/m·K)')
    ax.set_title(rf'T = {temperature} K')
    beautify_axes(ax)
    ok, secax = add_secondary_ssa_axis(ax, matched_times, matched_ssa)
    if not ok:
        print(f"[warn] SSA not strictly monotonic in {subfolder}; top axis uses approximate tick labels.")
    fig.tight_layout()
    if save_figures:
        fig_path = os.path.join(folder_path, 'kyy_time_ssa_dualx.png')
        fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)
    else:
        plt.show()

    # Dual X-Axis combined plot for k_xx (Time bottom, SSA top)
    fig, ax = plt.subplots(figsize=(11.28, 4.0))
    me = sparse_markevery(scaled_times)
    line, = ax.plot(scaled_times, matched_kxx, lw=2.0, marker='o', markevery=me, ms=4, label=r'$k_{xx}$')
    ax.set_xlabel(f'Time ({time_unit_label})')
    ax.set_ylabel(r'$k_{xx}/k_{xx,0}$' if normalize_values else r'$k_{xx}$ (W/m·K)')
    ax.set_title(rf'T = {temperature} K')
    beautify_axes(ax)
    ok, secax = add_secondary_ssa_axis(ax, matched_times, matched_ssa)
    if not ok:
        print(f"[warn] SSA not strictly monotonic in {subfolder}; top axis uses approximate tick labels.")
    fig.tight_layout()
    if save_figures:
        fig_path = os.path.join(folder_path, 'kxx_time_ssa_dualx.png')
        fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)
    else:
        plt.show()

    # Combined k_xx and k_yy vs SSA (dual y-axes)
    fig, ax1 = plt.subplots(figsize=(11.28, 4.0))

    valid_xx = np.isfinite(matched_ssa) & np.isfinite(matched_kxx)
    valid_yy = np.isfinite(matched_ssa) & np.isfinite(matched_kyy)

    any_series = False
    if np.count_nonzero(valid_xx) > 0:
        ln1 = ax1.plot(matched_ssa[valid_xx], matched_kxx[valid_xx], 'o-', lw=1.5, label=r'$k_{xx}$')
        any_series = True
    else:
        ln1 = []

    ax1.set_xlabel(r'SSA/SSA$_0$' if normalize_values else 'Specific Surface Area')
    if normalize_values:
        ax1.set_ylabel(r'$k_{xx}/k_{xx,0}$')
        ax1.set_title(rf'T = {temperature} K')
    else:
        ax1.set_ylabel(r'$k_{xx}$ (W/m·K)')
        ax1.set_title(rf'T = {temperature} K')
    beautify_axes(ax1)

    ax2 = ax1.twinx()
    if np.count_nonzero(valid_yy) > 0:
        ln2 = ax2.plot(matched_ssa[valid_yy], matched_kyy[valid_yy], 's-', lw=1.5, label=r'$k_{yy}$')
        any_series = True
    else:
        ln2 = []

    if normalize_values:
        ax2.set_ylabel(r'$k_{yy}/k_{yy,0}$')
    else:
        ax2.set_ylabel(r'$k_{yy}$ (W/m·K)')

    if any_series:
        lines = (ln1 if ln1 else []) + (ln2 if ln2 else [])
        labels = [l.get_label() for l in lines]
        if lines:
            place_rect = place_legend_outside(ax1, lines, labels)
        fig.tight_layout(rect=place_rect)
        if save_figures:
            fig_path = os.path.join(folder_path, 'kxx_kyy_vs_ssa_dualy.png')
            fig.savefig(fig_path, dpi=dpi, bbox_inches='tight', transparent=True)
            plt.close(fig)
        else:
            plt.show()
    else:
        plt.close(fig)
        print(f"No valid SSA↔k_xx/k_yy pairs to plot in {subfolder} (combined).")

# Early exit if nothing was collected (no matching folders or missing files)
if len(collected_data) == 0:
    print("\n⚠️  No simulations collected.\n  • Checked parent_dir:", parent_dir,
          "\n  • Pattern:", folder_pattern,
          "\n  • Ensure each folder contains both SSA_evo.dat and k_eff.csv.")
    import sys
    sys.exit(0)

# --- Summary Output ---
print(f"\n✅ Done collecting data from {len(collected_data)} simulations.")
for entry in collected_data:
    label = 'k/k0' if normalize_values else 'k (W/m·K)'
    print(f"T={entry['temperature']} K | timesteps={len(entry['time'])} | range({label})=({np.nanmin(entry['k_yy']):.3g}, {np.nanmax(entry['k_yy']):.3g})")

# --- Step 7: Build masked heatmaps (k_xx and k_yy) with per-simulation extents ---

# Sort by temperature for consistent vertical ordering
collected_data.sort(key=lambda d: d['temperature'])
temperatures = np.array([d['temperature'] for d in collected_data])

# Build a common time grid only over the overlap of all runs (so we don't extrapolate wildly)
final_times = np.array([d['time'][-1] for d in collected_data])
start_times = np.array([d['time'][0] for d in collected_data])
t_min = max(start_times)
t_max = min(final_times)
if not np.isfinite(t_min) or not np.isfinite(t_max) or t_max <= t_min:
    print("\n⚠️  Not enough overlap in time across simulations to build heatmaps.")
    Nt = 0
else:
    Nt = 200

if Nt > 0:
    common_time = np.linspace(t_min, t_max, Nt)
    common_scaled_time, common_time_unit = autoscale_time(common_time)

    # Interpolate k_yy onto common_time with masking
    kyy_rows = []
    kyy_mask_rows = []
    for d in collected_data:
        f = interp1d(d['time'], d['k_yy'], bounds_error=False, fill_value=np.nan)
        vals = f(common_time)
        kyy_rows.append(vals)
        kyy_mask_rows.append(np.isfinite(vals))
    kyy_map = np.ma.array(kyy_rows, mask=~np.array(kyy_mask_rows))
    kyy_map = np.ma.masked_invalid(kyy_map)

    # Interpolate k_xx onto common_time with masking (if present)
    kxx_rows = []
    kxx_mask_rows = []
    has_kxx = any('k_xx' in d and d['k_xx'] is not None for d in collected_data)
    if has_kxx:
        for d in collected_data:
            if 'k_xx' in d and d['k_xx'] is not None:
                f = interp1d(d['time'], d['k_xx'], bounds_error=False, fill_value=np.nan)
                vals = f(common_time)
            else:
                vals = np.full_like(common_time, np.nan, dtype=float)
            kxx_rows.append(vals)
            kxx_mask_rows.append(np.isfinite(vals))
        kxx_map = np.ma.array(kxx_rows, mask=~np.array(kxx_mask_rows))
        kxx_map = np.ma.masked_invalid(kxx_map)

    # ---- Plot heatmap for k_yy ----
    fig, ax = plt.subplots(figsize=(11.28, 3.8))
    c = ax.pcolormesh(temperatures, common_scaled_time, kyy_map.T, shading='auto', cmap='plasma')
    cb_label = r'$k_{yy}/k_{yy,0}$' if normalize_values else r'$k_{yy}$ (W/m·K)'
    cb = fig.colorbar(c, ax=ax, label=cb_label)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(f'Time ({common_time_unit})')
    ax.set_title(r'$k_{yy}/k_{yy,0}$ vs Temperature and Time' if normalize_values else r'$k_{yy}$ vs Temperature and Time')
    beautify_axes(ax)
    fig.tight_layout()
    out_path = os.path.join(parent_dir, 'kyy_colormap_masked.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)

    # ---- Plot heatmap for k_xx (if available) ----
    if has_kxx:
        fig, ax = plt.subplots(figsize=(11.28, 3.8))
        c = ax.pcolormesh(temperatures, common_scaled_time, kxx_map.T, shading='auto', cmap='plasma')
        cb_label = r'$k_{xx}/k_{xx,0}$' if normalize_values else r'$k_{xx}$ (W/m·K)'
        cb = fig.colorbar(c, ax=ax, label=cb_label)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel(f'Time ({common_time_unit})')
        ax.set_title(r'$k_{xx}/k_{xx,0}$ vs Temperature and Time' if normalize_values else r'$k_{xx}$ vs Temperature and Time')
        beautify_axes(ax)
        fig.tight_layout()
        out_path = os.path.join(parent_dir, 'kxx_colormap_masked.png')
        fig.savefig(out_path, dpi=300, bbox_inches='tight', transparent=True)
        plt.close(fig)
else:
    print("\n⚠️  Skipping heatmaps due to insufficient common time range across runs.")

#
# --- Step 8: Two-panel figure for ALL temperatures: (Left) vs Time, (Right) vs SSA ---
# Plots include BOTH k_xx (solid) and k_yy (dashed) for every temperature, single legend.

fig, (axL, axR) = plt.subplots(1, 2, figsize=(14.5, 4.6), sharey=False, gridspec_kw={
    'wspace': 0.18  # tighter space between subplots
})

# Left: vs Time (auto-scaled units)
units = [d['time_unit'] for d in collected_data]
unit = max(set(units), key=units.count) if units else 'seconds'
for idx, d in enumerate(collected_data):
    col = color_for(idx)
    t = d['time_scaled']
    # k_xx (solid) if available
    if 'k_xx' in d and d['k_xx'] is not None:
        valid_xx = np.isfinite(t) & np.isfinite(d['k_xx'])
        if np.count_nonzero(valid_xx) > 0:
            axL.plot(t[valid_xx], d['k_xx'][valid_xx], lw=2.0, marker='o', ms=4,
                     markevery=2, color=col)
    # k_yy (dashed)
    valid_yy = np.isfinite(t) & np.isfinite(d['k_yy'])
    if np.count_nonzero(valid_yy) > 0:
        axL.plot(t[valid_yy], d['k_yy'][valid_yy], lw=2.0, linestyle='--', marker='s', ms=4,
                 markevery=2, color=col)

axL.set_xlabel(f'Time ({unit})')
axL.set_ylabel(r'$k/k_0$' if normalize_values else r'$k$ (W/m·K)')
axL.set_title(r'$k_{xx}$ & $k_{yy}$ vs Time' if normalize_values else 'k vs Time')
beautify_axes(axL)

# Right: vs SSA (normalized SSA if normalize_values=True)
for idx, d in enumerate(collected_data):
    if 'ssa' not in d or d['ssa'] is None or d['ssa'].size == 0:
        continue
    col = color_for(idx)
    ssa = d['ssa']
    # k_xx (solid) if available
    if 'k_xx' in d and d['k_xx'] is not None:
        valid_xx = np.isfinite(ssa) & np.isfinite(d['k_xx'])
        if np.count_nonzero(valid_xx) > 0:
            axR.plot(ssa[valid_xx], d['k_xx'][valid_xx], lw=2.0, marker='o', ms=4,
                     markevery=2, color=col)
    # k_yy (dashed)
    valid_yy = np.isfinite(ssa) & np.isfinite(d['k_yy'])
    if np.count_nonzero(valid_yy) > 0:
        axR.plot(ssa[valid_yy], d['k_yy'][valid_yy], lw=2.0, linestyle='--', marker='s', ms=4,
                 markevery=2, color=col)

axR.set_xlabel(r'SSA/SSA$_0$' if normalize_values else 'SSA')
axR.set_ylabel(r'$k/k_0$' if normalize_values else r'$k$ (W/m·K)')
axR.set_title(r'$k_{xx}$ & $k_{yy}$ vs SSA' if normalize_values else 'k vs SSA')
beautify_axes(axR)

# Build a single figure-level legend: first entries describe components, then temperatures (colors)
from matplotlib.lines import Line2D
component_handles = [
    Line2D([0], [0], color='black', lw=2.0, linestyle='-', marker='o', ms=6, label=r'$k_{xx}$'),
    Line2D([0], [0], color='black', lw=2.0, linestyle='--', marker='s', ms=6, label=r'$k_{yy}$')
]
temp_handles = [
    Line2D([0], [0], color=color_for(i), lw=2.0, linestyle='-', label=f"{d['temperature']} K")
    for i, d in enumerate(collected_data)
]
all_handles = component_handles + temp_handles

fig.legend(all_handles, [h.get_label() for h in all_handles], loc='center left', bbox_to_anchor=(0.995, 0.5))
fig.subplots_adjust(right=0.92)  # reserve a slimmer margin for the legend
fig.tight_layout(rect=[0.0, 0.0, 0.92, 1.0], pad=0.2)

out_path = os.path.join(parent_dir, 'kxx_kyy_vs_ssa_and_time_all.png')
fig.savefig(out_path, dpi=dpi, bbox_inches='tight', transparent=True)
plt.close(fig)

print("\n✅ Generated outputs:")
print("  •", os.path.join(parent_dir, 'kyy_colormap_masked.png'))
if Nt > 0 and has_kxx:
    print("  •", os.path.join(parent_dir, 'kxx_colormap_masked.png'))
print("  •", os.path.join(parent_dir, 'kxx_kyy_vs_ssa_and_time_all.png'))