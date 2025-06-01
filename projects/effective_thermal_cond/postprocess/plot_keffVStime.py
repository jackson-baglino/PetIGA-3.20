import matplotlib.pyplot as plt
import csv

# ----------------------------
# File paths
# ----------------------------
keff_file = "/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/outputs/homog/ThermalSim_homog_2025-05-28__12.56.48/keff_summary.csv"
ssa_file = "/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/NASAv2_10G_2D_T-20.0_hum0.70_2025-03-13__14.20.59/SSA_evo.dat"

# ----------------------------
# Load SSA evolution data (no header, 4 columns)
# ----------------------------
ssa_by_step = {}  # step -> (ice_vol, ssa, time_sec)
with open(ssa_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        ice_vol, ssa, time_sec, step = map(float, parts)
        ssa_by_step[int(step)] = (ice_vol, ssa, time_sec)

# ----------------------------
# Load k_eff and match with SSA and time
# ----------------------------
matched_data = []  # list of (time_sec, ssa, k_xx, k_yy)

with open(keff_file, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        fname = row['filename']
        step = int(fname.split('_')[-1].split('.')[0])
        if step in ssa_by_step:
            _, ssa_val, time_sec = ssa_by_step[step]
            k_xx = float(row['k_eff_xx'])
            k_yy = float(row['k_eff_yy'])
            matched_data.append((time_sec, ssa_val, k_xx, k_yy))

# Unzip for plotting by time
times, ssas, kxxs, kyys = zip(*matched_data)

# ----------------------------
# Plot k_eff vs time
# ----------------------------
plt.figure(figsize=(10, 5))
plt.plot(times, kxxs, label="k_eff_xx")
plt.plot(times, kyys, label="k_eff_yy")
plt.xlabel("Time (s)")
plt.ylabel("Effective conductivity")
plt.title("Effective Thermal Conductivity vs. Time")
plt.legend()
plt.grid(True)
plt.tight_layout()

# ----------------------------
# Plot k_eff vs SSA (sorted)
# ----------------------------
# Sort by SSA
sorted_data = sorted(zip(ssas, kxxs, kyys))
ssas_sorted, kxx_sorted, kyy_sorted = zip(*sorted_data)

plt.figure(figsize=(10, 5))
plt.plot(ssas_sorted, kxx_sorted, label="k_eff_xx")
plt.plot(ssas_sorted, kyy_sorted, label="k_eff_yy")
plt.xlabel("Specific Surface Area (SSA)")
plt.ylabel("Effective conductivity")
plt.title("Effective Thermal Conductivity vs. SSA")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.show()