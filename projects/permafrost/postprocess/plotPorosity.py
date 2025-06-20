import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot SSA and porosity evolution from simulation results."
    )
    parser.add_argument("ssa_file", help="Path to SSA_evo.dat file")
    parser.add_argument("grain_file", help="Path to grainReadFile.csv")
    parser.add_argument("--outdir", default=".", help="Directory to save plots (default: current dir)")
    return parser.parse_args()

def load_ssa_and_porosity(filename):
    try:
        data = np.genfromtxt(filename, delimiter=' ', dtype=float, comments='#', filling_values=np.nan)
        data = data[~np.isnan(data).any(axis=1)]  # Drop rows with NaNs
        ssa = data[200:-1, 0]
        porosity = data[200:-1, 1]
        return ssa, porosity
    except Exception as e:
        print(f"[ERROR] Could not read SSA file: {e}")
        sys.exit(1)

def load_grain_data(filename):
    try:
        grain_data = np.loadtxt(filename, delimiter=',', usecols=(2,))
        return grain_data
    except Exception as e:
        print(f"[ERROR] Could not read grain file: {e}")
        sys.exit(1)

def normalize_ssa(ssa_raw, grain_data):
    # Normalize SSA based on analytical surface area from grains
    effective_grain_radius = np.sqrt(2 * (2e-3)**2)
    scaled_grains = grain_data * effective_grain_radius / np.sqrt(2 * 200**2)
    SSA0 = np.sum(2 * np.pi * scaled_grains)
    scale_factor = ssa_raw[0] / SSA0
    return ssa_raw / scale_factor / 10  # Final units in [m]

def normalize_porosity(porosity_raw):
    domain_volume = 2e-3 * 1.583e-3
    return porosity_raw / domain_volume

def main():
    args = parse_args()
    ssa_raw, porosity_raw = load_ssa_and_porosity(args.ssa_file)
    grain_data = load_grain_data(args.grain_file)

    normalized_ssa = normalize_ssa(ssa_raw, grain_data)
    normalized_porosity = normalize_porosity(porosity_raw)

    time = np.linspace(1, 3, len(normalized_ssa))  # Placeholder: adjust if needed

    os.makedirs(args.outdir, exist_ok=True)

    # Create dual-axis plot
    plt.figure(figsize=(10, 6))
    ax1 = plt.gca()
    line1, = ax1.plot(time, normalized_ssa, 'b-', label='Surface Area Evolution')
    ax1.set_xlabel('Time [days]', fontsize=18)
    ax1.set_ylabel('Surface Area [m]', fontsize=18, color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(False)
    ax1.set_title('Evolution of Surface Area and Porosity', fontsize=24)

    ax2 = ax1.twinx()
    line2, = ax2.plot(time, normalized_porosity, 'r-', label='Porosity Evolution')
    ax2.set_ylabel('Porosity', fontsize=18, color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    ax1.tick_params(axis='both', labelsize=14)
    ax2.tick_params(axis='both', labelsize=14)

    # Save figure
    fig_path = os.path.join(args.outdir, "combined_evolution_dual_y_axis.png")
    plt.savefig(fig_path, bbox_inches="tight")
    print(f"[INFO] Plot saved to: {fig_path}")
    plt.show()

if __name__ == "__main__":
    main()