import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def read_environment():
    """Fetch required environment variables."""
    dim = os.getenv("dim")
    input_file = os.getenv("inputFile")
    title = os.getenv("title", "untitled")
    output_folder = os.getenv("outputfolder", ".")
    
    if not dim or not input_file:
        print("[ERROR] Required environment variables 'dim' or 'inputFile' not set.")
        sys.exit(1)
    
    print(f"[INFO] Input grain file: {input_file}")
    return dim, input_file, title, output_folder

def load_grain_data(filepath):
    """Load grain radius data."""
    try:
        return np.loadtxt(filepath, delimiter=' ', usecols=(2,))
    except Exception as e:
        print(f"[ERROR] Failed to load grain data: {e}")
        sys.exit(1)

def load_ssa_data(filepath):
    """Load SSA and volume evolution data."""
    try:
        raw_data = np.genfromtxt(filepath, comments='#')
        if raw_data.ndim != 2 or raw_data.shape[1] < 3:
            raise ValueError("SSA data must have at least 3 columns.")
        area = raw_data[:, 0]
        volume = raw_data[:, 1]
        time = raw_data[:, 2] / 3600  # Convert seconds to hours
        print(f"[INFO] Loaded {len(time)} time steps from SSA data.")
        return area, volume, time
    except Exception as e:
        print(f"[ERROR] Failed to load SSA data: {e}")
        sys.exit(1)

def compute_normalized_ssa(area_data, volume_data, grain_data, dim):
    """Compute and normalize SSA data from grain and simulation values."""
    if dim == "2":
        area0 = np.sum(2 * np.pi * grain_data)
        volume0 = np.sum(np.pi * grain_data**2)
    elif dim == "3":
        area0 = np.sum(4 * np.pi * grain_data**2)
        volume0 = np.sum((4 / 3) * np.pi * grain_data**3)
    else:
        raise ValueError(f"Invalid dimension: {dim}. Expected '2' or '3'.")
    
    area_norm = area_data / area0
    volume_norm = volume_data / volume0
    ssa = area_norm / volume_norm
    return ssa / ssa[0]

def plot_ssa(time, ssa, title, output_path):
    """Plot and save SSA evolution figure."""
    plt.figure(figsize=(10, 6))
    plt.plot(time, ssa, label='Normalized SSA', color='blue')
    
    plt.xlabel('Time [hours]', fontsize=18)
    plt.ylabel('Normalized SSA', fontsize=18)
    plt.title('Surface Area Evolution', fontsize=24)
    
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(False)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches="tight")
    print(f"[INFO] Plot saved to {output_path}")

    plt.show(block=False)
    plt.pause(5)
    plt.close()

def main():
    ssa_file = "./SSA_evo.dat"
    
    dim, input_file, title, output_folder = read_environment()
    grain_data = load_grain_data(input_file)
    area_data, volume_data, time_data = load_ssa_data(ssa_file)
    ssa = compute_normalized_ssa(area_data, volume_data, grain_data, dim)

    output_path = os.path.join(output_folder, f"ssa_evolution_plot_{title}.png")
    plot_ssa(time_data, ssa, title, output_path)

if __name__ == "__main__":
    main()