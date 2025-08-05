import numpy as np
import matplotlib.pyplot as plt
import os

# ----------------------
# Parameters
# ----------------------
seed = 10
np.random.seed(seed)

plotFlag = True
datFlag = True

N = 35
Lx, Ly = 1.0e-3, 1.0e-3
lx, ly = 1.0e-3, 1.0e-3
rad = 90e-6

# ----------------------
# Initialization
# ----------------------
radius = np.zeros(N)
xC = np.zeros(N)
yC = np.zeros(N)

# First circle
a = rad - 0.5 * rad
b = rad + 0.5 * rad
xC[0] = Lx / 2
yC[0] = Ly / 2
radius[0] = np.random.uniform(a, b)

# ----------------------
# Helper Functions
# ----------------------
def in_domain(x, y, R, Lx, Ly):
    return 0 <= x <= Lx and 0 <= y <= Ly

def overlaps(x, y, R, xC, yC, radius, num_grain):
    for j in range(num_grain):
        dist = np.sqrt((x - xC[j])**2 + (y - yC[j])**2)
        if dist < R + radius[j]:
            return True
    return False

# ----------------------
# Generate touching circles
# ----------------------
counts = 0
counts_max = int(1e8)
num_grain = 1
ref_grain = 0
num_attempts = 0
neighbor = np.zeros(N, dtype=int)

while counts < counts_max and num_grain < N:
    theta = 2 * np.pi * np.random.rand()
    R1 = np.random.uniform(a, b)
    centX = xC[ref_grain] + (R1 + radius[ref_grain]) * np.cos(theta)
    centY = yC[ref_grain] + (R1 + radius[ref_grain]) * np.sin(theta)

    if in_domain(centX, centY, R1, Lx, Ly) and not overlaps(centX, centY, R1, xC, yC, radius, num_grain):
        xC[num_grain] = centX
        yC[num_grain] = centY
        radius[num_grain] = R1
        neighbor[num_grain] = ref_grain
        ref_grain = num_grain
        num_grain += 1
        num_attempts = 0
    else:
        num_attempts += 1
        if num_attempts % 20 == 0:
            ref_grain = (ref_grain - 1) % num_grain

    counts += 1

if num_grain < N:
    print(f"⚠️ Only {num_grain} particles initialized (target was {N})")

# ----------------------
# Post-processing
# ----------------------
radius = radius[:num_grain]
xC = xC[:num_grain] + (lx - Lx) / 2
yC = yC[:num_grain] + (ly - Ly) / 2

Lz = 2 * np.max(radius) + (lx - Lx)
porosity = 1 - np.pi * np.sum(radius**2) / (Lx * Ly)

# ----------------------
# Plotting
# ----------------------
if plotFlag:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    for (x, y, r) in zip(xC, yC, radius):
        circle = plt.Circle((x, y), r, color="#D95319", fill=False)
        ax.add_patch(circle)
    for i, (x, y) in enumerate(zip(xC, yC)):
        ax.text(x, y, str(i+1), fontsize=8, ha='center', va='center')
    ax.set_title(f"$\\phi$ = {porosity:.3f}, N = {num_grain}", fontsize=16)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

# ----------------------
# Save to file
# ----------------------
if datFlag:
    import subprocess

    os.makedirs("geomFiles", exist_ok=True)
    filename = f"grainReadFile-{num_grain}_s1-{seed}"
    output = np.column_stack([xC, yC, (Lz / 2) * np.ones_like(radius), radius])
    filepath = f"inputs/{filename}.dat"
    with open(filepath, 'w') as f:
        f.write(f"{Lx:.6e} {Ly:.6e} {Lz:.6e}\n")
        np.savetxt(f, output, fmt='%.6e')

    input_file = filepath
    output_env_path = f"configs/{filename}.env"

    cmd = [
        'python',
        'scripts/generate_env_from_input.py',
        input_file,
        output_env_path,
        str(Lx),
        str(Ly)
    ]
    print("[INFO] Running generate_env_from_input.py...")
    subprocess.run(cmd, check=True)