import numpy as np
import matplotlib.pyplot as plt

# Define radii
R1 = 1.01e-4
R2 = 7.3e-5

# Define padding
a = 0.25 * R1
b = a

# Define domain size
Lx = 2 * (b + R1)
Ly = 2 * (a + R1 + R2)
Lz = Lx

print(f"Domain size: Lx = {Lx}, Ly = {Ly}")

# Define circle centers
x1 = Lx / 2
x2 = Lx / 2
y1 = a + R1
y2 = a + 2 * R1 + R2
z1 = x1
z2 = x2

# Write the geometry to a file (.dat)
file = 'grainReadFile-2G_Molaro_0p25R1'
filename = f'./inputs/{file}.dat'

with open(filename, 'w') as f:
  f.write(f"{Lx} {Ly} {Lz}\n")
  f.write(f"{x1} {y1} {z1} {R1}\n")
  f.write(f"{x2} {y2} {z2} {R2}\n")

import subprocess

# Optional: call environment generation script
input_file = filename
output_env_path = f'./configs/{file}.env'

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
