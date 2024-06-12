import matplotlib.pyplot as plt
import numpy as np
import os

# Read data from the file
filename = "SSA_evo.dat"

# Read in environment variables
dim = os.getenv("dim")

# inputFile = os.getenv("inputFile")
inputFile = "/Users/jacksonbaglino/PetIGA-3.20/demo/input/grainReadFile-165_s1-10_s2-30.dat"

print(f"inputFile = {inputFile}")

# Read in the grain data
grain_data = np.loadtxt(inputFile, delimiter=' ', usecols=(2,))

with open(filename, 'r') as file:
  input_data = file.read()
  # Parse the input data into a numpy array
  input_array = np.array([line.split() for line in input_data.split('\n') if line.strip()])
  # Convert the array elements to float
  input_array = input_array.astype(float)

ssa_data = input_array[:, 0]

SSA0 = np.sum(2*np.pi*grain_data)

# Normalize SSA data by the third column of grain data
c = ssa_data[0] / SSA0
if c == 0:
    c = 1

ssa_data = ssa_data / c

# print(f"ssa_data = {ssa_data}")

normalized_ssa_data = ssa_data
normalized_ssa_data = normalized_ssa_data[1:]

# Create a time array (assuming data points are equally spaced)
# time = range(1, len(normalized_ssa_data) + 1)
time = np.linspace(1, 18, len(normalized_ssa_data))

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(time, normalized_ssa_data, label='Surface Area Evolution')
plt.xlabel('Time [days]', fontsize=18)
plt.ylabel('Surface Area', fontsize=18)
plt.title('Surface Area Evolution', fontsize=24)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(False)


# Save the plot as an image file
output_file = "ssa_evolution_plot.png"
plt.savefig(output_file)

# Display the plot
plt.show()
plt.close()
