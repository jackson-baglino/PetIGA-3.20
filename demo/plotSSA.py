import matplotlib.pyplot as plt
import numpy as np

# Read data from the file
filename = "SSA_evo.dat"
grainFile = "/Users/jacksonbaglino/PetIGA/demo/input/grainReadFile.csv"

try:
    with open(filename, 'r') as file:
        # Assuming the file contains one column of numbers representing SSA values over time
        ssa_data = [float(line.strip()) for line in file.readlines()]
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    exit()
except Exception as e:
    print(f"Error reading data from '{filename}': {e}")
    exit()

try:
    # Assuming grainFile is a CSV file with 150 rows and 4 columns
    grain_data = np.loadtxt(grainFile, delimiter=',', usecols=(2,))
except FileNotFoundError:
    print(f"Error: File '{grainFile}' not found.")
    exit()
except Exception as e:
    print(f"Error reading data from '{grainFile}': {e}")
    exit()

grain_data = grain_data*np.sqrt(2*(2e-3)**2)/np.sqrt(2*(200)**2)

SSA0 = np.sum(2*np.pi*grain_data)
print(f"SSA0 = {SSA0}")

# Normalize SSA data by the third column of grain data
c = ssa_data[0] / SSA0
ssa_data = ssa_data / c

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
