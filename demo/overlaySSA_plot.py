import matplotlib.pyplot as plt
import numpy as np
import os

# Read data from the file
filename2D = os.getenv("filename2D")
filename3D = os.getenv("filename3D")

# Read in environment variables
dim = os.getenv("dim")
inputFile = "/Users/jacksonbaglino/PetIGA-3.20/demo/input/grainReadFile-2.dat"
print("inputFile: ", inputFile)

# Read in the grain data
grain_data = np.loadtxt(inputFile, delimiter=' ', usecols=(2,))

# Read in the SSA data for 2D
with open(filename2D, 'r') as file:
  input_data2D = file.read()
  # Parse the input data into a numpy array
  input_array2D = np.array([line.split() for line in input_data2D.split('\n') if line.strip()])
  # Convert the array elements to float
  input_array2D = input_array2D.astype(float)

ssa_data2D = input_array2D[:, 0]

# Read in the SSA data for 3D
with open(filename3D, 'r') as file:
  input_data3D = file.read()
  # Parse the input data into a numpy array
  input_array3D = np.array([line.split() for line in input_data3D.split('\n') if line.strip()])
  # Convert the array elements to float
  input_array3D = input_array3D.astype(float)

ssa_data3D = input_array3D[:, 0]

SSA02D = np.sum(2*np.pi*grain_data)
SSA03D = np.sum(4*np.pi*grain_data)

# Normalize the data (2D)
ssa_data2D = ssa_data2D/ssa_data2D[0]

normalized_ssa_data2D = ssa_data2D
normalized_ssa_data2D = normalized_ssa_data2D[0:]

# Normalize the data (3D)
ssa_data3D = ssa_data3D/ssa_data3D[0]

normalized_ssa_data3D = ssa_data3D
normalized_ssa_data3D = normalized_ssa_data3D[0:]

# Create a time array (assuming data points are equally spaced)
time2D = input_array2D[:, 2]/60/60
time3D = input_array3D[:, 2]/60/60

# Plotting the data
plt.figure(figsize=(10, 6))
# plt.loglog(time2D, normalized_ssa_data2D, label='Surface Area Evolution (2D)')
# plt.loglog(time3D, normalized_ssa_data3D, label='Surface Area Evolution (3D)')
plt.plot(time2D, normalized_ssa_data2D, label='Surface Area Evolution (2 Grain - 3D)')
plt.plot(time3D, normalized_ssa_data3D, label='Surface Area Evolution (10 Grain - 3D)')
plt.xlabel('Time [hours]', fontsize=18)
plt.ylabel('Surface Area', fontsize=18)
plt.title('Surface Area Evolution', fontsize=24)
plt.legend(fontsize=14)

# print("The log slope for 2D is: ", np.polyfit(np.log(time2D), np.log(normalized_ssa_data2D), 1)[0])
# print("The log slope for 3D is: ", np.polyfit(np.log(time3D), np.log(normalized_ssa_data3D), 1)[0])

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(False)

# Save the plot as an image file
output_file = "ssa_evolution_plot.png"
plt.savefig(output_file)

# Display the plot
# plt.show(block=False)
plt.show()
# plt.pause(10)
plt.close()
