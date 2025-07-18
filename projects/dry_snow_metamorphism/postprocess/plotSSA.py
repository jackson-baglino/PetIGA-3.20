import matplotlib.pyplot as plt
import numpy as np
import os

# Read data from the file
filename = "./SSA_evo.dat"

# Read in environment variables
dim = os.getenv("dim")
inputFile = os.getenv("inputFile")
title = os.getenv("title")
print("inputFile: ", inputFile)

# Read in the grain data
grain_data = np.loadtxt(inputFile, delimiter=' ', usecols=(2,))

with open(filename, 'r') as file:
  input_data = file.read()
  # Parse the input data into a numpy array
  input_array = np.array([line.split() for line in input_data.split('\n') if line.strip()])

  # Convert the array elements to float
  input_array = input_array.astype(float)


  area_data   = input_array[:, 0]
  volume_data = input_array[:, 1]
  time_data_s = input_array[:, 2]  # time in seconds

  # Decide on time units
  max_time = np.max(time_data_s)

  if max_time < 3600:
      time_data = time_data_s / 60          # seconds → minutes
      time_units = "Minutes"
  elif max_time < 4 * 24 * 3600:
      time_data = time_data_s / 3600        # seconds → hours
      time_units = "Hours"
  else:
      time_data = time_data_s / (24 * 3600) # seconds → days
      time_units = "Days"

print(f"The number of time_data points is: {len(time_data)}")

# Calculate the initial area value
if dim == "2":
    area0 = np.sum(2*np.pi*grain_data)
    volume0 = np.sum(np.pi*grain_data**2)
elif dim == "3":
    area0 = np.sum(4*np.pi*grain_data**2)
    volume0 = np.sum(4/3*np.pi*grain_data**3)


# Normalize the area and volume data by the initial area and volume values
area_data = area_data / area0
volume_data = volume_data / volume0

# Compute the specific surface area (SSA) data
ssa_data = area_data / volume_data
ssa_data = ssa_data / ssa_data[0]

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(time_data, ssa_data, label='Surface Area Evolution')
# plt.loglog(time_data, ssa_data, label='Surface Area Evolution')
# Plot the slow
plt.xlabel('Time [hours]', fontsize=18)
plt.ylabel('Surface Area', fontsize=18)
plt.title('Surface Area Evolution', fontsize=24)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(False)

# Save the plot as an image file
outputfolder = os.getenv("outputfolder")
output_file = "ssa_evolution_plot_"+title+".png"
plt.savefig(output_file)

# Display the plot
plt.show(block=False)
plt.pause(5)
plt.close()