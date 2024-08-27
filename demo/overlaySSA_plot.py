import matplotlib.pyplot as plt
import numpy as np
import os

# Read data from the file
filename2D = os.getenv("filename2D")
filename3D = os.getenv("filename3D")

inputFile = os.getenv("inputFile")
output_file = os.getenv("output_file")

print("filename2D: ", filename2D)
print("filename3D: ", filename3D)


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


# Read in the SSA data for 3D
with open(filename3D, 'r') as file:
  input_data3D = file.read()
  # Parse the input data into a numpy array
  input_array3D = np.array([line.split() for line in input_data3D.split('\n') if line.strip()])

  # Convert the array elements to float
  input_array3D = input_array3D.astype(float)


# Unpack the input array into separate arrays
area_data2D       = input_array2D[:, 0]
volume_data2D     = input_array2D[:, 1]
time_data2D       = input_array2D[:, 2]/60/60   # Convert time from seconds to hours

# Unpack the input array into separate arrays
area_data3D       = input_array3D[:, 0]
volume_data3D     = input_array3D[:, 1]
time_data3D       = input_array3D[:, 2]/60/60   # Convert time from seconds to hours

# Calculate the initial area value
area02D = np.sum(2*np.pi*grain_data)
volume02D = np.sum(np.pi*grain_data**2)

area03D = np.sum(4*np.pi*grain_data**2)
volume03D = np.sum(4/3*np.pi*grain_data**3)

# Normalize the area and volume data by the initial area and volume values
area_data2D = area_data2D / area02D
volume_data2D = volume_data2D / volume02D

area_data3D = area_data3D / area03D
volume_data3D = volume_data3D / volume03D

# Compute the specific surface area (SSA) data
ssa_data2D = area_data2D / volume_data2D
ssa_data2D = ssa_data2D / ssa_data2D[0]

ssa_data3D = area_data3D / volume_data3D
ssa_data3D = ssa_data3D / ssa_data3D[0]


# Plotting the data
plt.figure(figsize=(10, 6))
# plt.loglog(time_data2D, ssa_data2D, label='Surface Area Evolution (2D)', linestyle='dashed')
# plt.loglog(time_data3D, ssa_data3D, label='Surface Area Evolution (3D)', linestyle='dotted')
plt.plot(time_data2D, ssa_data2D, label='Surface Area Evolution (2 Grain - 2D)', linestyle='dashed')
plt.plot(time_data3D, ssa_data3D, label='Surface Area Evolution (2 Grain - 3D)', linestyle='dotted')

# indices = np.arange(0, len(time_data2D), 5)
# plt.scatter(time_data2D[indices], ssa_data2D[indices], marker='o')


plt.xlabel('Time [hours]', fontsize=18)
plt.ylabel('Normalized Specific Surface Area [-]', fontsize=18)
plt.title('Specific Surface Area Evolution', fontsize=24)
plt.legend(fontsize=14)

print("The log slope for 2D is: ", np.polyfit(np.log(time_data2D[1:]), np.log(ssa_data2D[1:]), 1)[0])
print("The log slope for 3D is: ", np.polyfit(np.log(time_data3D[1:]), np.log(ssa_data3D[1:]), 1)[0])

plt.tick_params(axis='both', which='major', labelsize=14)
plt.grid(False)

# Save the plot as an image file
plt.savefig(output_file)

# Display the plot
plt.show()
plt.close()
