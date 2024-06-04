#!/bin/zsh

echo " "
echo "compiling"
echo " "
make NASAv2

# add name folder accordingly:    <---------------------------------------------
name=res_$(date +%Y-%m-%d__%H.%M.%S)
dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2
folder=$dir/$name

if [ ! -d "$dir" ]; then
  mkdir -p "$dir"
fi

mkdir $folder/

# Define variable names to be exported  <---------------------------------------
  # File names
input_dir="/Users/jacksonbaglino/PetIGA-3.20/demo/input/"
filename=$input_dir"grainReadFile-165_s1-10_s2-30.dat"


# Simulation parameters
# Domain sizes
Lx=3.2e-03                    # Domain size X
Ly=3.2e-03                    # Domain size Y
Lz=1.35e-04                   # Domain size Z

# Convert scientific notation to decimal using bc
Lx=$(echo "$Lx" | bc -l)
Ly=$(echo "$Ly" | bc -l)
Lz=$(echo "$Lz" | bc -l)

# Number of elements
Nx=440                        # Number of elements in X
Ny=440                        # Number of elements in Y
Nz=75                         # Number of elements in Z

# Time parameters
delt_t=1.0e-4                 # Time step
t_final=1.0e-3                # Final time

# Convert scientific notation to decimal using bc
delt_t=$(echo "$delt_t" | bc -l)
t_final=$(echo "$t_final" | bc -l)

# Other parameters
humidity=0.98                 # Relative humidity
temp=-30.0                    # Temperature

# Initial temperature gradients
grad_temp0X=0.0               # Initial temperature gradient X
grad_temp0Y=3.0               # Initial temperature gradient Y
grad_temp0Z=0.0               # Initial temperature gradient Z

# Convert scientific notation gradients to decimal using bc if needed
grad_temp0X=$(echo "$grad_temp0X" | bc -l)
grad_temp0Y=$(echo "$grad_temp0Y" | bc -l)
grad_temp0Z=$(echo "$grad_temp0Z" | bc -l)

# Define dimensions
dim=3

# Converty scientic notation to decimal using bc if needed
dim=$(echo "$dim" | bc -l)

export folder input_dir filename Lx Ly Lz Nx Ny Nz delt_t t_final humidity \
  temp grad_temp0X grad_temp0Y grad_temp0Z dim


cp NASAv2.c$folder
cp run_NASAv2.sh $folder

echo " "
echo "Calling ./NASAv2"
echo " "


mpiexec -np 12 ./NASAv2 -initial_PFgeom -temp_initial -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 -ksp_gmres_restart 150 -ksp_max_it 1000  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor -snes_linesearch_type basic | tee /Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/outp.txt


echo " "
echo "making directory" $folder
echo " "


# mv $dir/*.dat $folder
mv $dir/outp.txt $folder
cp NASAv2.c $folder
cp run_plotNASAv2.sh $folder


echo "Queing plotNASA.py"
./run_plotNASAv2.sh $name


echo "-------------------------------------------------------------------------"
echo "Done!"
echo " "

