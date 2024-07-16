#!/bin/zsh

echo " "
echo "compiling"
echo " "
make NASAv2

# add name folder accordinglu --------------------------------------------------
title=getWet_GEOM_DIM_TIME
name=$title_$(date +%Y-%m-%d__%H.%M.%S)
dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2
folder=$dir/$name

if [ ! -d "$dir" ]; then
  mkdir -p "$dir"
fi

mkdir $folder/


# Define variable names to be exported -----------------------------------------
  # File names
input_dir="/Users/jacksonbaglino/PetIGA-3.20/demo/input/"

# inputFile=$input_dir"grainReadFile-2.dat"                 # 2 grains
inputFile=$input_dir"grainReadFile-5_s1-10.dat"             # 5 grains
# inputFile=$input_dir"grainReadFile-10_s1-10.dat"          # 10 grains
# inputFile=$input_dir"grainReadFile-37_s1-10_s2-21.dat"    # 37 grains
# inputFile=$input_dir"grainReadFile-88_s1-10_s2-21.dat"    # 88 grains
# inputFile=$input_dir"grainReadFile-135_s1-10_s2-21.dat"   # 135 grains
# inputFile=$input_dir"grainReadFile-165_s1-10_s2-30.dat"   # 165 grains

# Define simulation parameters -------------------------------------------------
# Define dimensions
dim=2

# Converty scientic notation to decimal using bc if needed
dim=$(echo "$dim" | bc -l)

# Domain sizes
# Lx=488.4e-6                    # Domain size X -- 2 Grain
# Ly=244.2e-6                    # Domain size Y -- 2 Grain
# Lz=244.2e-6                    # Domain size Z -- 2 Grain

Lx=0.35e-03                    # Domain size X -- 5 Grain
Ly=0.35e-03                    # Domain size Y -- 5 Grain
Lz=2.6424e-04                  # Domain size Z -- 5 Grain

# Lx=0.5e-03                    # Domain size X -- 10 Grain
# Ly=0.5e-03                    # Domain size Y -- 10 Grain
# Lz=1.101e-04                  # Domain size Z -- 10 Grain

# Lx=0.5e-03                    # Domain size X -- 10 Grain
# Ly=0.5e-03                    # Domain size Y -- 10 Grain
# Lz=2.202e-04                  # Domain size Z -- 10 Grain

# Lx=2.0e-3                     # Domain size X -- 88 Grain
# Ly=2.0e-3                     # Domain size Y -- 88 Grain
# Lz=2.509e-04                  # Domain size Z -- 88 Grain

# Lx=3.2e-3                     # Domain size X -- 135/165 Grain
# Ly=3.2e-3                     # Domain size Y -- 135/165 Grain
# Lz=1.0e-3                     # Domain size Z -- 135/165 Grain

# Convert scientific notation to decimal using bc
Lx=$(echo "$Lx" | bc -l)
Ly=$(echo "$Ly" | bc -l)
Lz=$(echo "$Lz" | bc -l)

# Number of elements
# Nx=264                        # Number of elements in X -- 2 Grain
# Ny=132                        # Number of elements in Y -- 2 Grain
# Nz=132                        # Number of elements in Z -- 2 Grain

Nx=190                        # Number of elements in X -- 5 Grain
Ny=190                        # Number of elements in Y -- 5 Grain
Nz=143                        # Number of elements in Z -- 5 Grain

# Nx=275                        # Number of elements in X -- 10 Grain
# Ny=275                        # Number of elements in Y -- 10 Grain
# Nz=122                        # Number of elements in Z -- 10 Grain

# Nx=264                        # Number of elements in X -- 2 Grain
# Ny=132                        # Number of elements in Y -- 2 Grain
# Nz=132                        # Number of elements in Z -- 2 Grain

# Nx=1100                       # Number of elements in X -- 88 Grain
# Ny=1100                       # Number of elements in Y -- 88 Grain
# Nz=138                        # Number of elements in Z -- 88 Grain

# Nx=1760                       # Number of elements in X -- 135/165 Grain
# Ny=1760                       # Number of elements in Y -- 135/165 Grain
# Nz=100                        # Number of elements in Z -- 135/165 Grain

# Time parameters
delt_t=1.0e-4                 # Time step
t_final=2*24*60*60            # Final time
n_out=10                      # Number of output files


# Convert scientific notation to decimal using bc
delt_t=$(echo "$delt_t" | bc -l)
t_final=$(echo "$t_final" | bc -l)
n_out=$(echo "$n_out" | bc -l)

# Other parameters
humidity=0.98                 # Relative humidity
temp=-3.0                    # Temperature

# Initial temperature gradients
grad_temp0X=0.0               # Initial temperature gradient X
grad_temp0Y=3.0            # Initial temperature gradient Y
grad_temp0Z=0.0               # Initial temperature gradient Z

# Convert scientific notation gradients to decimal using bc if needed
grad_temp0X=$(echo "$grad_temp0X" | bc -l)
grad_temp0Y=$(echo "$grad_temp0Y" | bc -l)
grad_temp0Z=$(echo "$grad_temp0Z" | bc -l)

# Export variables
export folder input_dir inputFile Lx Ly Lz Nx Ny Nz delt_t t_final n_out \
    humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim


# Copy files to folder ---------------------------------------------------------
cp NASAv2.c$folder
cp run_NASAv2.sh $folder

echo " "
echo "Calling ./NASAv2"
echo " "


# Run the simulation -----------------------------------------------------------
mpiexec -np 12 ./NASAv2 -initial_PFgeom -temp_initial -snes_rtol 1e-3 \
-snes_stol 1e-6 -snes_max_it 7 -ksp_gmres_restart 150 -ksp_max_it 1000  \
-ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
-snes_linesearch_type basic | \
tee /Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/outp.txt


# ------------------------------------------------------------------------------ 
# Move output file to folder ---------------------------------------------------
echo " "
echo "making directory" $folder
echo " "

mv $dir/outp.txt $folder
cp NASAv2.c $folder
cp run_plotNASAv2.sh $folder

cp plotSSA.py $folder
cp plotPorosity.py $folder

# Plot the results -------------------------------------------------------------
echo "Queing plotNASA.py"
./run_plotNASAv2.sh $name

# Create descriptive file ------------------------------------------------------
echo "----- SIMULATION PARAMETERS -----" > $folder/sim_params.dat
echo "Input file: $inputFile" >> $folder/sim_params.dat

echo " " >> $folder/sim_params.dat

echo "Dimensions:" >> $folder/sim_params.dat
echo "dim = $dim" >> $folder/sim_params.dat

echo " " >> $folder/sim_params.dat

echo "Domain sizes:" >> $folder/sim_params.dat
echo "Lx = $Lx" >> $folder/sim_params.dat
echo "Ly = $Ly" >> $folder/sim_params.dat
echo "Lz = $Lz" >> $folder/sim_params.dat

echo " " >> $folder/sim_params.dat


echo "Number of elements:" >> $folder/sim_params.dat
echo "Nx = $Nx" >> $folder/sim_params.dat
echo "Ny = $Ny" >> $folder/sim_params.dat
echo "Nz = $Nz" >> $folder/sim_params.dat

echo " " >> $folder/sim_params.dat

echo "Time parameters:" >> $folder/sim_params.dat
echo "delt_t = $delt_t" >> $folder/sim_params.dat
echo "t_final = $t_final" >> $folder/sim_params.dat

echo " " >> $folder/sim_params.dat

echo "State parameters:" >> $folder/sim_params.dat
echo "humidity = $humidity" >> $folder/sim_params.dat
echo "temp = $temp" >> $folder/sim_params.dat

echo " " >> $folder/sim_params.dat

echo "Initial temperature gradients:" >> $folder/sim_params.dat
echo "grad_temp0X = $grad_temp0X" >> $folder/sim_params.dat
echo "grad_temp0Y = $grad_temp0Y" >> $folder/sim_params.dat
echo "grad_temp0Z = $grad_temp0Z" >> $folder/sim_params.dat

echo "-------------------------------------------------------------------------"
echo "Done!"
echo " "

