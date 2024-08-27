#!/bin/zsh

echo " "
echo "compiling"
echo " "

make Metamorph

# add name folder accordingly --------------------------------------------------
title=dry_2G_2D_24h_
name=$title$(date +%Y-%m-%d__%H.%M.%S)
dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/wetamorph
folder=$dir/$name

mkdir -p "$dir"
mkdir $folder/


# Define variable names to be exported -----------------------------------------
  # File names

  # Simulation parameters
  angl=0                      # contact angle: 0->real  1->120degrees
  aa=0		                    # 1->alpha_j=0
  mm=0		                    # 1->constant mobility

  # Define domain sizes

  # Number of elements

  # Define simulation parameters

  # Export variables
  export folder angl aa mm


# Copy files to folder
cp ~/PetIGA/demo/Metamorph.c $folder/
cp ~/PetIGA/demo/run_metam.sh $folder/


# Run simulation
echo " "
echo "running Metamorph"
echo " "

mpirun -np 12 ./Metamorph -initial_PFgeom -temp_initial -snes_rtol 1e-3 \
-snes_stol 1e-6 -snes_max_it 7 -ksp_gmres_restart 150 -ksp_max_it 1000 \
-ksp_converged_maxits 1 -ksp_converged_reason -snes_converged_reason  \
-snes_linesearch_monitor -snes_linesearch_type basic | tee $folder/outp.txt


# Move files to folder
cp Metamorph.c $folder/
cp run_metam.sh $folder/
# cp run_plotMetam.sh $folder/                # Need to create this file

cp plotSSA.py $folder/
cp plotPorosity.py $folder/

# echo "Queing plotMetam.py"                  # Need to create this file 
# ./run_plotMetam.sh $name


# Create descriptive file
  # Filename
  echo "----- SIMULATION PARAMETERS -----" > $folder/sim_params.dat
  echo "Input file: $inputFile" >> $folder/sim_params.dat
  echo " " >> $folder/sim_params.dat

  # Dimensions
  echo "Dimensions:" >> $folder/sim_params.dat
  echo "dim = $dim" >> $folder/sim_params.dat
  echo " " >> $folder/sim_params.dat

  # Interface width
  echo "Interface wiedth:" >> $folder/sim_params.dat
  echo "eps = $eps" >> $folder/sim_params.dat
  echo " " >> $folder/sim_params.dat

  # Domain sizes
  echo "Domain sizes:" >> $folder/sim_params.dat
  echo "Lx = $Lx" >> $folder/sim_params.dat
  echo "Ly = $Ly" >> $folder/sim_params.dat
  echo "Lz = $Lz" >> $folder/sim_params.dat
  echo " " >> $folder/sim_params.dat

  # Number of elements
  echo "Number of elements:" >> $folder/sim_params.dat
  echo "Nx = $Nx" >> $folder/sim_params.dat
  echo "Ny = $Ny" >> $folder/sim_params.dat
  echo "Nz = $Nz" >> $folder/sim_params.dat
  echo " " >> $folder/sim_params.dat

  # Time parameters
  echo "Time parameters:" >> $folder/sim_params.dat
  echo "delt_t = $delt_t" >> $folder/sim_params.dat
  echo "t_final = $t_final" >> $folder/sim_params.dat
  echo " " >> $folder/sim_params.dat

  # State parameters
  echo "State parameters:" >> $folder/sim_params.dat
  echo "humidity = $humidity" >> $folder/sim_params.dat
  echo "temp = $temp" >> $folder/sim_params.dat
  echo " " >> $folder/sim_params.dat

  # Initial temperature gradients
  echo "Initial temperature gradients:" >> $folder/sim_params.dat
  echo "grad_temp0X = $grad_temp0X" >> $folder/sim_params.dat
  echo "grad_temp0Y = $grad_temp0Y" >> $folder/sim_params.dat
  echo "grad_temp0Z = $grad_temp0Z" >> $folder/sim_params.dat

echo "-------------------------------------------------------------------------"
echo "Done!"
echo " "
