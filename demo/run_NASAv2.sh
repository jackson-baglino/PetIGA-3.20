#!/bin/zsh

# Function to create a timestamped results folder
create_folder() {
    name="$title$(date +%Y-%m-%d__%H.%M.%S)"
    dir="/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2"
    folder="$dir/$name"

    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
    fi

    mkdir -p "$folder"
}

# Function to compile the code
compile_code() {
    echo "Compiling..."
    make NASAv2
}

# Function to set simulation parameters
set_parameters() {
    input_dir="/Users/jacksonbaglino/PetIGA-3.20/demo/input/"

    # Select input file
    inputFile="${input_dir}${filename}"  # Default file

    # Set the domain sizes and number of elements based on the input file ------
    if [[ $inputFile == *"grainReadFile-2.dat"* ]]; then
        Lx=488.4e-6
        Ly=244.2e-6
        Lz=244.2e-6

        Nx=269
        Ny=135
        Nz=135

        eps=9.096e-07

    elif [[ $inputFile == *"grainReadFile-2_Molaro.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        Lx=3.0300e-04
        Ly=3.8280e-04
        Lz=3.0300e-04

        Nx=167
        Ny=211
        Nz=167

        eps=9.096e-07

    elif [[ $inputFile == *"grainReadFile-10_s1-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        Lx=0.5e-03
        Ly=0.5e-03
        Lz=2.202e-04

        Nx=275
        Ny=275
        Nz=122

        eps=9.096e-07

    elif [[ $inputFile == *"grainReadFile-27_MOLARO_s2-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        Lx=0.75e-03
        Ly=0.75e-03 
        Lz=0.000242175903182621 

        Nx=413
        Ny=413
        Nz=134

        eps=9.096e-07

    elif [[ $inputFile == *"grainReadFile_3D-30_s1-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        # Check that dim = 3
        if [[ $dim -ne 3 ]]; then
            echo "Error: Dimension mismatch. Expected dim = 3 for input file: $inputFile"
            exit 1
        fi
        
        Lx=0.5e-03
        Ly=0.5e-03
        Lz=0.5e-03

        Nx=275
        Ny=275
        Nz=275

        eps=9.096e-07

    else
        echo "Error: Unknown input file."
        exit 1
    fi

    # Exporting variables after determining the configuration
    export folder input_dir inputFile title Lx Ly Lz Nx Ny Nz delt_t t_final n_out \
        humidity temp grad_temp0X grad_temp0Y grad_temp0Z dim eps
}

# Function to run the simulation
run_simulation() {
    echo "Running simulation..."
    mpiexec -np 12 ./NASAv2 -initial_PFgeom -temp_initial -snes_rtol 1e-3 \
    -snes_stol 1e-6 -snes_max_it 7 -ksp_gmres_restart 150 -ksp_max_it 1000 \
    -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor \
    -snes_linesearch_type basic | tee $folder/outp.txt
}

# Function to copy files and generate the descriptive file
finalize_results() {
    echo "Finalizing results..."
    
    # Copy necessary files to results folder
    cp NASAv2.c run_NASAv2.sh plotNASA.py plotSSA.py plotPorosity.py $folder
    
    # Save simulation parameters
    cat << EOF > $folder/sim_params.dat
----- SIMULATION PARAMETERS -----
Input file: $inputFile

Dimensions:
dim = $dim

Interface width:
eps = $eps

Domain sizes:
Lx = $Lx
Ly = $Ly
Lz = $Lz

Number of elements:
Nx = $Nx
Ny = $Ny
Nz = $Nz

Time parameters:
delt_t = $delt_t
t_final = $t_final

State parameters:
humidity = $humidity
temp = $temp

Initial temperature gradients:
grad_temp0X = $grad_temp0X
grad_temp0Y = $grad_temp0Y
grad_temp0Z = $grad_temp0Z
EOF
}

# Function to run plotting scripts
run_plotting() {
    echo "Queuing plotNASA.py"
    ./run_plotNASAv2.sh $name
}

# Main execution starts here
echo " "
echo "Starting NASAv2 simulation workflow"
echo " "

# Define filename and title
filename="grainReadFile_3D-30_s1-10.dat"
title="NASAv2-TEST__"


compile_code
create_folder

# Define default time and physical parameters here
delt_t=1.0e-4
t_final=1
n_out=2

t_final=$(echo "$t_final" | bc -l)

humidity=0.98
temp=-20.0
grad_temp0X=0.0
grad_temp0Y=3.0
grad_temp0Z=0.0

dim=3
filename="grainReadFile-2.dat"
title="NASAv2-TEST__"

set_parameters
run_simulation
finalize_results
run_plotting

echo "-------------------------------------------------------------------------"
echo "Done!"
echo " "