#!/bin/zsh

################################################################################
############################### DEFINE FUNCTIONS ###############################
################################################################################

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

# Function to write parameters to CSV file
write_parameters_to_csv() {
    csv_file="$folder/simulation_parameters.csv"
    echo "Variable,Value" > "$csv_file"
    echo "folder,$folder" >> "$csv_file"
    echo "input_dir,$input_dir" >> "$csv_file"
    echo "inputFile,$inputFile" >> "$csv_file"
    echo "title,$title" >> "$csv_file"
    echo "Lx,$Lx" >> "$csv_file"
    echo "Ly,$Ly" >> "$csv_file"
    echo "Lz,$Lz" >> "$csv_file"
    echo "Nx,$Nx" >> "$csv_file"
    echo "Ny,$Ny" >> "$csv_file"
    echo "Nz,$Nz" >> "$csv_file"
    echo "delt_t,$delt_t" >> "$csv_file"
    echo "t_final,$t_final" >> "$csv_file"
    echo "n_out,$n_out" >> "$csv_file"
    echo "humidity,$humidity" >> "$csv_file"
    echo "temp,$temp" >> "$csv_file"
    echo "grad_temp0X,$grad_temp0X" >> "$csv_file"
    echo "grad_temp0Y,$grad_temp0Y" >> "$csv_file"
    echo "grad_temp0Z,$grad_temp0Z" >> "$csv_file"
    echo "dim,$dim" >> "$csv_file"
    echo "eps,$eps" >> "$csv_file"
}

# Function to set simulation parameters
set_parameters() {
    input_dir="/Users/jacksonbaglino/PetIGA-3.20/demo/input/"

    # Select input file
    inputFile="${input_dir}${filename}"  # Default file

    # Copy inputFile to results folder
    cp $inputFile $folder

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
        Lx=0.0002424
        Ly=0.0003884
        Lz=0.0002424

        Nx=134
        Ny=214
        Nz=134

        eps=9.096e-07
    elif [[ $inputFile == *"grainReadFile-2_Molaro_0p05.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "
        Lx=0.0002121;
        Ly=0.0003581;
        Lz=0.0002121;

        # Nx=117
        # Ny=196
        # Nz=117

        Nx=234
        Ny=392
        Nz=234

        eps=9.096e-07
    elif [[ $inputFile == *"grainReadFile-2_Molaro_0p2.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "
        Lx=0.0002424;
        Ly=0.0003884;
        Lz=0.0002424;

        # Nx=134
        # Ny=214
        # Nz=134

        Nx=268
        Ny=428
        Nz=268

        eps=9.096e-07
    elif [[ $inputFile == *"grainReadFile-2_Molaro_0p3.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "
        Lx=0.0002828;
        Ly=0.0004288;
        Lz=0.0002828;

        Nx=290
        Ny=450
        Nz=290

        eps=9.096e-07

    elif [[ $inputFile == *"grainReadFile-2_Molaro_0p5.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "
        Lx=0.000303;
        Ly=0.000449;
        Lz=0.000303;

        # Nx=167
        # Ny=247
        # Nz=167

        Nx=334
        Ny=494
        Nz=334

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

    elif [[ $inputFile == *"grainReadFile_3D-42_s1-10.dat"* ]]; then
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

        Nx=160
        Ny=160
        Nz=160

        eps=1.56979924263831e-06

    elif [[ $inputFile == *"grainReadFile_3D-175_s1-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        # Check that dim = 3
        if [[ $dim -ne 3 ]]; then
            echo "Error: Dimension mismatch. Expected dim = 3 for input file: $inputFile"
            exit 1
        fi
        
        Lx=0.75e-03
        Ly=2.0e-03
        Lz=0.75e-03

        Nx=239
        Ny=678
        Nz=239

        eps=1.56979924263831e-06

    elif [[ $inputFile == *"grainReadFile_3D-500_s1-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        # Check that dim = 3
        if [[ $dim -ne 3 ]]; then
            echo "Error: Dimension mismatch. Expected dim = 3 for input file: $inputFile"
            exit 1
        fi
        
        Lx=1.5e-03
        Ly=1.5e-03
        Lz=1.5e-03

        Nx=478
        Ny=478
        Nz=478

        eps=1.56979924263831e-06

    elif [[ $inputFile == *"grainReadFile-118_s1-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        # Check that dim = 2
        if [[ $dim -ne 2 ]]; then
            echo "Error: Dimension mismatch. Expected dim = 2 for input file: $inputFile"
            exit 1
        fi

        Lx=2.0e-03; Ly=2.0e-03; Lz=2.00e-03
        Nx=1100; Ny=1100; Nz=1100
        eps=9.09629658751972e-07
    elif [[ $inputFile == *"grainReadFile-96_s1-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        # Check that dim = 2
        if [[ $dim -ne 2 ]]; then
            echo "Error: Dimension mismatch. Expected dim = 2 for input file: $inputFile"
            exit 1
        fi

        Lx=2.03e-03; Ly=2.03e-03; Lz=2.03e-03
        Nx=1100; Ny=1100; Nz=1100
        eps=9.09629658751972e-07    
    elif [[ $inputFile == *"grainReadFile-30_s1-10.dat"* ]]; then
        echo "Using file $inputFile"
        echo " "

        # Check that dim = 2
        if [[ $dim -ne 2 ]]; then
            echo "Error: Dimension mismatch. Expected dim = 2 for input file: $inputFile"
            exit 1
        fi

        Lx=2.0e-03; Ly=2.0e-03; Lz=2.0e-03
        Nx=1100; Ny=1100; Nz=1100
        eps=9.09629658751972e-07

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

# Function to finalize results
finalize_results() {
    echo "Finalizing results..."
    cp NASAv2.c run_NASAv2.sh plotNASA.py plotSSA.py plotPorosity.py $folder
    write_parameters_to_csv

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

################################################################################
########################## Main execution starts here ##########################
################################################################################

echo " "
echo "Starting NASAv2 simulation workflow"
echo " "

delt_t=1.0e-4
t_final=14*24*60*60
n_out=40 #100
t_final=$(echo "$t_final" | bc -l)
humidity=1.00
temp=-80.0
grad_temp0X=0.0
grad_temp0Y=3.0e-4
grad_temp0Z=0.0
dim=2
filename="grainReadFile-30_s1-10.dat"
title="NASAv2_30G-${dim}D_T${temp}_hum${humidity}_"

compile_code
create_folder
set_parameters
finalize_results
run_simulation
run_plotting

echo "-------------------------------------------------------------------------"
echo "Done!"
echo " "
