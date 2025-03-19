#!/bin/zsh

# Function to copy files to the target directory
copy_files_to_directory() {
    local target_dir=$1
    echo "Copying necessary files to: $target_dir"
    cp writeNASA2CSV.py $target_dir
    cp ~/PetIGA-3.20/demo/rename_STL_files.py $target_dir
}

# Function to change to the target directory and print status
change_to_directory() {
    local target_dir=$1
    echo " "
    echo "Changing to directory: $target_dir"
    cd $target_dir || { echo "Failed to change directory to $target_dir"; exit 1; }
    echo "Directory changed to: $(pwd)"
}

# Function to execute Python scripts for plotting
execute_python_scripts() {
    local dir=$1
    echo "Executing plotNASA.py in directory: $dir"
    python3.11 $dir/plotNASA.py
    
    echo "Plotting SSA and Porosity"
    python3.11 $dir/plotSSA.py
    # Uncomment the line below if you want to include Porosity plotting
    # python3.11 $dir/plotPorosity.py
}

# Main script logic
if [[ -n $1 ]]; then
    echo "Starting process for directory: $1"
    dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/$1

    # Copy files and switch to the target directory
    copy_files_to_directory $dir
    change_to_directory $dir

    echo "Creating vtkOut and stlOut directories"
    mkdir -p vtkOut stlOut

    # Execute plotting scripts
    execute_python_scripts $dir
else
    echo "No inputs provided. Assuming we are already in the results folder."
    
    echo "Creating vtkOut and stlOut directories"
    mkdir -p vtkOut stlOut

    # Execute Python scripts in the current directory
    echo "Executing plotNASA.py and plotSSA.py in the current directory"
    python3.11 plotNASA.py
    python3.11 plotSSA.py
    # Uncomment if you need to execute plotPorosity.py in the current directory
    # python3.11 plotPorosity.py
fi