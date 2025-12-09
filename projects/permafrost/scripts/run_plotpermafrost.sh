#!/bin/zsh

# Function to copy files to the target directory
copy_files_to_directory() {
    local target_dir=$1
    echo "Copying necessary files to: $target_dir"
    cp ./scripts/writepermafrost2CSV.py $target_dir
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
    echo "Executing plotpermafrost.py in directory: $dir"
    python $dir/plotpermafrost.py
    
    # echo "Plotting SSA and Porosity"
    # python $dir/plotSSA.py
    # Uncomment the line below if you want to include Porosity plotting
    # python $dir/plotPorosity.py
}

# Main script logic
if [[ -n $1 ]]; then
    echo "Starting process for directory: $1"
    dir=/Users/jacksonbaglino/SimulationResults/permafrost/scratch/$1

    # Copy files and switch to the target directory
    copy_files_to_directory $dir
    change_to_directory $dir

    echo "Creating vtkOut directory"
    mkdir -p vtkOut

    # Execute plotting scripts
    execute_python_scripts $dir
else
    echo "No inputs provided. Assuming we are already in the results folder."
    
    echo "Creating vtkOut and directory"
    mkdir -p vtkOut

    # Execute Python scripts in the current directory
    echo "Executing plotpermafrost.py and plotSSA.py in the current directory"
    python plotpermafrost.py
    # python plotSSA.py
    # Uncomment if you need to execute plotPorosity.py in the current directory
    # python plotPorosity.py
fi