#!/bin/zsh

# Locate plotpermafrost.py: prefer the staged postprocess/ subdirectory,
# fall back to the root of the results folder for older result dirs.
find_plot_script() {
    local dir=$1
    if [[ -f "$dir/postprocess/plotpermafrost.py" ]]; then
        echo "$dir/postprocess/plotpermafrost.py"
    elif [[ -f "$dir/plotpermafrost.py" ]]; then
        echo "$dir/plotpermafrost.py"
    else
        echo ""
    fi
}

# Function to execute Python scripts for plotting
execute_python_scripts() {
    local dir=$1
    local plot_py
    plot_py=$(find_plot_script "$dir")

    if [[ -z "$plot_py" ]]; then
        echo "⚠️  plotpermafrost.py not found in $dir or $dir/postprocess/ — skipping."
        return 1
    fi

    echo "Executing: $plot_py"
    python3 "$plot_py"
}

# Main script logic
if [[ -n $1 ]]; then
    echo "Starting process for directory: $1"
    dir=/Users/jacksonbaglino/SimulationResults/permafrost/scratch/$1

    echo "Creating vtkOut directory"
    mkdir -p "$dir/vtkOut"

    execute_python_scripts "$dir"
else
    echo "No inputs provided. Assuming we are already in the results folder."

    echo "Creating vtkOut directory"
    mkdir -p vtkOut

    plot_py=$(find_plot_script "$(pwd)")
    if [[ -z "$plot_py" ]]; then
        echo "⚠️  plotpermafrost.py not found — skipping."
        exit 1
    fi
    echo "Executing: $plot_py"
    python3 "$plot_py"
fi