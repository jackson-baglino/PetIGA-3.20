#!/bin/zsh

# Locate plot_fields.py: prefer the staged postprocess/ subdirectory,
# fall back to the root of the results folder for older result dirs.
find_plot_script() {
    local dir=$1
    if [[ -f "$dir/postprocess/plot_fields.py" ]]; then
        echo "$dir/postprocess/plot_fields.py"
    elif [[ -f "$dir/plot_fields.py" ]]; then
        echo "$dir/plot_fields.py"
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
        echo "⚠️  plot_fields.py not found in $dir or $dir/postprocess/ — skipping."
        return 1
    fi

    echo "Executing: $plot_py --dir $dir"
    python3 "$plot_py" --dir "$dir"
}

# Main script logic
if [[ -n $1 ]]; then
    # Accept a full absolute path or a bare run-folder name (legacy).
    if [[ "$1" == /* ]]; then
        dir="$1"
    else
        dir=/Users/jacksonbaglino/SimulationResults/lunar_regolith_DSM/scratch/$1
    fi
    echo "Starting process for directory: $dir"

    echo "Creating vtkOut directory"
    mkdir -p "$dir/vtkOut"

    execute_python_scripts "$dir"
else
    echo "No inputs provided. Assuming we are already in the results folder."

    echo "Creating vtkOut directory"
    mkdir -p vtkOut

    dir="$(pwd)"
    plot_py=$(find_plot_script "$dir")
    if [[ -z "$plot_py" ]]; then
        echo "⚠️  plot_fields.py not found — skipping."
        exit 1
    fi
    echo "Executing: $plot_py --dir $dir"
    python3 "$plot_py" --dir "$dir"
fi