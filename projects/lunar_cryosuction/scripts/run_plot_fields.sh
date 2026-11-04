#!/bin/zsh

# Pick an interpreter that can actually import igakit -- plot_fields.py needs it
# to read the PetIGA .dat files, and a bare `python3` usually cannot. Silently
# failing here leaves the run with no vtkOut/, which every downstream
# measurement (neck_width.py, plot_neck_convergence.py) then finds empty.
find_python() {
    local root="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"
    local c
    for c in "$VIRTUAL_ENV/bin/python3" \
             "$root"/venv_*/bin/python3 \
             "$root"/../enceladus_DSM/venv_*/bin/python3 \
             "$(command -v python3)"; do
        [[ -x "$c" ]] || continue
        if "$c" -c 'import igakit' >/dev/null 2>&1; then echo "$c"; return 0; fi
    done
    return 1
}

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

    local PY
    if ! PY=$(find_python); then
        echo "⚠️  no interpreter with igakit found — cannot write vtkOut/."
        echo "    Install it into the project venv:"
        echo "      ./venv_lunar/bin/pip install \\"
        echo "        'igakit @ https://github.com/dalcinl/igakit/archive/refs/heads/master.zip'"
        return 1
    fi
    echo "Executing: $PY $plot_py --dir $dir"
    "$PY" "$plot_py" --dir "$dir"
}

# Main script logic
if [[ -n $1 ]]; then
    # Accept a full absolute path or a bare run-folder name (legacy).
    if [[ "$1" == /* ]]; then
        dir="$1"
    else
        dir=/Users/jacksonbaglino/SimulationResults/enceladus_DSM/scratch/$1
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
    local PY
    if ! PY=$(find_python); then
        echo "⚠️  no interpreter with igakit found — cannot write vtkOut/."
        echo "    Install it into the project venv:"
        echo "      ./venv_lunar/bin/pip install \\"
        echo "        'igakit @ https://github.com/dalcinl/igakit/archive/refs/heads/master.zip'"
        return 1
    fi
    echo "Executing: $PY $plot_py --dir $dir"
    "$PY" "$plot_py" --dir "$dir"
fi