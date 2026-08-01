#!/bin/zsh

# Pick an interpreter that can actually import igakit -- plot_fields.py needs it
# to read the PetIGA .dat files, and a bare `python3` usually cannot. Silently
# failing here leaves the run with no vtkOut/, which every downstream
# measurement (neck_width.py, plot_neck_convergence.py) then finds empty.
# Resolve the script's own directory HERE, at file scope, not inside the
# function. In zsh $0 inside a function is the FUNCTION NAME, not the script
# path (FUNCTION_ARGZERO is on by default), so `dirname $0` returned "." and the
# search root became the parent of the caller's cwd -- projects/ rather than
# projects/enceladus_DSM/. BASH_SOURCE does not exist in zsh either, so the
# fallback never helped. ${0:A:h} is zsh for "absolute path, head".
SCRIPT_DIR="${0:A:h}"
DEFAULT_ROOT="${SCRIPT_DIR:h}"      # scripts/ -> project root

find_python() {
    # null_glob: an unmatched pattern expands to nothing instead of aborting the
    # function with "no matches found", which is zsh's default (nomatch) and
    # which killed the search at the FIRST candidate directory that happened not
    # to exist -- before any valid interpreter was ever tried.
    setopt local_options null_glob

    local root="${PROJECT_ROOT:-$DEFAULT_ROOT}"
    local c
    for c in "$VIRTUAL_ENV/bin/python3" \
             "$root"/venv_*/bin/python3 \
             "$root"/../lunar_regolith_DSM/venv_*/bin/python3 \
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
        echo "      ./venv_enceladus/bin/pip install \\"
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
        echo "      ./venv_enceladus/bin/pip install \\"
        echo "        'igakit @ https://github.com/dalcinl/igakit/archive/refs/heads/master.zip'"
        return 1
    fi
    echo "Executing: $PY $plot_py --dir $dir"
    "$PY" "$plot_py" --dir "$dir"
fi