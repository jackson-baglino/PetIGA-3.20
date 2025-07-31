#!/bin/zsh

# Function to change to the target directory and print status
change_to_directory() {
  echo " "
  echo "Changing to directory: $folder"
  cd $folder || { echo "Failed to change directory to $folder"; exit 1; }
  echo "Directory changed to: $(pwd)"
}

# Function to execute Python scripts for plotting
execute_python_scripts() {
  local dir=$folder
  echo "Executing plotDSM.py in directory: $folder"
  python $folder/plotDSM.py

  echo "Plotting SSA and Porosity"
  python $folder/plotSSA.py
  python $folder/plotPorosity.py
}

# Main script logic
echo "Running plotting scripts in the current directory."

echo "Creating vtkOut and stlOut directories"
mkdir -p vtkOut stlOut

# Execute Python scripts in the current directory
change_to_directory

echo " "
echo "Current directory: $(pwd)"
echo " "

echo "Executing plotDSM.py and plotSSA.py in the current directory"

mkdir -p vtkOut

echo "Python version is $(python3 --version 2>&1)"

python plotDSM.py
python plotSSA.py