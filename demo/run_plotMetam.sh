#!/bin/zsh

if [[ -n $1 ]]; then
    echo "Copying run_plotMetam.py to: $1"
    dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/wetamorph/$1
    cp plotMetam.py $dir

    exec_dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/wetamorph/$1/

    echo " "
    echo $exec_dir
    echo " "

    cd $exec_dir
    echo " "
    echo "Directory changed"
    echo " "
    echo "Working directory: $(pwd)"
    echo " "
    echo "Executing python script"
    echo " "

    mkdir $dir/vtkOut

    python3.11 ~/SimulationResults/DrySed_Metamorphism/wetamorph/$1/plotMetam.py

    echo "Plotting SSA and Porosity"
    echo " "
    echo "Calling from: $1"

    python3.11 ~/SimulationResults/DrySed_Metamorphism/wetamorph/$1/plotSSA.py
else
    echo "No inputs are given. Assume we are already in the results folder"
    echo " "

    echo "Creating vtkOut and stlOut directories"
    echo " "
    
    mkdir -p vtkOut
    mkdir -p stlOut

    python3.11 plotMetam.py
    
    python3.11 plotSSA.py

    cp -r -p stlOut stlOut-copy
    python3.11 rename_STL_files.py
fi

