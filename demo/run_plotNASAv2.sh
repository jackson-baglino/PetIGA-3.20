#!/bin/zsh

if [[ -n $1 ]]; then
    echo "Copying run_plotNASAv2.py to: $1"
    dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/$1
    echo " "
    echo " "
    echo " --- Current directory is: $(pwd)"
    echo " "
    echo " "
    cp writeNASA2CSV.py $dir
 
    exec_dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/$1/

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

    python3.11 ~/SimulationResults/DrySed_Metamorphism/NASAv2/$1/plotNASA.py
    # python3.11 ~/SimulationResults/DrySed_Metamorphism/NASAv2/$1/writeNASA2CSV.py

    echo "Plotting SSA and Porosity"
    echo " "
    echo "Calling from: $1"

    python3.11 ~/SimulationResults/DrySed_Metamorphism/NASAv2/$1/plotSSA.py
    # python3.11 ~/SimulationResults/DrySed_Metamorphism/NASAv2/$1/plotPorosity.py


    cp ~/PetIGA-3.20/demo/rename_STL_files.py $dir
else
    echo "No inputs are given. Assume we are already in the results folder"
    echo " "

    echo "Creating vtkOut and stlOut directories"
    echo " "

    python3.11 plotNASA.py
    
    python3.11 plotSSA.py
fi