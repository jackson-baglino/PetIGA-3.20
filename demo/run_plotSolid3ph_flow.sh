#!/bin/zsh

if [[ -n $1 ]]; then
    echo "Copying run_plotSolid3ph_flow.py to: $1"
    dir=/Users/jacksonbaglino/SimulationResults/junning/$1
    cp plotSolid3ph_flow.py $dir
    # cp writeNASA2CSV.py $dir

    exec_dir=/Users/jacksonbaglino/SimulationResults/junning/$1/

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

    python3 ~/SimulationResults/junning/$1/plotSolid3ph_flow.py
    # python3 ~/SimulationResults/junning/$1/writeNASA2CSV.py
else
    echo "No inputs are given. Assume we are already in the results folder"
    echo " "
    # python3 plotNASA.py
    # python3 writeNASA2CSV.py
fi

