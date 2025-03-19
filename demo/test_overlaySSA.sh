#!/bin/zsh

echo " "
echo "Testing 'plotSSA.py' script"
echo " "

# Set the necessary environment variables
dir="/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2"
datFile="SSA_evo.dat"
filename2D=$dir"/2024-08-27__08.46.56/"$datFile
filename3D=$dir"/2024-08-26__15.00.09/"$datFile
inputfile="./input/grainReadFile-2.dat"
output_file="/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/ssa_evolution_plot-2Grain.png"


# Export variables
export dir filename2D filename3D output_file inputfile

# Call the script
python3.11 overlaySSA_plot.py