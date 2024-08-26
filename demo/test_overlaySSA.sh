#!/bin/zsh

echo " "
echo "Testing 'plotSSA.py' script"
echo " "

# Set the necessary environment variables
dir="./input/"

dim=3
inputFile=$dir"grainReadFile-2.dat"
outputfolder=$dir"test/"

# Export variables
export dim inputFile

# Call the script
python3.11 plotSSA.py