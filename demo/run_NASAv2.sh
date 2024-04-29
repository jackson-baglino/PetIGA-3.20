#!/bin/zsh

echo " "
echo "compiling"
echo " "
make NASAv2

# add name folder accordingly:    <---------------------------------------------
name=res_$(date +%Y-%m-%d__%H.%M.%S)
dir=/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2
folder=$dir/$name

if [ ! -d "$dir" ]; then
  mkdir -p "$dir"
fi

mkdir $folder/

export folder

cp NASAv2.c$folder
cp run_NASAv2.sh $folder

echo " "
echo "Calling ./NASAv2"
echo " "

mpiexec -np 6 ./NASAv2 -initial_PFgeom -temp_initial -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 -ksp_gmres_restart 150 -ksp_max_it 1000  -ksp_converged_reason -snes_converged_reason -snes_linesearch_monitor -snes_linesearch_type basic | tee /Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/outp.txt

echo " "
echo "making directory" $folder
echo " "

# mv $dir/*.dat $folder
mv $dir/outp.txt $folder
cp NASAv2.c $folder
cp run_plotNASAv2.sh $folder

echo "Queing plotNASA.py"
./run_plotNASAv2.sh $name

echo "-------------------------------------------------------------------------"
echo "Done!"
echo " "

