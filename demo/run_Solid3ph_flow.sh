#!/bin/zsh

echo " "
echo "compiling"
echo " "
make Solid3ph_flow
 
# add name folder accordingly:    <---------------------------------------------
name=res_$(date +%Y-%m-%d__%H.%M.%S)
dir=/Users/jacksonbaglino/SimulationResults/junning
folder=$dir/$name

mkdir $folder/

export folder

cp NASA_panda.c  $folder
cp run_NASA_panda.sh $folder

echo " "
echo "Calling ./Solid3ph_flow"
echo " "

mpiexec -np 4 ./Solid3ph_flow  -initial_PFgeom \
-temp_initial -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 7 \
-ksp_gmres_restart 150 -ksp_max_it 1000  -ksp_converged_reason \
-snes_converged_reason  -snes_linesearch_monitor \
-snes_linesearch_type basic | tee /Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASA_panda/outp.txt


echo " "
echo "making directory" $folder
echo " "

# mv $dir/*.dat $folder
# mv $dir/outp.txt $folder

echo "Queing plotNASA.py"
./run_plotSolid3ph_flow.sh $name

echo "-------------------------------------------------------------------------"
echo "Done!"
echo " "

