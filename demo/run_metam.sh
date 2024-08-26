#!/bin/bash

#SBATCH -J metam
#SBATCH -t 36:00:00  
##SBATCH --nodes=1
#SBATCH --ntasks=32
##SBATCH --exclusive
#SBATCH -o "%x.o%j"
#SBATCH -e "%x.e%j"
#SBATCH --export=ALL
##SBATCH --mem-per-cpu=1G   # memory per CPU core
##SBATCH -p                 #general partition

### /SBATCH -o slurm.%N.%j.out # STDOUT
### /SBATCH -e slurm.%N.%j.err # STDERR

 echo " "
 echo "compiling Metamorph"
 echo " "
 make Metamorph

echo ------------------------------------------------------
echo 'Job is running on node '; srun hostname
echo ------------------------------------------------------
echo  qsub is running on $SLURM_SUBMIT_HOST
echo  executing queue is $SLURM_JOB_ACCOUNT
echo  working directory is $SLURM_SUBMIT_DIR
echo  partition/queue is $SLURM_JOB_PARTITION
echo  job identifier is $SLURM_JOB_ID
echo  job name is $SLURM_JOB_NAME
echo  node file is $SLURM_JOB_NODELIST
echo  cluster  $SLURM_CLUSTER_NAME
echo  total nodes $SLURM_JOB_NUM_NODES
echo ------------------------------------------------------

echo " "
echo "setting up things"
echo " "

 id=${SLURM_JOB_ID:0:9}
 echo $id

 nam=mm_$id
 folder=/central/scratch/amoure/$nam
 echo $nam
 echo $folder
 export folder

 angl=0         # contact angle: 0->real  1->120degrees
 aa=0		# 1->alpha_j=0
 mm=0		# 1->constant mobility
 echo $angl
 echo $aa
 echo $mm
 export angl
 export aa
 export mm

 export I_MPI_PMI_LIBRARY=/path/to/slurm/pmi/library/libpmi.so
 
 mkdir $folder
 scp ~/PetIGA/demo/Metamorph.c $folder/
 scp ~/PetIGA/demo/run_metam.sh $folder/

 cd $SLURM_SUBMIT_DIR
 
echo " "
echo "running Metamorph"
echo " "

mpirun ./Metamorph -initial_cond "sol3035.dat" -initial_PFgeom -snes_rtol 1e-3 -snes_stol 1e-6 -snes_max_it 5 -ksp_gmres_restart 150 -ksp_max_it 500 -ksp_converged_maxits 1 -ksp_converged_reason -snes_converged_reason  -snes_linesearch_monitor -snes_linesearch_type basic | tee $folder/outp.txt

##sleep 1000
echo "done"
