#!/bin/bash

#SBATCH --job-name=run_prod
#SBATCH --partition=gpu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-task=1
#SBATCH --output=run_prod.out
#SBATCH --error=run_prod.error
#SBATCH --time=48:00:00

echo "Loading modules..."
module load amber  
source /opt/calstatela/amber-22/amber22/amber.sh


mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod1.out \
-p 6m0j_noHet_solvated.prmtop -c equil.rst -r prod1.rst -x prod1.mdcrd

echo "prod 1 done ..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod2.out \
-p 6m0j_noHet_solvated.prmtop -c prod1.rst -r prod2.rst -x prod2.mdcrd

echo "prod 2 done..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod3.out \
-p 06m0j_noHet_solvated.prmtop -c prod2.rst -r prod3.rst -x prod3.mdcrd

echo "prod 3 done..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod4.out \
-p 6m0j_noHet_solvated.prmtop -c prod3.rst -r prod4.rst -x prod4.mdcrd


