#!/bin/bash
#SBATCH --job-name=MD_pipeline
#SBATCH --partition=gpu
#SBATCH --exclusive
#SBATCH --output=MD_pipeline.out
#SBATCH --error=MD_pipeline.error
#SBATCH --time=72:00:00

echo "Loading modules..."
module load amber  
source /opt/calstatela/amber-22/amber22/amber.sh

$AMBERHOME/bin/pmemd.cuda -O -i min.in -o min.out -p 6m0j_noHetASN21_solvated.prmtop -c 6m0j_noHetASN21_solvated.inpcrd -r min.rst 
echo "min done"
$AMBERHOME/bin/pmemd.cuda -O -i heat.in -o heat.out -p 6m0j_noHetASN21_solvated.prmtop -c min.rst -r heat.rst -x heat.mdcrd -ref min.rst

echo "heat done"
$AMBERHOME/bin/pmemd.cuda -O -i density.in -o density.out -p 6m0j_noHetASN21_solvated.BOX.prmtop -c heat.rst -r density.rst -x density.mdcrd -ref heat.rst

echo "density done"
$AMBERHOME/bin/pmemd.cuda -O -i equil.in -o equil.out -p 6m0j_noHetASN21_solvated.prmtop -c density.rst -r equil.rst -x equil.mdcrd

mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod1.out -p 6m0j_noHetASN21_solvated.prmtop -c equil.rst -r prod1.rst -x prod1.mdcrd

echo "prod 1 done ..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod2.out -p 6m0j_noHetASN21_solvated.prmtop -c prod1.rst -r prod2.rst -x prod2.mdcrd

echo "prod 2 done..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod3.out -p 6m0j_noHetASN21_solvated.prmtop -c prod2.rst -r prod3.rst -x prod3.mdcrd

echo "prod 3 done..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod4.out -p 6m0j_noHetASN21_solvated.prmtop -c prod3.rst -r prod4.rst -x prod4.mdcrd
