#!/bin/bash

#SBATCH --job-name=run_mmm
#SBATCH --partition=gpu
#SBATCH --output=run_mmpbsa.out
#SBATCH --error=run_mmpbsa.error
#SBATCH --time=48:00:00

echo "Loading modules..."
module load amber  
source /opt/calstatela/amber-22/amber22/amber.sh

$AMBERHOME/bin/MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA_tleap_wild.dat -sp 6m0j_noHet_solvated.prmtop -cp 6m0j_noHet.prmtop -rp 6m0j_noHet_recpt.prmtop -lp6m0j_noHet_cov.prmtop  -y *.mdcrd


