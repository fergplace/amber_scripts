#!/bin/bash

#SBATCH --job-name=run_66_mut
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --output=run_mmpbsa_66.out
#SBATCH --error=run_mmpbsa_66.error
#SBATCH --time=48:00:00

echo "Loading modules..."
module load amber  
source /opt/calstatela/amber-22/amber22/amber.sh

$AMBERHOME/bin/MMPBSA.py -O -i mmpbsa_mut_66.in -o FINAL_RESULTS_MMPBSA_tleap_COV_MUT_66.dat -sp 6m0j_noHet_solvated.prmtop -cp 6m0j_noHet.prmtop -rp 6m0j_noHet_recpt.prmtop -lp 6m0j_noHet_cov.prmtop  -y *.mdcrd -mc 6m0j_E484A.prmtop -ml 6m0j_E484A_cov.prmtop

