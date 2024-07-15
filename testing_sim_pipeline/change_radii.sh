#!/bin/bash
#SBATCH --job-name=run_66_mut
#SBATCH --partition=cpu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --output=run_mmpbsa_66.out
#SBATCH --error=run_mmpbsa_66.error
#SBATCH --time=72:00:00
echo "Loading modules..."
module load amber 
source /opt/calstatela/amber-22/amber22/amber.sh

tleap -s -f tleap_mut.in > tleap_mut.out
c:\Users\13108\Documents\GitHub\amber_scripts\testing_sim_pipeline/change_radii_to_opt.py 6m0j_noHet_solvated.prmtop
c:\Users\13108\Documents\GitHub\amber_scripts\testing_sim_pipeline/change_radii_to_opt.py 6m0j_noHet.prmtop
c:\Users\13108\Documents\GitHub\amber_scripts\testing_sim_pipeline/change_radii_to_opt.py 6m0j_noHet_ligand.prmtop
c:\Users\13108\Documents\GitHub\amber_scripts\testing_sim_pipeline/change_radii_to_opt.py 6m0j_noHet_recpt.prmtop
