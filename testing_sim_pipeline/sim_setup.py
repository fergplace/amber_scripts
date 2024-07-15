import os
import sys
import subprocess
import numpy as np 

def pdb_split(pdb_data, option) -> list:
    '''
    pdb_data    : pdb file we want to split
    option      : option will determine if we want to get the receptor or the ligand,
    receptor =0, ligand =1 
     
    returns     : split as a list 
    '''
    #ignore HET and other line starts: 
    #NOTE: no need to use the no_HET source files, this will strip the files of the HET
    ter_state = 0 
    records = ('ATOM', 'ANISOU', 'TER')
    data = []
    
    for line in pdb_data:
        if line.startswith(records):
            if (option == 0) and (ter_state==0)  : 
                data.append(line)
                
                if line.startswith('TER') :
                    return data #break once we get to first Ter as option 0
            
            #need to check for Ter after store line starting with Ter due to structure 
            #of pdb files, TER line belongs to structure. 
            if line.startswith('TER') and (ter_state==0) :
                ter_state =1
                continue 
            if (option == 1 ) and (ter_state==1): 
                data.append(line)
                
    return data


def tleap_in_gen( pdbfh_base_name ): 
   
	tleap_mut_in = tleap_gen(pdbfh_base_name )
	with open("tleap_mut.in", "w+") as tleap : 
		for line in tleap_mut_in : 
			tleap.write(f"{line}\n")
		tleap.close()
	os.system(f"dos2unix tleap_mut.in") #not sure if needed. 
	tleap_file_name ="tleap_mut.in"
    
	return  tleap_file_name

def split_and_mut(pdbfh, pdbfh_base_name) :
    with open(pdbfh, "r") as f :
        pdb_data = f.readlines()
    #splits
    struct_pdb_data = pdb_split(pdb_data, 0 )
    file_handle_structure = pdbfh_base_name + "_recpt.pdb"
    with open(file_handle_structure, "w+") as pdb_file : 
        for line in struct_pdb_data : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    ligand_pdb = pdb_split(pdb_data, 1 )
    file_handle_ligand = pdbfh_base_name + "_ligand.pdb"
    with open(file_handle_ligand, "w+") as pdb_file : 
        for line in ligand_pdb : 
            pdb_file.write(f"{line}")
        pdb_file.close()
        

def tleap_gen(pdbfh_base_name ) -> list:
    '''
    pdbfh_base_name     : base name of the pdb file
    returns             : tleap file as a list 
 
    
    '''
    #standard leap.in for mut files 
    #TODO add options for radii, box, FF
    #TODO try modern FF, 19SB ; 
    #OPC and FF19SB protein.leaprc.ff19SB
    #leaprc.protein.ff19SB  
    #f"source oldff/leaprc.ff99",
    # f"source leaprc.water.tip3p"
    #source leaprc.protein.ff19SB",
       # f"source leaprc.water.opc",
    tleap_wild_in = [f"source leaprc.protein.ff19SB",
        f"source leaprc.water.opc",
        f"set default PBRadii mbondi2\n",
        f"com = loadpdb {pdbfh_base_name}.pdb"  ,
        f"ligand = loadpdb {pdbfh_base_name}_ligand.pdb" ,
        f"rcp = loadpdb {pdbfh_base_name}_recpt.pdb\n",
        f"saveamberparm com {pdbfh_base_name}.prmtop {pdbfh_base_name}.inpcrd",
        f"saveamberparm ligand {pdbfh_base_name}_ligand.prmtop {pdbfh_base_name}_ligand.inpcrd",
        f"saveamberparm rcp {pdbfh_base_name}_recpt.prmtop {pdbfh_base_name}_recpt.inpcrd",
        f"solvatebox com TIP3PBOX 12.0",
        f"saveamberparm com {pdbfh_base_name}_solvated.prmtop {pdbfh_base_name}_solvated.inpcrd\n",
        f"quit"]
    return tleap_wild_in

def change_radii_sh( pdbfh_base_name, cwd  ) :
    '''
	bash to change radii 
    '''
    
    radii_sh = [f"#!/bin/bash",
        f"#SBATCH --job-name=run_66_mut",
        f"#SBATCH --partition=cpu",
        f"#SBATCH --ntasks=4",
        f"#SBATCH --cpus-per-task=1",
        f"#SBATCH --mem=10000",
        f"#SBATCH --output=run_mmpbsa_66.out",
        f"#SBATCH --error=run_mmpbsa_66.error",
        f"#SBATCH --time=72:00:00",
        f'''echo "Loading modules..."'''  ,  
        f"module load amber " ,
        f"source /opt/calstatela/amber-22/amber22/amber.sh",
        f"",
        f"tleap -s -f tleap_mut.in > tleap_mut.out"
        f"",
        f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}_solvated.prmtop",
        f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}.prmtop",
        f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}_ligand.prmtop",
        f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}_recpt.prmtop"]

    
    return radii_sh


def main():
    
    pdbfh ="6m0j_noHet.pdb"
    pdbfh_base_name = os.path.basename(pdbfh).split(".")[0]
    ############################## splitting and mutations ##############################
    split_and_mut(pdbfh, pdbfh_base_name)
    
    
    #################################### tleap gen ######################################
    tleap_file_name = tleap_in_gen(pdbfh_base_name)

    
    cwd =  os.getcwd()
    mut_bash_file = change_radii_sh( pdbfh_base_name, cwd) 
            
    with open("change_radii.sh", "w+") as mut_bash_sh : 
        for line in mut_bash_file : 
            mut_bash_sh.write(f"{line}\n")
        mut_bash_sh.close()
        run_MMPBSA_sh_name = "change_radii.sh"
        
if __name__ == '__main__':
    main()