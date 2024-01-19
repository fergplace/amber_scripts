import os
import sys
#split our file, call with 0 for struc, 1 for cov 

def pdb_split(pdb_data, option) -> list:

    #ignore HET and other line starts: 
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

def mutations(pdb_data, name_from, name_to, idx) -> list:

    counter = 0 
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    mutated_cov = [] 
    
    for line in pdb_data:
        if line.startswith(records):
            if (line[17:20].strip() in name_from) and (line[22:26].strip() in idx):
                #add cases here for longer ones... 
                if counter <= 4 : #count for ALA:  
                    new_line= line[:17] + name_to.rjust(3) + line[20:]
                    mutated_cov.append(new_line)
                    counter = counter + 1 
                    continue
                else :
                    continue
        mutated_cov.append(line)

    return mutated_cov


def tleap_gen(pdbfh_base_name,file_handle_mut_all ):
    #standard leap.in for mut files 
    #TODO add options for radii, box, FF
    tleap_wild_in = [f"source /oldff/leaprc.ff99",
        f"source leaprc.water.tip3p",
        f"set default PBRadii mbondi2",
        f"com = loadpdb {pdbfh_base_name}.pdb"  ,
        f"cov = loadpdb{ pdbfh_base_name}_cov.pdb" ,
        f"rcp = loadpdb {pdbfh_base_name}_recpt.pdb",
        f"saveamberparm com {pdbfh_base_name}.prmtop {pdbfh_base_name}.inpcrd",
        f"saveamberparm cov {pdbfh_base_name}_cov.prmtop {pdbfh_base_name}_cov.inpcrd",
        f"saveamberparm rcp {pdbfh_base_name}_recpt.prmtop {pdbfh_base_name}_recpt.inpcrd",
        f"solvatebox com TIP3PBOX 12.0",
        f"saveamberparm com {pdbfh_base_name}_solvated.prmtop {pdbfh_base_name}_solvated.inpcrd",
        f"com_mut = loadpdb {file_handle_mut_all}.pdb",
        f"cov_mut = loadpdb {file_handle_mut_all}_cov.pdb",
        f"saveamberparm com_mut {file_handle_mut_all}.prmtop {file_handle_mut_all}.inpcrd",
        f"saveamberparm cov_mut {file_handle_mut_all}_cov.prmtop {file_handle_mut_all}_cov.inpcrd",
        f"quit"]
    return tleap_wild_in

def mut_bash( pdbfh_base_name, file_handle_mut_all) :
    
    
    mut_bash_sh = [f"#!/bin/bash",
        f"#SBATCH --job-name=run_66_mut",
        f"#SBATCH --partition=gpu",
        f"#SBATCH --output=run_mmpbsa_66.out",
        f"#SBATCH --error=run_mmpbsa_66.error",
        f"#SBATCH --time=48:00:00",
        f'''echo "Loading modules..."'''  ,  
        f"module load amber " ,
        f"source /opt/calstatela/amber-22/amber22/amber.sh",
        f"""$AMBERHOME/bin/MMPBSA.py.MPI -O -i\
mmpbsa_mut_66.in -o\
FINAL_RESULTS_MMPBSA_tleap_{file_handle_mut_all}.dat\
 -sp {pdbfh_base_name}_solvated.prmtop\
 -cp {pdbfh_base_name}.prmtop\
 -rp {pdbfh_base_name}_recpt.prmtop\
 -lp{pdbfh_base_name}_cov.prmtop \
 -y *.mdcrd\
 -mc {file_handle_mut_all}.prmtop\
 -ml {file_handle_mut_all}_cov.prmtop"""
        ]
    return mut_bash_sh
def main():
    
    inputs = sys.argv[1:]
    pdbfh = inputs[0]
    pdbfh_base_name = pdbfh.split(".")[0]
    name_from_char, idx, name_to_char = inputs[1].split(":")
    #dict for char to code conversion 
    amino_acid_dict = {"A": "ALA", "V": "VAL", \
                    "I":"LLE", "L":"LEU", "M":"MET", \
                    "F":"PHE", "Y":"TYR", "W":"TRP" , \
                    "S":"SER", "T":"THR", "N":"ASN", \
                    "Q":"GLN", "C":"CYS","U":"SEC", \
                    "G":"GLY", "P":"PRO", "R":"ARG", \
                    "H":"HIS", "K":"LYS", "D":"ASP", \
                    "E":"GLU"}
    #convert to 3 letter code
    name_from   = amino_acid_dict[name_from_char]
    name_to     = amino_acid_dict[name_to_char]
    
    naming_conv = name_from_char+ idx + name_to_char
    #open the file:
    with open(pdbfh, "r") as f :
        pdb_data = f.readlines()

    #splits
    struct_pdb_data = pdb_split(pdb_data, 1 )
    file_handle_structure = pdbfh_base_name + "_struct.pdb"
    with open(file_handle_structure, "w+") as pdb_file : 
        for line in struct_pdb_data : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    
    cov_pdb = pdb_split(pdb_data, 0 )
    file_handle_covid = pdbfh_base_name + "_cov.pdb"
    with open(file_handle_covid, "w+") as pdb_file : 
        for line in cov_pdb : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    
    #mutations:
    #cov_mutation
    mutation_pdb_data = mutations(cov_pdb, name_from, name_to, idx)
    file_handle_mut_base = pdbfh_base_name +"_" + naming_conv + ".pdb"
    file_handle_mut = file_handle_mut_base + "_cov.pdb"
    with open(file_handle_mut, "w+") as pdb_file : 
        for line in mutation_pdb_data : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    #full file mutation
    mutation_pdb_data_all = mutations(pdb_data, name_from, name_to, idx)
    file_handle_mut_all_base = pdbfh_base_name +"_" + naming_conv 
    file_handle_mut_all = file_handle_mut_all_base + ".pdb"
    with open(file_handle_mut_all, "w+") as pdb_file : 
        for line in mutation_pdb_data_all : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    
    
    #tleap gen
    tleap_mut_in = tleap_gen(pdbfh_base_name, file_handle_mut_all_base)
    with open("tleap_mut.in", "w+") as tleap : 
        for line in tleap_mut_in : 
            tleap.write(f"{line}\n")
        tleap.close()

    mut_bash_file = mut_bash(pdbfh_base_name, file_handle_mut_all_base)
    with open("run_MMPBSA.sh", "w+") as mut_bash_sh : 
        for line in mut_bash_file : 
            mut_bash_sh.write(f"{line}\n")
        mut_bash_sh.close()
        
if __name__ == '__main__':
    main()