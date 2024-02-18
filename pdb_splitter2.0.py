import os
import sys
import subprocess

def input_args_check( input_arg_path) -> dict :
    cwd = os.getcwd()


    #TODO add new inputs 
    '''
    want option for default names, and dynamic naming: 
        -require dynamic naming if not provided all in paths

    wnat to have inputs for the leap.in fatures
    want mmbpsa input opitons (include step stuff) :
        inputs: 
        start_frame :
        end_frame
        interval : 

    '''
    input_fields={"WILD_TYPE": [], 
                "MUTATIONS":[],
                "*MDCRD_DIRECTORY": cwd, 
                "LEAP.IN_PATH" : [], 
                "MMPBSA.IN_PATH": [],
                "MMPBSA.SH_PATH": []
                }

    with open("tmp_input_file.txt", "r") as input_file:
        for line in input_file:
            if line.startswith("#input"):
                tmp_key = line.split() #split
                if tmp_key[2:] != [] : #check for null input 
                    input_fields[tmp_key[1]] = tmp_key[2:] #for args of len >1 e.g. MUTATIONS
    
    
    ##TODO error without an input pdb and a mutation 
    return input_fields

def split_and_mut(pdbfh, pdbfh_base_name, name_from, name_to, idx, naming_conv) :
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
    #mutations:
    #ligand_mutation
    mutation_pdb_data = mutations(ligand_pdb, name_from, name_to, idx)
    file_handle_mut_base = pdbfh_base_name +"_" + naming_conv
    file_handle_mut = file_handle_mut_base + "_ligand.pdb"
    with open(file_handle_mut, "w+") as pdb_file : 
        for line in mutation_pdb_data : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    #full file mutation
    mutation_pdb_data_all = mutations(pdb_data, name_from, name_to, idx) 
    file_handle_mut_all = file_handle_mut_base + ".pdb"
    with open(file_handle_mut_all, "w+") as pdb_file : 
        for line in mutation_pdb_data_all : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    return  file_handle_mut_base

def tleap_in_gen( input_dict,  pdbfh_base_name, file_handle_mut_base ): 
    if input_dict["LEAP.IN_PATH" ] == [] :
        tleap_mut_in = tleap_gen(pdbfh_base_name, file_handle_mut_base)
        with open("tleap_mut.in", "w+") as tleap : 
            for line in tleap_mut_in : 
                tleap.write(f"{line}\n")
            tleap.close()
        os.system(f"dos2unix tleap_mut.in") #not sure if needed. 
        tleap_file_name ="tleap_mut.in"
    else : 
        tleap_file_name =input_dict["LEAP.IN_PATH" ]
    return  tleap_file_name


def mmbpsa_sh_gen(input_dict , pdbfh_base_name , file_handle_mut_base, cwd):
    if input_dict["MMPBSA.SH_PATH"] == []:
        mut_bash_file = mut_bash(pdbfh_base_name, file_handle_mut_base, cwd)
        
        with open("run_MMPBSA.sh", "w+") as mut_bash_sh : 
            for line in mut_bash_file : 
                mut_bash_sh.write(f"{line}\n")
            mut_bash_sh.close()
            run_MMPBSA_sh_name = "run_MMPBSA.sh"
    else :
        run_MMPBSA_sh_name = input_dict["MMPBSA.SH_PATH"]
    return run_MMPBSA_sh_name

def mmbpsa_in_gen(input_dict ) : 
    if input_dict["MMPBSA.IN_PATH"] == []:
        mmpbsa_in_list = mmpbsa_in()
        with open("mmpbsa.in", "w+") as mmpbsa:
            for line in mmpbsa_in_list :
                mmpbsa.write(f"{line}")
            mmpbsa.close()

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

def mutations(pdb_data, name_from, name_to, idx) -> list:
    '''
    pdb_data    : pdb file we want to mutate
    name_from   : three letter name for the initial amino acid 
    name_to     : three letter name for the final amino acid
    idx         : index of the mutation 
    
    returns     : mutated pdb as a list 
    '''
    counter = 0 
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    mutated_ligand = [] 
    
    for line in pdb_data:
        if line.startswith(records):
            
            
            if (line[17:20].strip() in name_from) and (line[22:26].strip() == idx):
                #add cases here for longer ones... 
                #add case for GLN which has AGLN and BGLN
                if (counter <= 4) and (line[16] == "A" or line[16]== " ") : #count for ALA:  
                    new_line= line[:16] + name_to.rjust(4) + line[20:]
                    mutated_ligand.append(new_line)
                    counter = counter + 1 
                    continue
                else :
                    continue
        mutated_ligand.append(line)

    return mutated_ligand


def tleap_gen(pdbfh_base_name,file_handle_mut_all ) -> list:
    '''
    pdbfh_base_name     : base name of the pdb file
    file_handle_mut_all : base bame of the mutated file 
    returns             : tleap file as a list 
    
    #TODO: add inputs for radii, leaprc, TIP3PBOX
    
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
        f"com_mut = loadpdb {file_handle_mut_all}.pdb",
        f"ligand_mut = loadpdb {file_handle_mut_all}_ligand.pdb\n",
        f"saveamberparm com_mut {file_handle_mut_all}.prmtop {file_handle_mut_all}.inpcrd",
        f"saveamberparm ligand_mut {file_handle_mut_all}_ligand.prmtop {file_handle_mut_all}_ligand.inpcrd",
        f"quit"]
    return tleap_wild_in

def mut_bash( pdbfh_base_name, file_handle_mut_all, cwd) :
    '''
    pdbfh_base_name     : base name of the pdb file
    file_handle_mut_all : base bame of the mutated file 
    cwd                 : cwd where script was called, this is parent dir afer chdir call
    returns             : .sh file as a list 

    TODO: consider using input for nodes, making use of queue on nodes to not take up all of the cluster
    TODO change mmpbsa.in into path to mmbpsa file : 
    '''
    
    mut_bash_sh = [f"#!/bin/bash",
        f"#SBATCH --job-name=run_66_mut",
        f"#SBATCH --partition=cpu",
        f"#SBATCH --ntasks=4",
        f"#SBATCH --cpus-per-task=1",
        f"#SBATCH --mem=8000"
        f"#SBATCH --output=run_mmpbsa_66.out",
        f"#SBATCH --error=run_mmpbsa_66.error",
        f"#SBATCH --time=72:00:00",
        f'''echo "Loading modules..."'''  ,  
        f"module load amber " ,
        f"source /opt/calstatela/amber-22/amber22/amber.sh",
        f"",
        f"tleap -s -f tleap_mut.in > tleap_mut.out"
        f"",
        f"""mpirun -np 4 $AMBERHOME/bin/MMPBSA.py -O -i \
{cwd}/mmpbsa.in -o \
FINAL_RESULTS_MMPBSA_tleap_{file_handle_mut_all}.dat\
 -sp {pdbfh_base_name}_solvated.prmtop\
 -cp {pdbfh_base_name}.prmtop\
 -rp {pdbfh_base_name}_recpt.prmtop\
 -lp {pdbfh_base_name}_ligand.prmtop\
 -y {cwd}/*.mdcrd\
 -mc {file_handle_mut_all}.prmtop\
 -ml {file_handle_mut_all}_ligand.prmtop"""
        ]
    return mut_bash_sh

def mmpbsa_in()->list:
    """
    #TODO make this custom, 
    inputs: 
    start_frame :
    end_frame
    interval : 
    
    gb:
    igb:
    saltcon
    
    pb:
    istring:
    
    return: mmpbsa_in_data as a list
    
    # use line.startswith(records): for case structure
    """
    """
    removed this for testing 
    /
&pb
  istrng=0.100
    
    #changed startfram =1; and interval=2 ; removed end frame 
    keep_files=0; don't want tmp files. 
    
    """
    mmpbsa_in_data = [
"""
sample input file for running alanine scanning
 &general
   startframe=1, interval=2,
   verbose=1, keep_files=0
/
&gb
  igb=66, saltcon=0.1

/
&alanine_scanning
/
"""]
    return mmpbsa_in_data

##TODO remove this 
def input_args_from_text( file_handle ) -> list :
    """
    #TODO: will want to be able to read inputs based on flags from input file
    
    can do option for file path to mmpbsa
    can do option for file path to tleap input file
    #input for *mdcrd input file generation 
    
    
    for full list of mmpbsa args: https://ambermd.org/doc12/Amber22.pdf
    """
    
    with open(file_handle) as input_file :
        for line in input_file : 
            input_arg_list = [] #add real input here
        
    return input_arg_list 
 
def general_method(input_dict, pdbfh, pdbfh_base_name, mutation) : 
    """
    general process: 
    input_dict: from input_args_check
    mut_num : range(len(input_dict["MUTATIONS"]))
    """
    
    name_from_char, idx, name_to_char = mutation.split(":")
    #dict for char to code conversion 
    amino_acid_dict = {"A": "ALA", "V": "VAL", \
                    "I":"ILE", "L":"LEU", "M":"MET", \
                    "F":"PHE", "Y":"TYR", "W":"TRP" , \
                    "S":"SER", "T":"THR", "N":"ASN", \
                    "Q":"GLN", "C":"CYS","U":"SEC", \
                    "G":"GLY", "P":"PRO", "R":"ARG", \
                    "H":"HIS", "K":"LYS", "D":"ASP", \
                    "E":"GLU"}
    #convert to 3 letter code
    name_from   = amino_acid_dict[name_from_char]
    name_to     = amino_acid_dict[name_to_char]
    
    #get cwd
    cwd = os.getcwd()
    #get the three letter code as a str e.g. E484A
    naming_conv = name_from_char+ idx + name_to_char
    #make a directory named: base_name_naming-conv_dir
    #string for dir name
    dir_name = pdbfh_base_name + "_" + naming_conv + "_dir"
    #path to dir
    dir_name_path = "./" + dir_name
    #path for other files
    os.system(f"mkdir {dir_name_path}")
    
    dir_name_path_full = dir_name_path + "/"
    #new name for pdb in the dir
    pdbfh_in_dir = dir_name_path_full + pdbfh
    
    #copy the base pdb into new dir
    os.system(f"cp {pdbfh} {pdbfh_in_dir}")
    #update base name to the file in subdir
    #pdbfh_base_name = pdbfh_base_name_in_dir
    os.chdir(dir_name_path) #note the change back use chdir("..")
    
    ##TODO make into smaller functions 
    
    #####################################################################################
    ############################## splitting and mutations ##############################
    #####################################################################################
    file_handle_mut_base = split_and_mut(pdbfh,
                                         pdbfh_base_name, 
                                         name_from, 
                                         name_to, 
                                         idx,
                                         naming_conv)
    #####################################################################################
    #################################### tleap gen ######################################
    #####################################################################################
    ##TODO add options for default
    tleap_file_name = tleap_in_gen( input_dict,
                                    pdbfh_base_name, 
                                    file_handle_mut_base,)
    
    #####################################################################################
    ############################### MMPBSA sh file gen ##################################
    #####################################################################################
    ##TODO add options for default 
    run_MMPBSA_sh_name = mmbpsa_sh_gen(input_dict , 
                                       pdbfh_base_name , 
                                       file_handle_mut_base,
                                       cwd)

    #####################################################################################
    ############################## Running tleap + MMBPSA ###############################
    #####################################################################################
    #run tleap to get solvated files, TODO: put into sbatch
    #os.system(f"tleap -s -f {tleap_file_name} > tleap_mut.out")
    ##TODO: finish the MMPBSA call, just need to have a source for the intial 
    #MMPBSA files, then the .sh should be fine. 
    #os.system(f"sbatch {run_MMPBSA_sh_name}")
    #return to parent dir. 
    os.chdir("..")
    return 


def main():
    
    inputs = sys.argv[1:]
    in_file = inputs[0]
    input_dict = input_args_check(in_file)
    
    pdbfh = input_dict["WILD_TYPE"][0]
    pdbfh_base_name = os.path.basename(pdbfh).split(".")[0] #getting the base name 

    #####################################################################################
    ################################ MMPBSA in file gen #################################
    #####################################################################################
    mmbpsa_in_gen(input_dict )

    ##TODO if mutations iterable then allow for iterable MMPBSA.sh file, user can generate them
    for i in range(len(input_dict["MUTATIONS"])) : 
        mutation = input_dict["MUTATIONS"][i]
        general_method(input_dict, pdbfh, pdbfh_base_name, mutation)
    
    
    ##TODO add creations of summary file: 
    
if __name__ == '__main__':
    main()