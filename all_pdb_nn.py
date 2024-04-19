import glob
import numpy as np 
import pandas as pd
import os 

def main()   :
    print("working")
    tmp_file_loc = "*_dir/6m0j_noHet_*.pdb"
    tmp_file_names = glob.glob(tmp_file_loc)
    dir_path = "all_pdb_nn"
    os.system(f"mkdir {dir_path}")
    
    for file in tmp_file_names :
        cp_path = dir_path + "/" + file
        os.system(f"cp {file} {cp_path}")

if __name__ == '__main__':
    main()