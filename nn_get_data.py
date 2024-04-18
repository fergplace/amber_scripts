import glob
import numpy as np 
import pandas as pd

def main()   :
    print("working")
    tmp_file_loc = "*_dir/FINAL_RESULTS*"
    tmp_file_names = glob.glob(tmp_file_loc)
    mut_cmplx = {}
    
    for file_name in tmp_file_names:
        with open(file_name) as f :
            all_bond = [] 
            all_angle = [] 
            all_dihed = [] 
            all_VDW = [] 
            all_eel = [] 
            all_1_4_VDW = [] 
            all_1_4_eel = [] 
            all_egb = [] 
            all_esurf = [] 
            all_gas = [] 
            all_solv = [] 
            all_e_total = [] 
            data = f.readlines()
            c = 0
            mut = file_name.split("_")[-1].split(".")[0]
            for line in data : 
                #gb-66 trigger 
                if ("GENERALIZED BORN (GBNSR6)" in line)  :
                    c +=1 

                if c >=2 :     
                    if line.startswith("BOND" ) :
                        all_bond.append(line.split()[1])
                    if line.startswith("ANGLE" ) :
                        all_angle.append(line.split()[1])
                    if line.startswith("DIHED" ) :
                        all_dihed.append(line.split()[1])
                    if line.startswith("VDWAALS") :
                        all_VDW.append(line.split()[1])
                    if line.startswith("EEL ") :
                        all_eel.append(line.split()[1])
                        
                    if line.startswith("1-4 VDW")  :
                        all_1_4_VDW.append(line.split()[2])
                        
                    if line.startswith("1-4 EEL") :
                        all_1_4_eel.append(line.split()[2])
                    if line.startswith("EGB")  :
                        all_egb.append(line.split()[1])
                    if line.startswith("ESURF")  :
                        all_esurf.append(line.split()[1])
                    if line.startswith("G gas")  :
                        all_gas.append(line.split()[2])
                    if line.startswith("G solv") :
                        all_solv.append(line.split()[2])
                    if  line.startswith("TOTAL"):
                        all_e_total.append(line.split()[1])
                        
                    if line.startswith("Differences (Complex - Receptor - Ligand):"):
                        break
                    
            data = np.array([all_e_total[0], all_1_4_eel[0], all_eel[0], all_egb[0], all_esurf[0], 
                all_e_total[1], all_1_4_eel[1], all_eel[1], all_egb[1], all_esurf[1],
                all_e_total[2], all_1_4_eel[2], all_eel[2], all_egb[2], all_esurf[2]])
            mut_cmplx[mut] = data
            
            


    cols =np.array( [ "gb-complex-etot","gb-complex-1-4-eel","gb-complex-eelec","gb-complex-egb","gb-complex-esurf",\
            "gb-protein-etot","gb-protein-1-4-eel","gb-protein-eelect","gb-protein-egb","gb-protein-esurf",\
            "gb-ligand-etot","gb-ligand-1-4-eel","gb-ligand-eelec","gb-ligand-egb","gb-ligand-esurf"])
  
  
    df_data = pd.DataFrame.from_dict(mut_cmplx, orient='index')
    df_data.columns= cols
    df_data.to_csv("all_amino_nn_data.csv")

if __name__ == '__main__':
    main()