import glob
import pandas as pd

def main()   :
    print("working")
    mutations = [] 
    delta_delta_g_gen = [] 
    delta_delta_g_gen_error = []
    delta_delta_g_gen_66 = [] 
    delta_delta_g_gen_66_error = []

    tmp_file_loc = "*_dir/FINAL_RESULTS*A_0*"
    tmp_file_names = glob.glob(tmp_file_loc)
    for file_name in tmp_file_names:
        with open(file_name) as f :
            data = f.readlines()
            mutations.append(file_name.split("_")[-2] + "_" + file_name.split("_")[-1].split(".")[-2] )
            c = 0
            for line in data : 
                if "RESULT OF ALANINE SCANNING" in line :
                    if c == 0 :
                        
                        delta_delta_g_gen.append(\
                            line.split("=")[1].split("+/-")[0].strip()
                            )
                        delta_delta_g_gen_error.append(\
                            line.split("=")[1].split("+/-")[1].strip()
                            )
                        c +=1 
                    else: 
                        delta_delta_g_gen_66.append(\
                            line.split("=")[1].split("+/-")[0].strip()
                            )
                        delta_delta_g_gen_66_error.append(\
                            line.split("=")[1].split("+/-")[1].strip()
                            )



    data= {"mutations" :mutations, "delta_delta_g_gen" :delta_delta_g_gen, 
        "delta_delta_g_gen_error" : delta_delta_g_gen_error, 
        "delta_delta_g_gen_66" : delta_delta_g_gen_66,
        "delta_delta_g_gen_66_error" :delta_delta_g_gen_66_error
        }
  
    df_data = pd.DataFrame.from_dict(data)
    df_data.to_csv("all_amino_salt_con.csv")

if __name__ == '__main__':
    main()