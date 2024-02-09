

import glob

def main()   :
    mutations = [] 
    delta_delta_g_gen = [] 
    delta_delta_g_gen_error = []
    delta_delta_g_gen_66 = [] 
    delta_delta_g_gen_66_error = []

    tmp_file_loc = "/*_dir/FINAL_RESULTS*"
    tmp_file_names = glob.glob(tmp_file_loc)

    for file_name in tmp_file_names:
        with open(file_name) as f :
            data = f.readlines()
            c = 0
            for line in data : 
                if "RESULT OF ALANINE SCANNING" in line :
                    if c == 0 :
                        mutations.append(line.split(":")[1].split()[0].replace("(" , "").replace(")", ""))
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

    
    for i in len(tmp_file_names) :
        print(mutations[i] )
        print(delta_delta_g_gen[i] )
        print(delta_delta_g_gen_error[i] )
        print(delta_delta_g_gen_66[i] )
        print(delta_delta_g_gen_66_error[i] )
    
    

if __name__ == '__main__':
    main()