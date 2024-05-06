import os
import argparse
import logging
import pandas as pd
import glob
import time
import shutil
from seaborn import color_palette

# for converting output of DALIouput parser to chimeraX scripts, input folder is folder of PDBs, output folder is where you want the pictures/scripts to go, query pdb is the original pdb file, all pdbs in input folder are aligned to query pdb, this allows for easier crossreferencing of residue location

def write_chimeraX_scripts(input_dir, output_dir, query_pdb):
    cols = color_palette("bright", 20).as_hex()
    # pdb names 
    input_name = os.path.basename(input_dir)
    pdb_paths = [f for f in glob.glob(os.path.join(input_dir, "*.pdb"))]
    folder_names = [os.path.basename(pdb_path) for pdb_path in pdb_paths]
    # pdb for pdb_names
    print(folder_names)
    
    for folder in folder_names:
        
        # folder_p = os.path.join(input_dir, folder)
        chain = folder.rsplit("-", 1)[1]
        # pic_dir = os.makedirs(os.path.join(folder, f"{folder}_chimeraX_rmsd_pics"), exist_ok=True)
        command_list = []
        
        # Move variable definitions outside of the for loop
        pdb_name = ""
        pdb_chain = ""
        pdb_code = ""
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            
        for i, pdb_p in enumerate(pdb_paths):
            print(pdb_p)
            random_color = cols[i % len(cols)]
            par_folder = os.path.dirname(pdb_p)
            gpar_folder = os.path.dirname(par_folder)
            pdb_name = os.path.basename(pdb_p).rsplit(".pdb", 1)[0]
            pic_dir = os.path.join(output_dir, f"{pdb_name}_chimeraX_rmsd_pics")
            if not os.path.exists(pic_dir):
                os.makedirs(pic_dir, exist_ok=True)
            pdb_chain = pdb_name.split("-", 1)[1]
            pdb_code = pdb_name.split("-", 1)[0]
            print(f'pic_dir: {pic_dir}\n pdb_name: {pdb_name} \n pdb_chain: {pdb_chain} \n pdb_code: {pdb_code}')
            
            open_model = f"open {pdb_p}; hide #2 target acs; hide solvent; hide ligand; show #2/{pdb_chain} cartoon ; SetUpSession; mmaker #2/{pdb_chain} to #1 ; hide #1 target acs"
            view_model = f"show #2/{pdb_chain} cartoon ; view #2/{pdb_chain}"
            four_styles = f"fourstyles #2/{pdb_chain} '{pic_dir}/{pdb_name}_{pdb_chain}'"
            # chain_show_string = f"sel /{pdb_chain} | #1 ; show sel cartoons ; hide ~sel target acs ; select clear "
            # styler = f"med_licorice; col #1 grey ; col #1/{pdb_chain} black; col #2 {random_color} ; transparency #1 50 target c"
            # pics1 = f"spin_pic {os.path.join(pic_dir, f'{pdb_code}_dali')}"
            # reset_style = "transparency 0 target c; show #1 surface ; col #1/A #75c376; col #1/B #a7cee2 ; col #1/C #2078b4 "
            # style2 =  f"med_licorice; hide #1 cartoon ; show #1/{pdb_chain} cartoon ; thin_licorice_chain #1/{pdb_chain}; col #1/{pdb_chain} black target c ; hide #2 target as; transparency #1 80 target s"
            # pics2 = f"spin_pic {os.path.join(pic_dir, f'{pdb_code}_dali_s2')}"
            # reset_style2 = "transparency 0 target sc; hide #1 surface ; show #1 cartoon; col #1/A #75c376; col #1/B #a7cee2 ; col #1/C #2078b4 "
            # del_pdb = "del #2 \n"
            close_model = "close #2 \n"
            
            commands = [open_model, view_model, four_styles]
            commands.append(close_model)
            # commands = [open_string, chain_show_string, styler,  reset_style, style2, reset_style2, del_pdb]

            # commands = [open_string, chain_show_string, styler, pics1, reset_style, style2, pics2, reset_style2, del_pdb]
          
            command_list.append("\n".join(commands))
        
        # TODO refactor this so that it doesnt execute mutiple times
        
        file_out = os.path.join(output_dir, f"{input_name}_Generate_DALI_HITs_photos.cxc")
        query_pic_dir = os.path.join(output_dir, f"query_{input_name}_pics")
        if not os.path.exists(query_pic_dir):
            os.makedirs(query_pic_dir, exist_ok=True)
        print(file_out)
        with open(file_out, "w") as file:
            file.write(f"# ChimeraX SetUp and Run RMSD photos for {input_name} \n")
            file.write("# Assumes you are in an open session with your query structure \n")
            file.write(f"cd ~/thesis_code \n open 'chimeraX_functions/SetUpSession-ChX.py' \n")
            file.write(f"open '{query_pdb}' ; hide solvent ; hide ligand ; hide atoms ; SetUpSession \n")
            # file.write(f"open {os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src/SetUpSession-ChX.py')} \n")
            file.write(f"fourstyles #1 '{query_pic_dir}/{input_name}' \n")
            
            for command in command_list:
                file.write(command)
    
    return 


def main():
    parser = argparse.ArgumentParser(description="Generate ChimeraX scripts for DALI hits")
    parser.add_argument("input_dir", help="Input directory containing PDB files")
    parser.add_argument("output_dir", help="Output directory for ChimeraX scripts")
    parser.add_argument("query_pdb", help="Original query PDB file")
    args = parser.parse_args()
    
    start_time = time.time()
    write_chimeraX_scripts(args.input_dir, args.output_dir, args.query_pdb)
    end_time = time.time()
    
    logging.info(f"Complete, Approximate Time Taken: {end_time - start_time:.2f} sec")
    return


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()