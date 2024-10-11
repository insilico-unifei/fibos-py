import pandas as pd
import fibos
import statistics

column_names = ['Resnum', 'Resname','OS','os*[1-raylen]','OSP']

def make():
    with open("table_4.txt", "r") as list_pdbs, open("all_paks_mean.dat","w") as means:
        for items in list_pdbs:
            fibos.occluded_surface(items,"os")
            fibos.respak()
            respak_pdb = pd.read_csv("prot.pak", delim_whitespace=True, names = column_names)
            respak_pdb = respak_pdb['OSP']
            mean_pak = statistics.mean(respak_pdb)
            means.write(mean_pak)
            name_pak = items+"_pak_py.dat"
            with open(name_pak,"w") as file_pak:
                file_pak.writelines(respak_pdb)