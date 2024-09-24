import os

import pandas as pd
import fibos
import statistics

column_names = ['Resnum', 'Resname','OS','os*[1-raylen]','OSP']

def make():

    with open("table_4.dat", "r") as list_pdbs, open("all_paks_mean.dat","w") as means:
        for items in list_pdbs:
            items = items.strip()
            fibos.occluded_surface(items,"fibos")
            fibos.respak("prot.srf")
            respak_pdb = pd.read_csv("prot.pak", delim_whitespace=True, names = column_names)
            respak_pdb = respak_pdb['OSP'][1:]
            aux = respak_pdb.astype(float)
            mean_pak = statistics.mean(aux)
            means.write(items+" "+str(mean_pak)+"\n")
            name_pak = items+"_pak_py.dat"
            with open(name_pak,"a") as file_pak:
                for it in respak_pdb:
                    file_pak.write(str(it)+"\n")

if __name__ == '__main__':
    make()
    os.remove("prot.pak")
    os.remove("prot.srf")
