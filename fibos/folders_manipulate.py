import os
import shutil

def create_folder(pdb):
    if ((os.path.exists("fibos_files") == False) and (os.path.basename(os.getcwd())!= "fibos_files")):
        os.mkdir("fibos_files")

def change_files(pdb_name):
    if(".pdb" in pdb_name):
        pdb_name = pdb_name.removesuffix(".pdb")
        if(len(pdb_name)>4):
            pdb_name = pdb_name[-4:]
    name_pdb = pdb_name
    name_raydist = pdb_name
    if(os.path.exists("prot.srf")):
        shutil.copy("prot.srf", "fibos_files")
        os.remove("prot.srf")
        name_pdb = "fibos_files/prot_"+name_pdb+".srf"
        os.rename("fibos_files/prot.srf", name_pdb)

    if (os.path.exists("raydist.lst")):
        shutil.copy("raydist.lst", "fibos_files")
        os.remove("raydist.lst")
        name_raydist = "fibos_files/raydist_" + name_raydist + ".lst"
        os.rename("fibos_files/raydist.lst", name_raydist)
    return (name_pdb)