import respak75
import os
import shutil
import pandas as pd
import re



def osp(prot_file):
    if not(os.path.exists(prot_file)):
        raise FileNotFoundError("File not Found: "+prot_file)
    if not(os.path.exists("fibos_files")):
        os.create_folder("fibos_files")
    else:
        if(prot_file!="prot.srf"):
            shutil.copy(prot_file,"prot.srf")
        respak75.respak()
        if(prot_file!="prot.srf"):
            os.remove("prot.srf")
        prot_name = prot_file.removesuffix(".srf")
        prot_name = "prot_"+prot_name[-4:]
        prot_name = prot_name+".pak"
        os.rename("prot.pak",prot_name)
        #os.rename("prot.pak",prot_name)
        shutil.copy2(prot_name,"fibos_files")
        os.remove(prot_name)
        prot_name = "fibos_files/"+prot_name
    return (pd.read_table(prot_name, header=0, sep=r'\s+'))

def read_osp(prot_file):
    if not(os.path.exists(prot_file)):
        raise FileNotFoundError("File not Found: "+prot_file)
    return (pd.read_table(prot_file, header=0, sep=r'\s+'))