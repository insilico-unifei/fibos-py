import respak75
import os
import shutil
import pandas as pd

def respack(prot_file):
    if not(os.path.exists(prot_file)):
        print("Prot File not found!")
    else:
        if(prot_file!="prot.srf"):
            shutil.copy(prot_file,"prot.srf")
        respak75.respak()
        if(prot_file!="prot.srf"):
            os.remove("prot.srf")
    return (pd.read_table("prot.pak", header=0, sep=r'\s+'))