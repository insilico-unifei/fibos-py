from Bio.PDB import PDBList
import os
import shutil
import renum75

def get_file(file):
    path = os.getcwd()
    if not(".pdb" in file):
        if not (os.path.exists(file+".pdb")):
            pdb = PDBList()
            pdb.retrieve_pdb_file(file, pdir = path, file_format='pdb')
            path = "pdb"+file+".ent"
            os.rename(path, file+".pdb")
        shutil.copy(file+".pdb", "temp.pdb")
    elif ".pdb" in file:
        shutil.copy(file, path)
        shutil.copy(file, "temp.pdb")
    last = clean("temp.pdb")
    return last

def clean(file):
    i = 0
    l1 = ["B", "2", "L"]
    l2 = ["A","B","C","D","E","F","G","H","I"]
    l3 = ["A","1","U"]
    l4 = ["HOH","PMS","FOR","ALK","ANI"]
    t = False
    with open(file,"r") as file, open("temp.cln","w") as new_file:
        for line in file:
            if line.startswith("ATOM"):
                if line[13]!="H" and line[12]!="H" and line[13:16]!="OXT" and line[13:15] !="2H" and line[13:15]!="3H" and line[13]!="D" and line[12]!="D":
                    if not(line[16] in l1) and not (line[26] in l2):
                        if (line[16] in l3):
                            line = line[:16]+" "+line[17:]
                        for item in l4:
                            if item in line:
                                t = True
                        if i == 0:
                            line = line[:23]+"  1"+line[26:]
                            i = 1
                        if line.__contains__("HSD"):
                            line = line.replace("HSD","HIS")
                        if line.__contains__("HSE"):
                            line = line.replace("HSE", "HIS")
                        if line.__contains__("OT1"):
                            line = line.replace("OT1", "O")
                        if line.__contains__("OT2"):
                            line = line.replace("OT2", "OXT")
                        if line.__contains__("CD ILE"):
                            line = line.replace("CD ILE", "CD1 ILE")
                        if t == False:
                            new_file.write(line)
        new_file.write("END")
    renum75.renum()
    with open("new.pdb","r") as file:
       for line in file:
           if not ("END" in line):
              last_line = line
    last_residue = int(last_line[22:26])
    os.remove("temp.pdb")
    os.rename("new.pdb", "temp.pdb")
    return(last_residue)