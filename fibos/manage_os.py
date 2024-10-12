import os
import cleaner
import main_intermediary
import pkgutil
import shutil
import tempfile
from pathlib import Path
from Bio.PDB import PDBParser
from Bio import PDB
from read_os import read_prot
from folders_manipulate import create_folder
from folders_manipulate import change_files
from respak import osp

MAX_RES = 10000
MAX_AT = 50000
IRESF = 1


def occluded_surface(pdb,method = "FIBOS"):
    remove_files()
    #create_folder(pdb)
    source_path = os.getcwd()
    change = False
    name_pdb = pdb
    if not os.path.exists("fibos_files"):
        os.makedirs("fibos_files")
    if not (pdb.endswith(".pdb")):
            arq_aux = pdb + ".pdb"
            if os.path.exists(arq_aux):
                os.remove(arq_aux)
    else:
            change = True
            name_pdb = name_pdb[-8:]
            if (os.path.exists(pdb) == False):
                raise FileNotFoundError(f"File not Found: {pdb}")
            #parser = PDB.PDBParser()
            #estrutura = parser.get_structure("protein", pdb)
    with tempfile.TemporaryDirectory() as temp_dir:
        if(change == True):
            #io = PDB.PDBIO()
            #io.set_structure(estrutura)
            #output_pdb = os.path.join(temp_dir, os.path.basename(name_pdb))
            #io.save(output_pdb)
            shutil.copy(pdb,temp_dir)
        os.chdir(temp_dir)
        pdb = pdb.lower()
        name_pack = "fibos"
        path_pack = pkgutil.get_loader(name_pack).get_filename()
        path_pack = os.path.dirname(path_pack)
        path_abs = os.path.abspath(path_pack)
        path_abs = path_abs+"/radii"
        shutil.copy(path_abs,".")
        iresl = cleaner.get_file(name_pdb)
        method = method.upper()
        meth = 0
        if method == "OS":
            meth = 1
        elif method == "FIBOS":
            meth = 2
        else:
            print("Wrong Method")
        file_remove = pdb
        if pdb.endswith(".pdb"):
            file_remove = pdb.replace(".pdb","")
        ray_remove = "raydist_"+file_remove+".lst"
        pack_remove = "prot_"+file_remove+".pak"
        file_remove = "prot_"+file_remove+".srf"
        if os.path.exists(file_remove):
            os.remove(file_remove)
        if os.path.exists(ray_remove):
            os.remove(ray_remove)
        if os.path.exists(pack_remove):
            os.remove(pack_remove)
        main_intermediary.call_main(IRESF,iresl,MAX_RES,MAX_AT,meth)
        remove_files()
        file_name = change_files(pdb)
        from_path = source_path+"/fibos_files"
        for arquivo in os.listdir(temp_dir):
            caminho_arquivo = os.path.join(temp_dir,arquivo)
            if os.path.isfile(caminho_arquivo) and not arquivo.endswith(".pdb") and "radii" not in arquivo:
                shutil.copy2(caminho_arquivo,from_path)
        os.chdir(source_path)
    #file_name = "fibos_files/"+file_name
    return read_prot(os.path.join("fibos_files", file_name))


def remove_files():
    path = os.getcwd()
    extensions = ['.ms','.txt','.inp']
    files = [file for file in os.listdir(path) if any(file.endswith(ext) for ext in extensions)]
    if len(files)>0:
        for file in files:
            os.remove(os.path.join(path, file))
    if os.path.exists(os.path.join(path,"fort.6")):
        os.remove(os.path.join(path,"fort.6"))
    if os.path.exists(os.path.join(path, "part_i.pdb")):
        os.remove(os.path.join(path, "part_i.pdb"))
    if os.path.exists(os.path.join(path, "part_v.pdb")):
        os.remove(os.path.join(path, "part_v.pdb"))
    #os.remove("temp.pdb")
    if os.path.exists("temp.cln"):
        os.remove("temp.cln")
