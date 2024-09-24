import os
import cleaner
import main_intermediary
import pkgutil
import shutil
from read_Os import read_prot
from folders_manipulate import  create_folder
from folders_manipulate import rename_file
from respak import respack

MAX_RES = 10000
MAX_AT = 50000
IRESF = 1


def occluded_surface(pdb,method = "FIBOS"):
    create_folder(pdb)
    pdb = pdb.lower()
    name_pack = "fibos"
    path_pack = pkgutil.get_loader(name_pack).get_filename()
    path_pack = os.path.dirname(path_pack)
    path_abs = os.path.abspath(path_pack)
    path_abs = path_abs+"/radii"
    shutil.copy(path_abs,".")
    iresl = cleaner.get_file(pdb)
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
    file_name = rename_file(pdb)
    pak_file = file_name
    file_name = "prot_"+file_name+".srf"
    os_file = read_prot(file_name)
    res_file = respack(file_name)
    pak_file = "prot_"+pak_file+".pak"
    os.rename("prot.pak", pak_file)
    return (os_file, res_file)

def remove_files():
    path = os.getcwd()
    extensions = ['.ms','.txt','.inp']
    files = [file for file in os.listdir(path) if any(file.endswith(ext) for ext in extensions)]
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
