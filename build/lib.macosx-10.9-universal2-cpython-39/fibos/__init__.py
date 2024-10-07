import os
import pkgutil
import sys
import platform
import shutil
name_pack = "fibos"


path_pack = pkgutil.get_loader(name_pack).get_filename()

path_pack = os.path.dirname(path_pack)

path_abs = os.path.abspath(path_pack)
sys.path.append(path_abs)

origem = path_abs +"\.libs"

if((platform.system() == "Windows") and os.path.exists(origem)):
       destino = path_abs
       #os.makedirs(destino)
       arquivos = os.listdir(origem)
       for arquivo in arquivos:
            caminho_origem = os.path.join(origem, arquivo)
            caminho_destino = os.path.join(destino, arquivo)
            shutil.move(caminho_origem, caminho_destino)

from .manage_os import occluded_surface
#from .read_os import read_prot
from .respak import osp
#from .respak import read_osp
