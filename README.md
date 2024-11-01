---
editor_options: 
  markdown: 
    wrap: sentence
---

# FIBOS Python (BETA)

The Occluded Surface (OS) algorithm is a widely used approach for analyzing atomic packing in biomolecules. 
Here, we introduce **fibos**, an R and Python package that extends the OS methodology with enhancements. 
It integrates efficient Fortran code from the original OS implementation and introduces an innovation: 
the use of Fibonacci spirals for surface point distribution. This modification reduces anisotropy and 
ensures a more uniform and even distribution of surface dots, improving the accuracy
of the algorithm.

R fibos version: https://github.com/insilico-unifei/fibos-R.git.

## Operating Systems

FIBOS was designed to be multiplatform and run on Linux, Windows and Mac.\
Tested on:

- **Linux**: Ubuntu ($\geq$ 20.04)
- **Windows**: Windows 11
- **Mac**: MacOS 15.0.1

## Compilers

- gfortran
- gcc

## Python versions

Tested on: 3.9, 3.10

## Instalations

### Preliminary:

Some preliminary actions according to OS:

#### Linux (Ubuntu)

Install gfortran, Python dev and venv:
```
$ sudo apt install gfortran
$ sudo apt install python3.x-dev python3.x-venv
```
where "x" is your Python version.

#### Windows
Install Desktop Development with C++ from Microsoft C++ Build Tools from:
```
https://visualstudio.microsoft.com/visual-cpp-build-tools/ 
```

Install gfortran (version 13.2+):
```
http://www.equation.com/servlet/equation.cmd?fa=fortran
```

Install git from:
```
https://git-scm.com/downloads
```

#### MacOS

Install Homebrew:
```
$ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/
HEAD/install.sh)”
```

In your shell, set the PATH to include the Homebrew bin folder by adding it into 
the .zshrc file

```
export PATH= "/path/to/homebrew/bin:$PATH"
```
where "/path/to/homebrew/bin" is the actual homebrew path in your system. So, reload it:

```
$ source ~/.zshrc
```

Install gfortran and gcc from:

```
$ brew install gfortran
$ brew install gcc

```

### Virtual environment (Venv) 

It is highly recommended to work with virtual environments (Venv or Conda) in Python. 
We show how create Venv:

```
# From shell terminal, in project directory:
# Create a virtual environment ".venv"
$ python3.x -m venv .venv
(where "x" is your Python version.)

# Activate the virtual environment:

# Mac/Linux
$ source .venv/bin/activate

# Windows
$ .venv\Scripts\activate

# The prompt will change to something like:
(.venv)$  
```

### Basic Instalations

Install packages in requirements.txt (in .venv):

```
(.venv)$ pip install -r requirements.txt
```

Install fibos:

```         
(.venv)$ pip install git+https://github.com/insilico-unifei/fibos-py
```

## Main functions:

1.  **`occluded_surface(pdb, method = "FIBOS")`**: Implements the Occluded Surface 
algorithm, generating points, areas, and normals for each atom. It accepts the path 
to a PDB file and a method selection — either the classic 'OS' or the default 'FIBOS'. 
The function returns the results as a table and creates a file named 
`prot_PDBid.srf` in the `fibos_file` directory.

2.  **`osp(file)`**: Implements the Occluded Surface Packing (OSP) metric for 
each residue. Accepts a path to an .srf file generated by `occluded_surface` as a 
parameter and returns the results as a table summarized by residue.

### Quickstart:

```Python
import fibos
import os

# Calculate FIBOS per atom and create .srf files in fibos_files folder
pdb_fibos = fibos.occluded_surface("1fib", method="FIBOS")

# Show first 3 rows of pdb_fibos table
print(pdb_fibos.head(3))

#                     ATOM NUMBER_POINTS   AREA RAYLENGTH DISTANCE
# 0  GLN 1@N___>HIS 3@NE2_             6  1.287     0.791     5.49
# 1  GLN 1@N___>HIS 3@CE1_             1  0.200     0.894     6.06
# 2  GLN 1@N___>HIS 3@CG__             1  0.160     0.991     6.27

# Calculate OSP metric per residue from .srf file in fibos_files folder
pdb_osp = fibos.osp(os.path.join("fibos_files","prot_1fib.srf"))

# Show first 3 rows of pdb_osp table
print(pdb_osp.head(3))

#    Resnum Resname     OS  os*[1-raylen]    OSP
# 0       1     GLN  36.81          21.94  0.157
# 1       2     ILE  49.33          36.13  0.317
# 2       3     HIS  64.14          43.17  0.335
```

### A more complex example:

```Python
import fibos
import os
from Bio.PDB import PDBList
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# Auxiliary function to calculate occluded surface
def occluded_surface_worker(pdb_path, method):
    return fibos.occluded_surface(pdb_path, method=method)

# Get PDB file from RCSB and put it into the PDB folder  
# Rename the file appropriately and return path to it 
# (i.e., PDB/prot_8rxn.ent -> PDB/8rxn.pdb)
def get_pdb(id, path = "."):
    pdbl = PDBList()
    new_path = os.path.join(path, f"{id.lower()}.pdb")
    if not os.path.exists(new_path):
        original_path = pdbl.retrieve_pdb_file(id.lower(), pdir=path, file_format='pdb')
        os.rename(original_path, new_path)
    return new_path

if __name__ == "__main__":

    # source of PDB files
    pdb_folder = "PDB"

    # fibos folder output
    fibos_folder = "fibos_files"

    # Create PDB folder if it does not exist
    os.makedirs(pdb_folder, exist_ok=True)
    
    # Prevent overwriting of fibos folder output
    if os.path.exists(fibos_folder):
        raise Exception(f"{fibos_folder} folder exists, rename or remove it!")

    # PDB ids list
    pdb_ids = ["8RXN", "1ROP"]

    # Get PDB files from RCSB and put them into the PDB folder  
    pdb_paths = list(map(lambda pdb_id: get_pdb(pdb_id, path=pdb_folder), pdb_ids))
    print(pdb_paths)
    
    # Detect number of physical cores and update cores according to pdb_ids size
    my_ideal_cores = min(os.cpu_count(), len(pdb_ids))

    # Calculate in parallel FIBOS per PDBid 
    # Create .srf files in fibos_files folder
    # Return FIBOS tables in pdb_fibos list

    worker_with_params = partial(occluded_surface_worker, method="FIBOS")
    with ProcessPoolExecutor() as executor:
        pdb_fibos = list(executor.map(worker_with_params, pdb_paths))

    # Show first 3 rows of first pdb_fibos table
    print(pdb_fibos[0].head(3))
    
    # Prepare paths for the generated .srf files in folder fibos_files
    srf_paths = list(map(lambda pdb_id: os.path.join(fibos_folder, f"prot_{pdb_id.lower()}.srf"), pdb_ids))
    print(srf_paths)
    
    # Calculate OSP metric by residue
    # Return OSP tables in pdb_osp list
    pdb_osp = list(map(lambda srf_path: fibos.osp(srf_path), srf_paths))
    
    # Show first 3 rows of the first pdb_osp table
    print(pdb_osp[0].head(3))
    
# OBS: If you need to run this example in a Jupyter Notebook, move the 
# occluded_surface_worker function to a "fun.py" file and import it as 
# "from fun import occluded_surface_worker"
```

### Case Study:
[Here](https://github.com/insilico-unifei/fibos-R-case-study-supp.git) we show a 
case study (currently only in R), aiming to compare the packing density between experimentally 
determined structures and the same structures predicted by AlphaFold (AF).

## Authors

-   Carlos Silveira ([carlos.silveira\@unifei.edu.br](mailto:carlos.silveira@unifei.edu.br))\
    Herson Soares ([d2020102075\@unifei.edu.br](mailto:d2020102075@unifei.edu.br))\
    Institute of Technological Sciences,\
    Federal University of Itajubá,\
    Campus Itabira, Brazil.

-   João Romanelli ([joaoromanelli\@unifei.edu.br](mailto:joaoromanelli@unifei.edu.br)) \
    Institute of Applied and Pure Sciences, \
    Federal University of Itajubá, \
    Campus Itabira, Brazil.

-   Patrick Fleming ([Pat.Fleming\@jhu.edu](mailto:Pat.Fleming@jhu.edu)) \
    Thomas C. Jenkins Department of Biophysics, \
    Johns Hopkins University, \
    Baltimore, MD, USA

## References

Fleming PJ, Richards FM. Protein packing: Dependence on protein size, secondary structure and amino acid composition. J Mol Biol 2000;299:487–98.

Pattabiraman N, Ward KB, Fleming PJ. Occluded molecular surface: Analysis of protein packing. J Mol Recognit 1995;8:334–44.

## Status

In Progress.
