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

FIBOS was designed to be multiplatform and run on Linux, Windows and Mac.
However, it has been tested on the following versions:

- **Linux**: Ubuntu ($\geq$ 20.04)
- **Windows**: Windows 11
- **Mac**: MacOS 15.0.1

## Python versions

Tested on: 3.9, 3.10

## Instalations

### Preliminary packages:

Some preliminary packages may be necessary:

```
# Ubuntu (where "x" is your Python version)
$ sudo apt install python3.x-dev
$ sudo apt install python3.x-venv
```

### Virtual environment (venv) 

It is highly recommended to work with virtual environments in Python

```
# From shell terminal, in project directory:
# Create a virtual environment ".venv"
$ python3.x -m venv .venv

# Activate the virtual environment:

# Mac/Linux
$ source .venv/bin/activate

# Windows
$ .venv\Scripts\activate

# The prompt will change to something like:
(.venv)$  
```

### Basic Instalations

These additional Python packages may be required (install in .venv):

```
(.venv)$ pip install biopython 
```

Install fibos:

```         
(.venv)$ pip install git+https://github.com/insilico-unifei/fibos-py
```

## Main functions:

1.  **`occluded_surface(pdb, method = "FIBOS")`**: Implements the Occluded Surface 
algorithm, generating points, areas, and normals for each atom. It accepts the path 
to a PDB file and a method selection — either the classic 'OS' or the default 'FIBOS'. 
The function returns the results as a table (tibble) and creates a file named 
prot_PDBid.srf in the fibos_file directory.

1.  **`osp(file)`**: Implements the Occluded Surface Packing (OSP) metric for 
each residue. Accepts a path to an SRF file generated by occluded_surface as a 
parameter and returns the results as a table summarized by residue.

### Quickstart:

```         
import fibos

# Calculate FIBOS per atom and create .srf files in fibos_files folder
pdb_fibos = fibos.occluded_surface("1fib", method="FIBOS")

# Show first 3 rows of pdb_fibos table
print(pdb_fibos.head(3))

#              INTERACTION NUMBER_POINTS   AREA RAYLENGTH DISTANCE
# 0  GLN 1@N___>HIS 3@NE2_             6  1.287     0.791     5.49
# 1  GLN 1@N___>HIS 3@CE1_             1  0.200     0.894     6.06
# 2  GLN 1@N___>HIS 3@CG__             1  0.160     0.991     6.27

# Calculate OSP metric per residue from .srf file in fibos_files folder
pdb_osp = fibos.osp("fibos_files/prot_1fib.srf")

# Show first 3 rows of pdb_osp table
print(pdb_osp.head(3))

#    Resnum Resname     OS  os*[1-raylen]    OSP
# 0       1     GLN  36.81          21.94  0.157
# 1       2     ILE  49.33          36.13  0.317
# 2       3     HIS  64.14          43.17  0.335
```

### A more complex example:

```     
import os
import fibos
from Bio.PDB import PDBList
from multiprocessing import Pool

# Auxiliary function to calculate occluded surface
def occluded_surface_worker(pdb_path, method):
    return fibos.occluded_surface(pdb_path, method=method)

if __name__ == "__main__":
    # Create PDB folder if it does not exist
    folder = "PDB"
    os.makedirs(folder, exist_ok=True)

    # PDB ids list
    pdb_ids = ["8RXN", "1ROP"]

    # Initialize PDBList object
    pdbl = PDBList()

    # Get PDB files from RCSB and put them into the PDB folder  
    # Rename files and return path to them 
    # (i.e., PDB/prot_8rxn.ent -> PDB/8rxn.pdb)
    pdb_paths = []
    for pdb_id in pdb_ids:
	    original_path = pdbl.retrieve_pdb_file(pdb_id, pdir=folder, file_format='pdb')
	    new_path = os.path.join(folder, f"{pdb_id.lower()}.pdb")
	    os.rename(original_path, new_path)
	    pdb_paths.append(new_path)

    print(pdb_paths)
    
    # Calculate in parallel FIBOS per atom per PDBid 
    # Create .srf files in fibos_files folder
    # Return FIBOS tables in pdb_fibos list
    with Pool(processes=os.cpu_count()) as pool:
        pdb_fibos = pool.starmap(occluded_surface_worker, [(pdb_path, "FIBOS") for pdb_path in pdb_paths])
        
    # Show first 3 rows of first pdb_fibos table
    print(pdb_fibos[0].head(3))
    
    # Prepare paths for the generated .srf files in folder fibos_files
    srf_paths = list(map(lambda pdb_id: os.path.join("fibos_files", f"prot_{pdb_id.lower()}.srf"), pdb_ids))
    print(srf_paths)
    
    # Calculate OSP metric by residue
    # Return OSP tables in pdb_osp list
    pdb_osp = list(map(lambda srf_path: fibos.osp(srf_path), srf_paths))
    
    # Show first 3 rows of the first pdb_osp table
    print(pdb_osp[0].head(3))
    
    # Rename the fibos_files folder to preserve it and prevent overwriting
    os.rename("fibos_files", "fibos_files_test")
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
