---
editor_options: 
  markdown: 
    wrap: sentence
---

# Phyton fibos (BETA)

The Occluded Surface (OS) algorithm is a widely used approach for analyzing atomic packing in biomolecules. 
Here, we introduce **fibos**, an R and Python package that extends the OS methodology with enhancements. 
It integrates efficient Fortran code from the original OS implementation and introduces an innovation: 
the use of Fibonacci spirals for surface point distribution. This modification reduces anisotropy and 
ensures a more uniform and even distribution of surface dots, improving the accuracy
of the algorithm.

R fibos version: https://github.com/insilico-unifei/fibos-R.git.

## Operating Systems

fibos was designed to be multiplatform and run on Linux, Windows and Mac.
However, it has been tested on the following versions:

- **Linux**: Ubuntu ($\geq$ 20.04)
- **Windows**: Windows 11
- **Mac**: MacOS 15.0.1

## Instalation

### Python versions

Tested on: 3.9, 3.10

### Preliminar packages:

Some preliminar packages may be necessary:

```
$ sudo apt-get install python3.x-dev # Ubuntu
```

### Virtual environment (venv) 

It is highly recommended to work with virtual environments in Python

```
# From shell terminal, in project directory:
# (creates a virtual environment ".venv")
$ python3 -m venv .venv

# Activate the virtual environment:

# Mac/Linux
$ source .venv/bin/activate

# Windows
$ .venv\Scripts\activate

# The prompt will change to something like:
(.venv)$  
```

### Instalations

These additional Python packages may be required (install in .venv):

```
(.venv)$ pip install testresources 
(.venv)$ pip install biopython 
```

Install fibos:

```         
pip install git+https://github.com/insilico-unifei/fibos-py
```

## Main functions:

1.  **`occluded_surface(pdb, method = "FIBOS")`**: Implements the Occluded Surface 
algorithm, generating points, areas, and normals for each atom. It accepts the path 
to a PDB file and a method selection—either the classic 'OS' or the default 'FIBOS'. 
The function returns the results as a table and creates a file named 
prot_PDBid.srf in the fibos_file directory.

1.  **`osp(file)`**: Implements the Occluded Surface Packing (OSP) metric for 
each residue. Accepts a path to an SRF file generated by occluded_surface as a 
parameter and returns the results as a table summarized by residue.

### A simple example:

```     
import os
import shutil
import pkgutil
import fibos
from Bio.PDB import PDBList
from concurrent.futures import ThreadPoolExecutor

# Create folder if it does not exist
folder = "PDB"
os.makedirs(folder, exist_ok=True)

# PDB ids
pdb_ids = ["8RXN", "1ROP"]

# Initialize PDBList object
pdbl = PDBList()

# Get PDB from RCSB, put it in folder and return path to it
pdb_paths = [pdbl.retrieve_pdb_file(pdb_id, pdir=folder, file_format='pdb') for pdb_id in pdb_ids]

with ThreadPoolExecutor(max_workers=2) as executor:
    pdb_fibos = list(executor.map(lambda x: fibos.occluded_surface(x, method="FIBOS"), pdb_paths))
```

### More complex example:
[Here](https://github.com/insilico-unifei/fibos-R-case-study-supp.git) we show a 
case study (in R only), aiming to compare the packing density between experimentally 
determinedstructures and the same structures predicted by AlphaFold (AF).

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