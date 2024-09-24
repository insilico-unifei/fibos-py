import pandas as pd
import re

def read_prot(prot_file):
    lines = []
    name_columns = ['INTERACTION', 'NUMBER_POINTS', 'AREA', 'RAYLENGTH', 'DISTANCE']
    with open(prot_file,"r") as pfile:
        for item in pfile:
            if("INF" in item):
                item = re.sub('A2', ' ', item)
                item = re.sub('Rlen', ' ', item)
                item = re.sub('pts', ' ', item)
                item = re.sub('Dxx', ' ', item)
                item = re.sub('INF', ' ', item)
                item = re.sub(r'\s+', ' ', item)
                lines.append(item)

    prot = pd.DataFrame(lines,columns = ["TEMP"])
    prot = prot['TEMP'].str.split(' ',expand=True)
    prot[1] = prot[1] + " " + prot[2] + " " + prot[3]
    prot = prot.drop(columns=[0])
    prot = prot.drop(columns=[8])
    prot = prot.drop(columns=[2])
    prot = prot.drop(columns=[3])
    prot.columns = name_columns
    return prot