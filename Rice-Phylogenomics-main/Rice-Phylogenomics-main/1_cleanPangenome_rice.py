#%%


import os
import pandas as pd

def prefix_cell(cell, header):
    if pd.notna(cell) and cell.strip():  # Check for non-empty and non-whitespace cells
        items = cell.split("|")
        items = [header + item for item in items]
        return "|".join(items)
    else:
        return "NA"


work_dir = "/path/to/working/directory"
os.chdir(work_dir)
pangenome = "pangenome_table"
pangenome2 = "cleaned_pangenome_table"

with open(pangenome, "r") as f1:
    table = []
    for line in f1:
        # read text data into a data table
        line = line.strip()
        line = line.split("\t")
        e = [line[0],line[1],line[4]] + line[9:]
        table.append(e)
    
    df = pd.DataFrame(table)
    species_code = df.iloc[0] + "@" # extract species code from the first raw
    df.iloc[1:, 3:] = df.iloc[1:, 3:].apply(lambda col: col.apply(lambda cell: prefix_cell(cell, species_code[col.name])))
    
    df.to_csv(pangenome2, header=False, index=False, sep='\t')
# %%
