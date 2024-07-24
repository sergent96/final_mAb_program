#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 16:15:47 2024

@author: jasonhwang
"""

import pandas as pd
from CifFile import ReadCif

g=open("currData.txt").read().split()
print(g)
chain=str(g[1])
name=str(g[2])

path=str(g[3])


# Path to your CIF file
cif_file_path = path
# Output CSV file path
csv_file_path = f'{name}_map.csv'

# Read the CIF file
cif_data = ReadCif(cif_file_path)

# Assuming the CIF file contains a single data block
data_block = list(cif_data.keys())[0]

# Extract relevant data from the CIF file
data = cif_data[data_block]

# Convert the data to a pandas DataFrame
# This example assumes you want to extract atomic positions (you may need to adjust this)
group=data["_atom_site.group_PDB"]
num=data["_atom_site.id"]
sym=data["_atom_site.type_symbol"]
atom_id=data["_atom_site.label_atom_id"]
sep=data["_atom_site.label_alt_id"]
residue=data["_atom_site.label_comp_id"]
residue_num=data["_atom_site.label_seq_id"]


# Create a DataFrame
df = pd.DataFrame({
    'atom num': num,
    'symbol': sym,
    'aa_id': atom_id,
    'residue': residue, 
    "residue num":residue_num
})

# Save the DataFrame to a CSV file
df.to_csv(csv_file_path, index=False)

print(f"CIF file has been converted to CSV and saved as {csv_file_path}")
