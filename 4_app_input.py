#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:35:01 2024

@author: jasonhwang
"""
import biopandas
import pandas as pd
from CifFile import ReadCif


g=open("currData.txt").read().split()
pdb=str(g[0])
chain=str(g[1])

altered_protein=[]
altered_protein.append(["A", "wt_rbd"])
altered_protein.append(["B", "vh"])
altered_protein.append(["C", "vl"])


path=f"/Users/jasonhwang/Documents/pdb_trim/{pdb}_trimmed.cif"
# Path to your CIF file
cif_file_path = path
# Output CSV file path
csv_file_path = f'{pdb}_map.csv'
# Read the CIF file
cif_data = ReadCif(cif_file_path)
# Assuming the CIF file contains a single data block
data_block = list(cif_data.keys())[0]
# Extract relevant data from the CIF file
data = cif_data[data_block]
group=data["_atom_site.group_PDB"]
num=data["_atom_site.id"]
sym=data["_atom_site.type_symbol"]
atom_id=data["_atom_site.label_atom_id"]
sep=data["_atom_site.label_alt_id"]
residue=data["_atom_site.label_comp_id"]
residue_num=data["_atom_site.label_seq_id"]


# Create a DataFrame
ppdb = pd.DataFrame({
    'atom num': num,
    'symbol': sym,
    'aa_id': atom_id,
    'residue': residue, 
    "residue num":residue_num
})

num=ppdb["atom num"]
resname=ppdb["residue"]
resnum=ppdb["residue num"]
pdbdata=pd.DataFrame({"atom_num":num, "res id":resnum, "residue":resname})

data=pd.read_csv(f"/Users/jasonhwang/Documents/pdb_trim/{pdb}_{chain}.csv", sep="\s+", header=None)
#adds columns to df
data.columns=["dist", "vh atoms", "rbd atoms"]
print(data)

#NOTE: This creates intial df
vh=[]
vhcode=[]
app_vhcode=[]
res=[]
dist=[]
#id=[]

#goes through pdb data and inputs residue and atom data given atom id
rbd=[]
rbdcode=[]
app_rbdcode=[]
res_rbd=[]

print(ppdb)
for index, row in data.iterrows():
    vh.append(ppdb.iloc[int(row["vh atoms"])-1, 0])
    vhcode.append(ppdb.iloc[int(row["vh atoms"])-1, 3])
    app_vhcode.append(ppdb.iloc[int(row["vh atoms"])-1, 3])
    res.append(ppdb.iloc[int(row["vh atoms"])-1, 4])
    y = ppdb.iloc[int(row["vh atoms"])-1, 4]
    #print(ppdb.iloc[int(row["vh atoms"])-1, 0])
    x=ppdb.iloc[int(row["vh atoms"])-1, 0]
    currPos = int(x)
print(res)
for index, row in data.iterrows():
    rbd.append(ppdb.iloc[int(row["rbd atoms"])-1, 0])
    rbdcode.append(ppdb.iloc[int(row["rbd atoms"])-1, 3])
    res_rbd.append(ppdb.iloc[int(row["rbd atoms"])-1, 4])
########################
#dict of amino acids SHM
amino_acids = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
#labels data df
data.columns=["dist", "vh atoms", "rbd atoms"]

#dict of amino acids AA code
aa_code={
    'ARG': ["CB", "CZ", "NH1", "NH2", "CD", "NE", "CG"],
    'ALA': ["CB"],
    'ASN': ["CB", "CG", "OD1", "ND2"],
    'ASP': ["CB", "CG", "OD1", "OD2"],
    'CYS': ["CB", "SG"],
    'GLN': ["CB", "CD", "OE1", "NE2", "CG", "OE2"],
    'GLU': ["CB", "CD", "OE1", "OE2"],
    'GLY': ["CA"],
    'HIS': ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    'ILE': ["CB", "CG1", "CG2", "CD"],
    'LEU': ["CB", "CG", "CD1", "CD2"],
    'LYS': ["CB", "CE", "NZ", "CD", "CG"],
    'MET': ["CB", "CG", "SD", "CE"],
    'PHE': ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    'PRO': ["CB", "CG", "CD"],
    'SER': ["CB", "OG"],
    'THR': ["CB", "OG1", "CG2"],
    'TRP': ["CB", "NE1", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    'TYR': ["CB", "OH", "CE1", "CE2", "CZ", "CG", "CD1", "CD2"],
    'VAL': ["CB", "CG1", "CG2"],

}

#assigns amino acid code for amino acids for every elt
for i in range(len(vhcode)):
  vhcode[i]=amino_acids[vhcode[i]]
for j in range(len(rbdcode)):
  rbdcode[j]=amino_acids[rbdcode[j]]

for i in range(len(app_vhcode)):
  app_vhcode[i]=aa_code[app_vhcode[i]]
for j in range(len(app_rbdcode)):
  app_rbdcode[j]=aa_code[app_rbdcode[j]]

#creates new df with combined data
output=pd.DataFrame({"vh aa code":app_vhcode, "res code":app_vhcode,  "res location": res})
appegio_input=[]
for i in range(len(output)):
  k=int(res[i])
  if chain==altered_protein[1][1]: #vh
    print("asdfasd")
    for j in range(len(app_vhcode[i])):
      appegio_input.append(f"/{altered_protein[1][0]}/{k}/{app_vhcode[i][j]}")
  if chain==altered_protein[2][1]: #vh
    for j in range(len(app_vhcode[i])):
      appegio_input.append(f"/{altered_protein[2][0]}/{k}/{app_vhcode[i][j]}")
for i in range(len(res_rbd)):
  l=int(res_rbd[i])
  appegio_input.append(f"/{altered_protein[0][0]}/{l}/")
no_duplicates_appegio = list(set(appegio_input))

file = open(f'appegio_output_{pdb}_{chain}.txt','w')
for line in no_duplicates_appegio:
	file.write(line+"\n")
file.close()

print("program completed")