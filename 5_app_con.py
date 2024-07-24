#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:35:28 2024

@author: jasonhwang
"""
path="/Users/jasonhwang/Downloads/alpha-twelve-pasta/alpha-twelve-pasta.contacts"


import biopandas
import pandas as pd
from CifFile import ReadCif


altered_protein=[]
altered_protein.append(["B", "vh"])
altered_protein.append(["C", "vl"])
altered_protein.append(["A", "wt_rbd"])



g=open("currData.txt").read().split()

pdb=str(g[0])
chain=str(g[1])
chained=chain
path_cif=f"/Users/jasonhwang/Documents/pdb_trim/{pdb}_trimmed.cif"
# Path to your CIF file
cif_file_path = path_cif
# Output CSV file path
csv_file_path = f'{pdb}_map.csv'
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
ppdb = pd.DataFrame({
    'atom num': num,
    'symbol': sym,
    'aa_id': atom_id,
    'residue': residue, 
    "residue num":residue_num
})
# Save the DataFrame to a CSV file
ppdb.to_csv(csv_file_path, index=False)


with open(f"arpeggio_data_{pdb}{chain}.txt", "w") as file:
    with open(path, "r") as f:
        x=f.read()
        file.write(x)

g=open(f"arpeggio_data_{pdb}{chain}.txt").read().split()
output=open(f"arpeggio_data_{pdb}{chain}_processed.csv", "w")

#group into 3, conditions for split
grouped_text_list=list(zip(*(iter(g),)*18))[1:]
for ele1, ele2, ele3, ele4, ele5, ele6, ele7, ele8, ele9, ele10, ele11, ele12, ele13, ele14, ele15, ele16, ele17, ele18 in grouped_text_list:
    if ele3 !=0 or ele4 !=0 or ele5 !=0 or ele6 !=0 or ele7 !=0 or ele8 !=0 or ele9 !=0 or ele10 !=0 or ele11 !=0 or ele12 !=0 or ele13 !=0 or ele14 !=0 or ele15 !=0 or ele16 !=0 or ele17 !=0:
      output.write(f"{ele1}|{ele2}|{ele3}|{ele4}|{ele5}|{ele6}|{ele7}|{ele8}|{ele9}|{ele10}|{ele11}|{ele12}|{ele13}|{ele14}|{ele15}|{ele16}|{ele17}|{ele18} \n")   #write data into file
output.close()
data=pd.read_csv(f"/Users/jasonhwang/Documents/pdb_trim/arpeggio_data_{pdb}{chain}_processed.csv", sep="|", header=None)
data.columns=["site 1", "site 2", "clash", "covalent", "Vdw clash", "Vdw", "proximal", "H-bond", "weak H-bond", "halogen bond", "ionic", "metal complex", "aromatic", "hydrophobic", "carbonyl", "polar", "weak polar", "interacting entities"]

data.to_csv(f"arpeggio_data_post_{pdb}_{chain}.csv")

print("appegcombo done")

#NOTE: change everything from content to proper pathway

#reads the csv file and seperates
data_new=pd.read_csv(f"/Users/jasonhwang/Documents/pdb_trim/{pdb}_{chain}.csv", sep="\s+", header=None)
#adds columns to df
data_new.columns=["dist", "vh atoms", "rbd atoms"]


#reads pdb file and creates dataframes
num=ppdb["atom num"]
resname=ppdb["residue"]
resnum=ppdb["residue num"]
#combines dataframes and gives it a column title
pdbdata=pd.DataFrame({"atom_num":num, "res id":resnum, "residue":resname})
#goes through pdb data and inputs residue and atom data given atom id
app_rbdcode=[]
app_vhcode=[]
res=[]
res_rbd=[]
ident=[]

for index, row in data_new.iterrows():
    ident.append(ppdb.iloc[int(row["vh atoms"])-1, 0])
    app_vhcode.append(ppdb.iloc[int(row["vh atoms"])-1, 3])
    res.append(pdbdata.iloc[int(row["vh atoms"])-1, 1])
    #print(ppdb.iloc[int(row["vh atoms"])-1, 0])
    x=ppdb.iloc[int(row["vh atoms"])-1, 0]
      
for index, row in data_new.iterrows():
    app_rbdcode.append(ppdb.iloc[int(row["rbd atoms"])-1, 3])
    res_rbd.append(ppdb.iloc[int(row["rbd atoms"])-1, 4])


coast=pd.DataFrame({"vh atom":ident, "vh residue": app_vhcode, "res location":res, "rbd res":app_rbdcode, "rbd loc":res_rbd})

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


for i in range(len(app_vhcode)):
  app_vhcode[i]=aa_code[app_vhcode[i]]
for j in range(len(app_rbdcode)):
  app_rbdcode[j]=aa_code[app_rbdcode[j]]

    
cooked=pd.read_csv(f"/Users/jasonhwang/Documents/pdb_trim/final_{pdb}_{chained}.csv")
cookies=pd.read_csv(f"/Users/jasonhwang/Documents/pdb_trim/arpeggio_data_post_{pdb}_{chained}.csv")

z=list(cookies.columns)
f=open(f"output_{pdb}_{chained}_arpegg_final.csv", "w")
chain=[]
for i in range(len(cooked)):
  for m in range(len(app_vhcode[i])):
    x=f"{altered_protein[0][0]}/{res[i]}/{app_vhcode[i][m]}"
    print(x)
    for n in range(len(app_rbdcode[i])):
      y=f"{altered_protein[2][0]}/{res_rbd[i]}/{app_rbdcode[i][n]}"
      print(y)
      conditions=["clash","covalent","Vdw clash",'Vdw','proximal','H-bond','weak H-bond','halogen bond','ionic','metal complex','aromatic','hydrophobic','carbonyl','polar','weak polar']
      pain=cookies[(cookies["site 1"]==x)|(cookies["site 2"]==x)]
      surgery=pain[(pain["site 1"]==y)|(pain["site 2"]==y)]
      if surgery.empty==False:
        for j in range(len(conditions)):
          g=surgery[["site 1", "site 2", conditions[j]]]
          k=surgery[conditions[j]].to_string(index=False)
          if g[conditions[j]].values==1:
            currlocation=surgery["site 1"].to_string(index=False)
            currtar=surgery["site 2"].to_string(index=False)
            chain.append(f"{currlocation}|{currtar}|{conditions[j]}\n")

pizza=list(set(chain))

for lines in pizza:
  f.write(lines)

print("organization done")
#organization done









