#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:10:04 2024

@author: jasonhwang
"""

import biopandas
import pandas as pd
from biopandas.pdb import PandasPdb

g=open("currData.txt").read().split()

pdb=str(g[0])
ppdb = PandasPdb().fetch_pdb(pdb)
chain=str(g[1])

altered_protein=[]
altered_protein.append(["A", "vh"])
altered_protein.append(["B", "vl"])
altered_protein.append(["E", "wt_rbd"])



#NOTE: ADD Heavy Chain seq found in fasta code above. put seq into clustal and find gl arrangement
seq="EVQLVESGGGLIQPGGSLRLSCAASGFIVSSNYMSWVRQAPGKGLEWVSIIYSGGSTFYADSVKGRFTISRDNSKNTLYLQMNSLRVEDTAVYYCARDLQELGSLDYWGQGTLVTVSSASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRVESKYGPPCPPCP"
gl= "EVQLVESGGGLIQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"


hc_timeline=[]
hc_timeline.append(26)
hc_timeline.append(34)
hc_timeline.append(51)
hc_timeline.append(58)
hc_timeline.append(96)
hc_timeline.append(98)
hc_timeline.append(118)
########################
#dict of residue locations
#NOTE: have the end "maker" for each range 1 larger then it. I.e. if it is from residues 1-25, the range would be (1, 26). The next range would be (26- end+1)
#this is because range of 1-10 will give us numbers 1-9 instead of 1-10. Hence, add 1 more to make it frm 1-10 (ignores 11)
#check ig blast. after the last annotated area, put the end range to some high, unreachable value
residue_vh={}
for i in range(1, hc_timeline[0]):
  residue_vh[i]="fr1"
for i in range(hc_timeline[0], hc_timeline[1]):
  residue_vh[i]="cdr1"
for i in range(hc_timeline[1], hc_timeline[2]):
  residue_vh[i]="fr2"
for i in range(hc_timeline[2], hc_timeline[3]):
  residue_vh[i]="cdr2"
for i in range(hc_timeline[3], hc_timeline[4]):
  residue_vh[i]="fr3"
for i in range(hc_timeline[4], hc_timeline[5]):
  residue_vh[i]="cdr3"
for i in range(hc_timeline[5], hc_timeline[6]):
  residue_vh[i]="fr4"
  
residue_vl={}
for i in range(1, 26):
  residue_vl[i]="fr1"
for i in range(27, 32):
  residue_vl[i]="cdr2"
for i in range(33, 49):
  residue_vl[i]="fr2"
for i in range(50, 52):
  residue_vl[i]="cdr2"
for i in range(53, 88):
  residue_vl[i]="fr3"
for i in range(89, 95):
  residue_vl[i]="gl"
for i in range(95, 100000):
  residue_vl[i]="vl"

#reads the csv file and seperates
data=pd.read_csv(f"/Users/jasonhwang/Documents/annotation/{pdb}_{chain}.csv", sep="|", header=None)
#adds columns to df
data.columns=["dist", "vh atoms", "rbd atoms"]

cut=seq.rfind("YYC")
spliced_input=seq[:cut+3]
spliced_input_list=[*spliced_input]


#reads pdb file and creates dataframes
num=ppdb.df['ATOM']["atom_number"]
resname=ppdb.df["ATOM"]["residue_name"]
resnum=ppdb.df["ATOM"]["residue_number"]
#combines dataframes and gives it a column title
pdbdata=pd.DataFrame({"atom_num":num, "res id":resnum, "residue":resname})


#just seeing if everything works correctly
#pdbdata.to_csv("/content/pdb.csv")
#data.to_csv("/content/raw.csv")


#NOTE: This creates intial df
vh=[]
vhcode=[]
app_vhcode=[]
res=[]
dist=[]
id=[]

seq_annotated=[]
for j in range(len(seq)):
  seq_annotated.append("n/a")
for i in range(len(spliced_input_list)):
    if seq[i]==gl[i] or seq[i]==gl[i]:
        seq_annotated[i]=(f"{spliced_input[i]} {i+1}")
    else:
        seq_annotated[i]=(f"{gl[i]} {i+1} {spliced_input[i]}")
seq_print=[]


#goes through pdb data and inputs residue and atom data given atom id
rbd=[]
rbdcode=[]
app_rbdcode=[]
res_rbd=[]
for index, row in data.iterrows():
  for i, r in pdbdata.iterrows():
    if row["vh atoms"]==r["atom_num"]:
      vh.append(row["vh atoms"])
      
      x=pdbdata.loc[pdbdata["atom_num"]==row["vh atoms"],"residue"]
      vhcode.append(x.to_string(index=False))
      app_vhcode.append(x.to_string(index=False))

      y=pdbdata.loc[pdbdata["atom_num"]==row["vh atoms"], "res id"]
      currPos=int(y.to_string(index=False))
      res.append(y.to_string(index=False))
      
      k=seq_annotated[int(y.to_string(index=False))-1]
      seq_print.append(k)

      dist.append(data["dist"].iloc[index])
      z=residue_vh[currPos]
      id.append(z)
      
for a, b in data.iterrows():
  for c, d in pdbdata.iterrows():
    if b["rbd atoms"]==d["atom_num"]:
      rbd.append(b["rbd atoms"])

      x=pdbdata.loc[pdbdata["atom_num"]==b["rbd atoms"],"residue"]
      rbdcode.append(x.to_string(index=False))
      app_rbdcode.append(x.to_string(index=False))

      y=pdbdata.loc[pdbdata["atom_num"]==b["rbd atoms"], "res id"]
      res_rbd.append(y.to_string(index=False))
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
#for k in range(len(res)):
 # print(res[k])
#  shmcode[i]=seq[int(res[k])]

for i in range(len(app_vhcode)):
  app_vhcode[i]=aa_code[app_vhcode[i]]
for j in range(len(app_rbdcode)):
  app_rbdcode[j]=aa_code[app_rbdcode[j]]

#creates new df with combined data
output=pd.DataFrame({"vh aa code":app_vhcode, "res code":res, "res location": id})
appegio_input=[]
for i in range(len(output)):
  if chain==altered_protein[0][1]: #vh
    k=int(res[i])+1
    for j in range(len(app_vhcode[i])):
      appegio_input.append(f"/{altered_protein[0][0]}/{k}/{app_vhcode[i][j]}")
  if chain==altered_protein[1][1]: #vh
    for j in range(len(app_vhcode[i])):
      appegio_input.append(f"/{altered_protein[1][0]}/{k}/{app_vhcode[i][j]}")
for i in range(len(res_rbd)):
  l=int(res_rbd[i])
  appegio_input.append(f"/{altered_protein[2][0]}/{l}/")
  
output.to_csv(f"final_{pdb}_{chain}.csv")

file = open(f'appegio_output_{pdb}_{chain}.txt','w')

no_duplicates_appegio = list(set(appegio_input))

for line in no_duplicates_appegio:
	file.write(line+"\n")
file.close()

print("program completed")