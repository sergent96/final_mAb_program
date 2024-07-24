#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 10:12:22 2024

@author: jasonhwang
"""

import biopandas
from biopandas.pdb import PandasPdb
import pandas as pd
from CifFile import ReadCif

seq="EVQLVESGGGLIQPGGSLRLSCTASEIIVSRNYMSWVRQAPGKGLEWVSLIYAGGSTFYADSVKGRFTISRDDSKNTLYLQMNSLRAEDTAVYYCARDLFEAGATDYWGQGTLVTVSSAS"
gl= "EVQLVESGGGLIQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"

# Read input data
g = open("currData.txt").read().split()
pdb = str(g[0])
ppdb = PandasPdb().fetch_pdb(pdb)
chain = str(g[1])
path=f"/Users/jasonhwang/Documents/annotation/{pdb}.cif"


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
alpha = pd.DataFrame({
    'atom num': num,
    'symbol': sym,
    'aa_id': atom_id,
    'residue': residue, 
    "residue num":residue_num
})

print(alpha)

altered_protein = [
    ["H", "vh"],
    ["L", "vl"],
    ["R", "wt_rbd"]
]

hc_timeline = [26, 34, 51, 58, 96, 98, 118]

# Annotate residues
residue_vh = {}
for i in range(1, hc_timeline[0]):
    residue_vh[i] = "fr1"
for i in range(hc_timeline[0], hc_timeline[1]):
    residue_vh[i] = "cdr1"
for i in range(hc_timeline[1], hc_timeline[2]):
    residue_vh[i] = "fr2"
for i in range(hc_timeline[2], hc_timeline[3]):
    residue_vh[i] = "cdr2"
for i in range(hc_timeline[3], hc_timeline[4]):
    residue_vh[i] = "fr3"
for i in range(hc_timeline[4], hc_timeline[5]):
    residue_vh[i] = "cdr3"
for i in range(hc_timeline[5], hc_timeline[6]):
    residue_vh[i] = "fr4"

residue_vl = {}
for i in range(1, 26):
    residue_vl[i] = "fr1"
for i in range(27, 32):
    residue_vl[i] = "cdr2"
for i in range(33, 49):
    residue_vl[i] = "fr2"
for i in range(50, 52):
    residue_vl[i] = "cdr2"
for i in range(53, 88):
    residue_vl[i] = "fr3"
for i in range(89, 95):
    residue_vl[i] = "gl"
for i in range(95, 100000):
    residue_vl[i] = "vl"

# Read the CSV file and separate data
data = pd.read_csv(f"/Users/jasonhwang/Documents/annotation/{pdb}_{chain}.csv", sep="|", header=None)
data.columns = ["dist", "vh atoms", "rbd atoms"]

cut = seq.rfind("YYC")
spliced_input = seq[:cut + 3]
spliced_input_list = [*spliced_input]

# Read PDB file and create dataframes
num = ppdb.df['ATOM']["atom_number"]
resname = ppdb.df["ATOM"]["residue_name"]
resnum = ppdb.df["ATOM"]["residue_number"]
pdbdata = pd.DataFrame({"atom_num": num, "res id": resnum, "residue": resname})

# Annotate sequences
seq_annotated = ["n/a"] * len(seq)
for i in range(len(spliced_input_list)):
    if seq[i] == gl[i]:
        seq_annotated[i] = f"{spliced_input[i]} {i + 1}"
    else:
        seq_annotated[i] = f"{gl[i]} {i + 1} {spliced_input[i]}"

# Process VH atoms
vh, vhcode, app_vhcode, res, dist, id = [], [], [], [], [], []
seq_print = []

for index, row in data.iterrows():
    for i, r in pdbdata.iterrows():
        if row["vh atoms"] == r["atom_num"]:
            vh.append(row["vh atoms"])
            x = pdbdata.loc[pdbdata["atom_num"] == row["vh atoms"], "residue"]
            vhcode.append(x.to_string(index=False))
            y = pdbdata.loc[pdbdata["atom_num"] == row["vh atoms"], "res id"]
            k = seq_annotated[int(y.to_string(index=False)) - 1]
            seq_print.append(k)
            currPos = int(y.to_string(index=False))
            res.append(y.to_string(index=False))
            dist.append(data["dist"].iloc[index])
            z = residue_vh[currPos]
            id.append(z)

# Process RBD atoms
rbd, rbdcode, app_rbdcode, res_rbd = [], [], [], []
for index, row in data.iterrows():
    for i, r in pdbdata.iterrows():
        if row["rbd atoms"] == r["atom_num"]:
            rbd.append(row["rbd atoms"])
            x = pdbdata.loc[pdbdata["atom_num"] == row["rbd atoms"], "residue"]
            rbdcode.append(x.to_string(index=False))
            y = pdbdata.loc[pdbdata["atom_num"] == row["rbd atoms"], "res id"]
            res_rbd.append(y.to_string(index=False))

# Map amino acid codes
amino_acids = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


# Assign amino acid codes
vhcode = [amino_acids[v] for v in vhcode]
rbdcode = [amino_acids[r] for r in rbdcode]

for index, row in data.iterrows():
    app_vhcode.append(alpha.iloc[int(row["vh atoms"])-1, 2])
    app_rbdcode.append(alpha.iloc[int(row["rbd atoms"])-1, 2])

# Create the output DataFrame
output = pd.DataFrame({
    "dist": dist, "shm": seq_print, "vh": vh, "vh code": vhcode,
    "vh aa code": app_vhcode, "rbd pos": rbd, "rbd code": rbdcode,
    "rbd aa code": app_rbdcode, "res code": res, "res location": id
})

# Read contact data
contact = pd.read_csv("/Users/jasonhwang/Documents/MAb_project/code/ConAlign - Sheet1.csv", header=None)
contact.columns = ["rbd residue", "rbd code", "rbd pos", "posAA", "RBM", "conserved", "ACE2", "V30V4", "P5A_3A1", "P22A-1D1"]

# Initialize columns in the final DataFrame
target = [""] * len(output)
curr = [""] * len(output)

typeBond = [""] * len(output)
rbm = [""] * len(output)
con = [""] * len(output)
ace2 = [""] * len(output)
v30 = [""] * len(output)
p5a = [""] * len(output)
p22a = [""] * len(output)

final_out = pd.DataFrame({
    "dist": dist, "shm": seq_print, "vh residue": res, "vh": vh,
    "vh aa code": app_vhcode, "vh code": vhcode, "res location": id,
    "rbd residue": res_rbd, "rbd": rbd, "rbd aa code": app_rbdcode,
    "rbd code": rbdcode, "site 1":curr, "site 2": target, "bond": typeBond, "rbm": rbm,
    "conserved": con, "ACE2": ace2, "v30v4": v30, "p5a_3a1": p5a, "p22a-1d1": p22a
})

# Update final DataFrame based on contact data
for a, b in output.iterrows():
    for c, d in contact.iterrows():
        if b["vh"] == d["rbd pos"] and b["vh code"] == d["rbd code"]:
            final_out.loc[contact["rbd pos"] == b["vh"], "rbm"] = contact.loc[contact["rbd pos"] == b["vh"], "RBM"]
            final_out.loc[contact["rbd pos"] == b["vh"], "conserved"] = contact.loc[contact["rbd pos"] == b["vh"], "conserved"]
            final_out.loc[contact["rbd pos"] == b["vh"], "ACE2"] = contact.loc[contact["rbd pos"] == b["vh"], "ACE2"]
            final_out.loc[contact["rbd pos"] == b["vh"], "v30v4"] = contact.loc[contact["rbd pos"] == b["vh"], "V30V4"]
            final_out.loc[contact["rbd pos"] == b["vh"], "p5a_3a1"] = contact.loc[contact["rbd pos"] == b["vh"], "P5A_3A1"]
            final_out.loc[contact["rbd pos"] == b["vh"], "p22a-1d1"] = contact.loc[contact["rbd pos"] == b["vh"], "P22A-1D1"]

appeg = pd.read_csv(f"/Users/jasonhwang/Documents/annotation/output_{pdb}_{chain}_arpegg_final.csv", sep="|", header=None)
appeg.columns = ["site 1", "site 2", "type of bond"]

if len(app_vhcode) != len(final_out):
    print("Length mismatch between app_vhcode and final_out.")
else:
        # Update final DataFrame based on ARPEGG data
    for i in range(len(final_out)):
        g = f"{altered_protein[0][0]}/{res[i]}/{app_vhcode[i]}"
        l = f"{altered_protein[2][0]}/{res_rbd[i]}/{app_rbdcode[i]}"
        
        # Filter rows where either (g, l) or (l, g) matches in 'site 1' and 'site 2'
        matched_rows = appeg[((appeg["site 1"] == g) & (appeg["site 2"] == l)) | ((appeg["site 1"] == l) & (appeg["site 2"] == g))]
        
        if not matched_rows.empty:
            curr = matched_rows["site 1"].values[0] if len(matched_rows) > 0 else None
            final_out.at[i, "site 1"] = curr
            
            target = matched_rows["site 2"].values[0] if len(matched_rows) > 0 else None
            final_out.at[i, "site 2"] = target
            
            typeBond = matched_rows["type of bond"].values[0] if len(matched_rows) > 0 else None
            final_out.at[i, "bond"] = typeBond
        else:
            # Handle case where no matches found (optional)
            final_out.at[i, "site 1"] = ""
            final_out.at[i, "site 2"] = ""
            final_out.at[i, "bond"] = ""
                
# Save the final output
final_out.to_csv(f"final_combined_{pdb}_{chain}_output.csv", index=False)