#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:35:51 2024

@author: jasonhwang
"""
import pandas as pd
from CifFile import ReadCif



currSele=open("currData_selection.txt").read().split()

PBD_code=str(currSele[0])
selection=str(currSele[1])
target=str(currSele[2])
shm=str(currSele[3])

file_path = f"/Users/jasonhwang/Documents/pdb_trim/{PBD_code}_fasta.txt"

with open(file_path, 'r') as file:
    lines = file.readlines()
chains = []
for line in lines:
    chains.append(line)
identify=[]
temp_seq=[]
for i in range(len(chains)):
    if i%2==0:
        identify.append(chains[i])
    else:
        temp_seq.append(chains[i])
dqf=[]
pos=[]
for i in range(len(identify)):
    words=identify[i]
    parts = words.split('|')
    chain_info = parts[2].strip()
    if shm in chain_info:
        dqf.append(chain_info)
        pos.append(i)
seq=temp_seq[pos[0]]
print(seq)
df=pd.DataFrame({"id":identify, "seq":temp_seq})

residue_vh = {}
def annote(seq_pure):
    pos=[]
    pos.append(["AAS", "fr1_end"])
    pos.append(["MSW", "fr2_start"])
    pos.append(["MNW", "fr2_start"])
    pos.append(["WVS", "fr2_end"])
    pos.append(["YYA", "fr3_start"])
    pos.append(["FYA", "fr3_start"])
    pos.append(["YYC", "fr3_end"])
    
    cookies=[]
    position=[]
    
    cut = seq.rfind("YYC")
    spliced_input = seq[:cut + 3]
    
    for i in range(len(pos)):
        poss = spliced_input.find(pos[i][0])
        if poss != -1:
            cookies.append(pos[i][0])
            position.append(poss)
        
    print(position)
    
    for i in range(1, position[0]):
        residue_vh[i] = "fr1"
    for i in range(position[0], position[1]):
        residue_vh[i] = "cdr1"
    for i in range(position[1], position[2]):
        residue_vh[i] = "fr2"
    for i in range(position[2], position[3]):
        residue_vh[i] = "cdr2"
    for i in range(position[3], position[4]):
        residue_vh[i] = "fr3"
    for i in range(position[4], len(seq)):
        residue_vh[i] = "cdr3"
    return residue_vh
seq_annotated = []
def cut_set(seq):
    global seq_annotated  # Use the global variable seq_annotated
    cut = seq.rfind("YYC")
    spliced_input = seq[:cut + 3]
    spliced_input_list = [*spliced_input]
    gl = "EVQLVESGGGLIQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
    
    # Annotate sequences
    seq_annotated = ["n/a"] * len(seq)
    for i in range(len(spliced_input_list)):
        if seq[i] == gl[i]:
            seq_annotated[i] = f"{spliced_input[i]} {i + 1}"
        else:
            seq_annotated[i] = f"{gl[i]} {i + 1} {spliced_input[i]}"
    return seq_annotated

seq=temp_seq[pos[0]]
annote(seq)
cut_set(seq)

# Read input data
g = open("currData.txt").read().split()
pdb = str(g[0])
chain = str(g[1])
path=f"/Users/jasonhwang/Documents/pdb_trim/{pdb}_trimmed.cif"


# Path to your CIF file
cif_file_path = path
# Read the CIF file
cif_data = ReadCif(cif_file_path)
# Assuming the CIF file contains a single data block
data_block = list(cif_data.keys())[0]
# Extract relevant data from the CIF file
data = cif_data[data_block]
# Convert the data to a pandas DataFrame
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


altered_protein = [
    ["H", "vh"],
    ["L", "vl"],
    ["R", "wt_rbd"]
]

# Read the CSV file and separate data
data = pd.read_csv(f"/Users/jasonhwang/Documents/pdb_trim/{pdb}_{chain}.csv", sep="\s+", header=None)
data.columns = ["dist", "vh atoms", "rbd atoms"]


# Read PDB file and create dataframes
num=alpha["atom num"]
resname=alpha["residue"]
resnum=alpha["residue num"]
pdbdata = pd.DataFrame({"atom_num": num, "res id": resnum, "residue": resname})

# Process VH atoms
vh, vhcode, app_vhcode, res, dist, id = [], [], [], [], [], []
seq_print = []

for index, row in data.iterrows():
    vh.append(alpha.iloc[int(row["vh atoms"])-1, 0])
    vhcode.append(alpha.iloc[int(row["vh atoms"])-1, 3])
    y = alpha.iloc[int(row["vh atoms"])-1, 4]

    k = seq_annotated[int(y) - 1]
    seq_print.append(k)
    currPos = int(y)
    res.append(y)
    dist.append(data["dist"].iloc[index])
    z = residue_vh.get(currPos)
    id.append(z)

# Process RBD atoms
rbd, rbdcode, app_rbdcode, res_rbd = [], [], [], []
for index, row in data.iterrows():
    rbd.append(alpha.iloc[int(row["rbd atoms"])-1, 0])
    rbdcode.append(alpha.iloc[int(row["rbd atoms"])-1, 3])
    res_rbd.append(alpha.iloc[int(row["rbd atoms"])-1, 4])

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
    "dist": dist, "shm": seq_print, "vh residue": res, "vh atoms": vh,
    "vh aa code": app_vhcode, "vh code": vhcode, "res location": id,
    "rbd residue": res_rbd, "rbd atoms": rbd, "rbd aa code": app_rbdcode,
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
appeg_file_path = f"/Users/jasonhwang/Documents/pdb_trim/output_{pdb}_{chain}_arpegg_final.csv"

try:
    # Try to read the CSV file
    appeg = pd.read_csv(appeg_file_path, sep="|", header=None)
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
except pd.errors.EmptyDataError:
    # Handle the case where the file is empty
    print(f"The file {appeg_file_path} is empty.")
    appeg = pd.DataFrame(columns=["site 1", "site 2", "type of bond"])
except FileNotFoundError:
    # Handle the case where the file is missing
    print(f"The file {appeg_file_path} does not exist.")
    appeg = pd.DataFrame(columns=["site 1", "site 2", "type of bond"])


                
# Save the final output
final_out.to_csv(f"final_combined_{pdb}_{chain}_output.csv", index=False)

