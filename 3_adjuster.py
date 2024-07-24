import pandas as pd
from CifFile import ReadCif

g=open(f"currData.txt").read().split()
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
print(f"CIF file has been converted to CSV")

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
ppdb = pd.DataFrame({
    'atom num': num,
    'symbol': sym,
    'aa_id': atom_id,
    'residue': residue, 
    "residue num":residue_num
})

#print(ppdb.head())
#print("Asdfasdfasdf")
#print(ppdb.iloc[10, 0]) #goes by row, col ->0th row, 1st col=atom num 1, symbol
#returns row+1->iloc[10] gives value for 11th atom


#starts at 0, 0
#atom num=0, symbol=1, aa_id=2, residue=3, residue num=4

#columns = ppdb.columns
#print(columns)

#NOTE: ADD Heavy Chain seq found in fasta code above. put seq into clustal and find gl arrangement
seq="EVQLVESGGGLIQPGGSLRLSCAASEFIVSRNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLNLQMNSLRAEDTAVYYCARDYGDYYFDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDK"
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
data=pd.read_csv(f"/Users/jasonhwang/Documents/MAb_Af3/{name}_{chain}.csv", sep="|", header=None)
#adds columns to df
data.columns=["dist", "vh atoms", "rbd atoms"]

cut=seq.rfind("YYC")
spliced_input=seq[:cut+3]
spliced_input_list=[*spliced_input]


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


ppdb = pd.DataFrame({
    'atom num': num,
    'symbol': sym,
    'aa_id': atom_id,
    'residue': residue, 
    "residue num":residue_num
})
vh_num=[]
vh_sym=[]
vh_aa_id=[]
vh_res=[]
vh_res_num=[]
dist=[]
for index, row in data.iterrows():
    dist.append(data["dist"].iloc[index])
    vh_num.append(ppdb.iloc[int(row["vh atoms"])-1, 0])
    vh_sym.append(ppdb.iloc[int(row["vh atoms"])-1, 1])
    vh_aa_id.append(ppdb.iloc[int(row["vh atoms"])-1, 2])
    vh_res.append(ppdb.iloc[int(row["vh atoms"])-1, 3])
    vh_res_num.append(ppdb.iloc[int(row["vh atoms"])-1, 4])
    #goes by row, col ->0th row, 1st col=atom num 1, symbol
    #returns row+1->iloc[10] gives value for 11th atom
rbd_num=[]
rbd_sym=[]
rbd_aa_id=[]
rbd_res=[]
rbd_res_num=[]
for index, row in data.iterrows():
    rbd_num.append(ppdb.iloc[int(row["rbd atoms"])-1, 0])
    rbd_sym.append(ppdb.iloc[int(row["rbd atoms"])-1, 1])
    rbd_aa_id.append(ppdb.iloc[int(row["rbd atoms"])-1, 2])
    rbd_res.append(ppdb.iloc[int(row["rbd atoms"])-1, 3])
    rbd_res_num.append(ppdb.iloc[int(row["rbd atoms"])-1, 4])
    
    
#creates new df with combined data
output=pd.DataFrame({"dist":dist, "vh":vh_num, "vh residue": vh_res, "vh code":vh_sym,"vh aa code":vh_aa_id, "rbd pos":rbd_num, "rbd code":rbd_sym, "rbd aa code":rbd_aa_id, "res code":vh_res, "res location": vh_res_num})
#output.to_csv("/content/processedoutput.csv")


rbm=[]
con=[]
ace2=[]
v30=[]
p5a=[]
p22a=[]

#reads csv of contact data and creates a df with label
contact=pd.read_csv("/Users/jasonhwang/Documents/MAb_project/code/ConAlign - Sheet1.csv", header=None)
contact.columns=["rbd residue", "rbd code", "rbd pos", "posAA", "RBM", "conserved", "ACE2", "V30V4", "P5A_3A1", "P22A-1D1"]
#initially, everything is empty
for i in range(len(output)):
  rbm.append("")
  con.append("")
  ace2.append("")
  v30.append("")
  p5a.append("")
  p22a.append("")
  
  
  
#adds these lists to the curr df
final_out=pd.DataFrame({"dist":dist, "vh residue":vh_res, "vh":vh_num, "vh aa code":vh_aa_id, "vh code":vh_sym,"res location":vh_res_num, "rbd residue":rbd_res, "rbd":rbd_num, "rbd aa code":rbd_aa_id, "rbd code":rbd_sym, "rbd res loc":rbd_res_num,  "rbm":rbm, "conserved":con, "ACE2":ace2, "v30v4":v30, "p5a_3a1":p5a, "p22a-1d1":p22a})
#chekcs to see if the vh code and vh pos corresponds to a value in contact data. if so, changes to reflect that
for a, b in output.iterrows():
  for c, d in contact.iterrows():
    if b["vh"]==d["rbd pos"]:
      #if statement below works
      #if b["vh code"]!= d["rbd code"]:
       # a1="1"
        #final_out.loc[final_out["vh"]==b["vh"], "rbm"]=a1

      ###########

      if b["vh code"]==d["rbd code"]:
        a1=contact.loc[contact["rbd pos"]==b["vh"], "RBM"]
        final_out.loc[contact["rbd pos"]==b["vh"], "rbm"]=a1

        a2=contact.loc[contact["rbd pos"]==b["vh"], "conserved"]
        final_out.loc[contact["rbd pos"]==b["vh"], "conserved"]=a2

        a3=contact.loc[contact["rbd pos"]==b["vh"], "ACE2"]
        final_out.loc[contact["rbd pos"]==b["vh"], "ACE2"]=a3

        a4=contact.loc[contact["rbd pos"]==b["vh"], "V30V4"]
        final_out.loc[contact["rbd pos"]==b["vh"], "V30V4"]=a4

        a5=contact.loc[contact["rbd pos"]==b["vh"], "P5A_3A1"]
        final_out.loc[contact["rbd pos"]==b["vh"], "P5A_3A1"]=a5

        a6=contact.loc[contact["rbd pos"]==b["vh"], "P22A-1D1"]
        final_out.loc[contact["rbd pos"]==b["vh"], "P22A-1D1"]=a6

final_out.to_csv(f"final_{name}_{chain}.csv")

