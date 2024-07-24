#NOTE: this is a script for pymol. This will NOT run on VS Code or in any IDEs, nor is it designed for it to run on those IDEs
#go into pymol and select "run script" then click this. 

import pymol
from pymol import cmd


#####################################

#put in the PDB code
PBD_code="v30v4_trimmed"
#input vh or vl (all lowercase)
selection="vh"

with open("currData.txt", "w+") as vader:
    vader.truncate(0)
    vader.write(f"{PBD_code} {selection}")



#note:[chain identifier, what to rename it as] add more as needed
altered_protein=[]
altered_protein.append(["C//", "vh"])
altered_protein.append(["D//", "vl"])
altered_protein.append(["A//", "wt_rbd"])

###############################
#cmd.fetch(PBD_code)
cmd.load("/Users/jasonhwang/Documents/annotation/v30v4_trimmed.pdb")
#colors the image based on chain
cmd.util.cbc()

#####################################
currProtein=[]

#NOTE: following code is used for labeling and idenfitying the stuff we want to examine->VERY IMPORTANT

#this gets list of every chain
for ch in cmd.get_chains():
    currProtein.append(ch)
#renames the chain
for i in range(len(altered_protein)):
    print(f"/{PBD_code}/{altered_protein[i][0]}/")
    cmd.alter(f'/{PBD_code}/{altered_protein[i][0]}', f'chain="{altered_protein[i][1]}"')
#select
for i in range(len(altered_protein)):
    cmd.select(f"{altered_protein[i][1]}", f"chain {altered_protein[i][1]}")
# epitotes
for i in range(len(altered_protein)):
    if altered_protein[i][1] !="wt_rbd":
        cmd.select("epi "+altered_protein[i][1], f"chain wt_rbd within 5.0 of chain {altered_protein[i][1]}")
# paratope
for i in range(len(altered_protein)):
    if altered_protein[i][1] != "wt_rbd":
        cmd.select("para "+altered_protein[i][1], f"chain {altered_protein[i][1]} within 5.0 of chain wt_rbd")

#get atom id of atoms in specified selection
epi_vh_p=cmd.index("epi_vh")
epi_vl_p=cmd.index("epi_vl")
para_vh_p=cmd.index("para_vh")
para_vl_p=cmd.index("para_vl")

#####################################
#be sure to check if checking for vl or vh

#in short, code:
#creates new file->for every aa in specfied chain and in rbd, get distance and write into file

def sele_exist(sele):
    return cmd.select(f"id {sele}")
print(para_vh_p)
#nodo=[4999, 5000, 5001, 5002, 5012, 5013, 5014, 5015, 5016, 5017, 5018, 5023, 5031, 5032, 5072]
with open(f"{PBD_code}_{selection}.txt", "w+") as kenobi:
    kenobi.truncate(0)
    if selection=="vl": 
        for o in range(len(epi_vl_p)):
            pizza=sele_exist(epi_vl_p[o][1])
            if pizza==1:
                l1=cmd.select("l1", f"id {epi_vl_p[o][1]}")
                for p in range(len(para_vl_p)):
                    tie_fighter=sele_exist(para_vl_p[p][1])
                    if tie_fighter==1:
                        print(para_vl_p[p][1])
                        l2=cmd.select("l2", f"id {para_vl_p[p][1]}")
                        dur="N/A"
                        print(l1)
                        dur=cmd.get_distance("l1", "l2")
                        kenobi.write(f"{dur} {para_vl_p[p][1]} {epi_vl_p[o][1]}\n")
                        print("completed")
    if selection=="vh":
        for i in range(len(epi_vh_p)): #sel1
            #pizza=cmd.select(f"id {epi_vh_p[i][1]}")
            cake=sele_exist(epi_vh_p[i][1])
            if cake==1:
                p1=cmd.select("p1", f"id {epi_vh_p[i][1]}")
                for j in range(len(para_vh_p)): #sel2
 #                   if para_vh_p[j][1] in nodo:
  #                      continue
                    pie=sele_exist(para_vh_p[j][1])
                    if pie==1:
                        p2=cmd.select("p2", f"id {para_vh_p[j][1]}")
                        print(para_vh_p[j][1])
                        dist="N/A"
                        dist=cmd.get_distance("p1", "p2")
                        kenobi.write(f"{dist} {para_vh_p[j][1]} {epi_vh_p[i][1]} \n" )





# import pyrotein as pr
# import os
# import pymolPy3


#put in the PDB code
# temp_pdb="6xc4.pdb"
# temp_list=pr.atom.read(temp_pdb)
# temp_dict=pr.atom.create_lookup_table(temp_list)
# print (temp_dict)
# pm=pymolPy3.pymolPy3(0)
# pm(f"fetch 6xc4")
# altered_protein=[]

# g=open("vl_7chb1.txt", "w")
# for i in range(len(epi_vl_p)):
#     for j in range(len(para_vl_p)):
#         x=cmd.get_distance(f"id {epi_vl_p[i][1]}", f"id {para_vl_p[j][1]}")
#         g.write(f"{x} {para_vl_p[j][1]} {epi_vl_p[i][1]} \n" )

# f=open(f"7chb_vh.txt", "w") 
# for i in range(len(epi_vh_p)): #sel1
#     p1=cmd.select("p1", f"id {epi_vh_p[i][1]}")

#     for j in range(len(para_vh_p)): #sel2
#         p2=cmd.select("p2", f"id {para_vh_p[j][1]}")
#         dist="N/A"
#         if p2==0 or p1==0:
#             break
#         else:
#             dist=cmd.get_distance("p1", "p2")
#             f.write(f"{dist} {para_vh_p[j][1]} {epi_vh_p[i][1]} \n" )




#if len(epi_vl_p)<len(para_vl_p):
    #for j in range(len(para_vl_p)): #sel2
        #p1=cmd.select("p1", f"id {epi_vl_p[i][1]}")

        #for i in range(len(epi_vl_p)): #sel1
       #     p2=cmd.select("p2", f"id {para_vl_p[j][1]}")
      #      print(para_vl_p[j][1])
     #       dist="N/A"
    #        if p1==1 and p2==1:
   #             dist=cmd.get_distance("p1", "p2")
  #          f.write(f"{dist} {para_vl_p[j][1]} {epi_vl_p[i][1]} \n" )
 #           if p2==0:
#                break

# g=open("7chb1234.txt", "w")
# for epi in range(len(epi_vl_p)):
#     for epivh in range(len(epi_vh_p)):
# #     for j in range(len(para_vh_p)):
# #         x=cmd.get_distance(f"id {epi_vh_p[i][1]}", f"id {para_vh_p[j][1]}")
# #         f.write(f"{x} {para_vh_p[j][1]} {epi_vh_p[i][1]} \n" )
#         for j in range(len(para_vl_p)):
#             x=cmd.get_distance(f"id {epi_vl_p[i][1]}", f"id {para_vl_p[j][1]}")
    
#####################################
#NOTE: below. this is one way to determine interactions. does NOT return list tho

# paratope
# for i in range(len(altered_protein)):
#     if altered_protein[i][1] != "wt_rbd":
#         cmd.select(f"wt+{altered_protein[i][1]}", f"chain {altered_protein[i][1]}+wt_rbd")
# for i in range(len(altered_protein)):
#     if altered_protein[i][1] != "wt_rbd":
#         x=cmd.util.interchain_distances(f"wt+{altered_protein[i][1]}_interchain_any", f"wt+{altered_protein[i][1]}", cutoff=5.0)
#         # f.write(x)