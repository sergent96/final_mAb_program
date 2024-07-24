#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:25:15 2024

@author: jasonhwang
"""

#NOTE: this is a script for pymol. This will NOT run on VS Code or in any IDEs, nor is it designed for it to run on those IDEs
#go into pymol and select "run script" then click this. 
    
import pymol
from pymol import cmd


#####################################

#put in the PDB code
PBD_code="7ch5"
#input vh or vl (all lowercase)
selection="vh"
name_now="fold_bd57_0394_wt_model_0"

name=f"/Users/jasonhwang/Documents/MAb_Af3/{name_now}.cif"

with open("currData.txt", "w+") as vader:
    vader.truncate(0)
    vader.write(f"{PBD_code} {selection} {name_now} {name}")


#note:[chain identifier, what to rename it as] add more as needed
altered_protein=[]
altered_protein.append(["A", "vh"])
altered_protein.append(["B", "vl"])
altered_protein.append(["C", "wt_rbd"])

###############################
cmd.load(name)

#colors the image based on chain
cmd.util.cbc()

#####################################
currProtein=[]

#NOTE: following code is used for labeling and idenfitying the stuff we want to examine->VERY IMPORTANT

#this gets list of every chain
for ch in cmd.get_chains():
    currProtein.append(ch)
#select
for i in range(len(altered_protein)):
    cmd.select(f"{altered_protein[i][1]}", f"chain {altered_protein[i][1]}")
# epitotes
for i in range(len(altered_protein)):
    if altered_protein[i][1] !="wt_rbd":
        cmd.select("epi "+altered_protein[i][1], f"chain C within 5.0 of chain {altered_protein[i][0]}")
# paratope
for i in range(len(altered_protein)):
    if altered_protein[i][1] != "wt_rbd":
        cmd.select("para "+altered_protein[i][1], f"chain {altered_protein[i][0]} within 5.0 of chain C")

#get atom id of atoms in specified selection
epi_vh_p=cmd.index("epi_vh")
epi_vl_p=cmd.index("epi_vl")

para_vh_p=cmd.index("para_vh")
para_vl_p=cmd.index("para_vl")

#####################################
#NOTE: adjust as see fit. code will result in error for some reason, and can only run one at a time, 
#be sure to check if checking for vl or vh

#in short, code:
#creates new file->for every aa in specfied chain and in rbd, get distance and write into file

def sele_exist(sele):
    return cmd.select(f"id {sele}")

with open(f"{name_now}_{selection}.txt", "w+") as kenobi:
    kenobi.truncate(0)
    if selection=="vl": 
        for o in range(len(epi_vl_p)):
            pizza=sele_exist(epi_vl_p[o][1])
            if pizza==1:
                l1=cmd.select("l1", f"id {epi_vl_p[o][1]}")
                for p in range(len(para_vl_p)):
                    tie_fighter=sele_exist(para_vl_p[p][1])
                    if tie_fighter==1:
                        l2=cmd.select("l2", f"id {para_vl_p[p][1]}")
                        dur="N/A"
                        print(para_vl_p[p][1])
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
                    pie=sele_exist(para_vh_p[j][1])
                    if pie==1:
                        p2=cmd.select("p2", f"id {para_vh_p[j][1]}")

                        dist="N/A"
                        dist=cmd.get_distance("p1", "p2")
                        print(dist)
                        para=para_vh_p[j][1]
                        epi=epi_vh_p[i][1]
                        print(epi)

                        kenobi.write(f"{dist} {para} {epi} \n")
                        #kenobi.write(f"{dist} {para_vh_p[j][1]} {epi_vh_p[i][1]} \n" )      
print(name_now)
cmd.reinitialize()

