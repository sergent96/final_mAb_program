#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:34:27 2024

@author: jasonhwang
"""

#NOTE: this is a script for pymol. This will NOT run on VS Code or in any IDEs, nor is it designed for it to run on those IDEs
#go into pymol and select "run script" then click this. 

#this will return trimmed PDB and CIF file

import pymol
from pymol import cmd

#####################################

#put in the PDB code

g=open("currData_selection.txt").read().split()

PBD_code=str(g[0])
selection=str(g[1])
target=str(g[2])

cmd.fetch(PBD_code)

with open("currData.txt", "w+") as vader:
    vader.truncate(0)
    vader.write(f"{PBD_code} {selection}")

#input every important chain coordinates
#A for RBD, B for vh, C for vl-> ONLY MODIFY first and last. keep middle ("save") constant
keep=[]
coods=[]
name=[]
x=open("/Users/jasonhwang/Documents/pdb_trim/currSeperation.txt").read().split()
print(len(x))
for i in range(len(x)):
    if i%2==0:
        print(x[i])
        coods.append(x[i])
    if i%2==1:
        print(x[i])
        name.append(x[i])
#print(coods)
print(coods)
print(name)
for i in range(len(coods)):
    keep.append([f"{coods[i]}/", "save", f"{name[i]}"])
print(keep)
#keep.append(["B//", "save", "A"])
#keep.append(["C//", "save", "A"])
#keep.append(["D/H/","save", "B"])
#keep.append(["E/L/", "save", "C"])



#from this point onwards, no code needs to be modified. 


cookies=[]
for i in range(len(keep)):
    x=keep[i][0].split("/", 10)
    cookies.append(keep[i][0])
 
pizza=[]
taco=[]
for i in range(len(keep)):
 #   print(i)
    x=keep[i][0].split("/", 10)
    if len(x[1])==1:
        pizza.append(x[1])
    else:
        pizza.append(x[0])
    taco.append(x[0])

#print(cookies)

for x in cmd.get_names():
    for ch in cmd.get_chains(x):
        if ch not in pizza:
            cmd.remove(f"chain {ch}")
            print(ch)
for x in cmd.get_names():
    for ch in cmd.get_chains(x):
        #print(ch)
        cmd.alter(f"/{PBD_code}/{ch}/", f'chain="delete"')
for i in range(len(keep)):
   g=f"/{PBD_code}/{keep[i][0]}"
   cmd.alter(g, f'chain="{keep[i][1]}"')
for x in cmd.get_names():
    for ch in cmd.get_chains(x):
        if ch != "save":
            cmd.remove(f"chain {ch}")

for i in range(len(keep)):
    #if taco[i][1]=="":
    g=f"/{PBD_code}/{taco[i]}/save/"
    cmd.alter(g, f'chain="{keep[i][2]}"')

cmd.util.cbc()
cmd.save(f"{target}_{PBD_code}_trimmed.pdb")
cmd.save(f"{target}_{PBD_code}_trimmed.cif")

#select
for i in range(len(keep)):
    cmd.select(f"{keep[i][2]}", f"chain {keep[i][2]}")
# epitotes
for i in range(len(keep)):
    if keep[i][2] !="A":
        cmd.select("epi "+keep[i][2], f"chain A within 5.0 of chain {keep[i][2]}")
# paratope
for i in range(len(keep)):
    if keep[i][2] != "A":
        cmd.select("para "+keep[i][2], f"chain {keep[i][2]} within 5.0 of chain A")

#get atom id of atoms in specified selection
epi_vh_p=cmd.index("epi_B")
epi_vl_p=cmd.index("epi_C")
para_vh_p=cmd.index("para_B")
para_vl_p=cmd.index("para_C")

#####################################
#be sure to check if checking for vl or vh

#in short, code:
#creates new file->for every aa in specfied chain and in rbd, get distance and write into file

def sele_exist(sele):
    return cmd.select(f"id {sele}")

with open(f"{PBD_code}_{selection}.csv", "w+") as kenobi:
    kenobi.truncate(0)
    if selection=="vl": 
        for o in range(len(epi_vl_p)):
            pizza=sele_exist(epi_vl_p[o][1])
            if pizza==1:
                l1=cmd.select("l1", f"id {epi_vl_p[o][1]}")
                for p in range(len(para_vl_p)):
                    tie_fighter=sele_exist(para_vl_p[p][1])
                    if tie_fighter==1:
                       # print(para_vl_p[p][1])
                        l2=cmd.select("l2", f"id {para_vl_p[p][1]}")
                        dur="N/A"
                        #print(l1)
                        dur=cmd.get_distance("l1", "l2")
                        kenobi.write(f"{dur} {para_vl_p[p][1]} {epi_vl_p[o][1]}\n")
                        #print("completed")
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
                        if dist<=10:
                            kenobi.write(f"{dist} {para_vh_p[j][1]} {epi_vh_p[i][1]} \n" )
print("completed")



