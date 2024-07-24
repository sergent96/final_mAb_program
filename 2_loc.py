#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 13:35:15 2024

@author: jasonhwang
"""
import pandas as pd
g=open("currData_selection.txt").read().split()

PBD_code=str(g[0])
selection=str(g[1])
target=str(g[2])

chains = []

def parse_chains(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for line in lines:
        if line.startswith('>'):
            parts = line.split('|')
            chain_info = parts[1].strip()
            chains.append(chain_info)
    return chains

def identify (file_path):
    iden = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for line in lines:
        if line.startswith('>'):
            parts = line.split('|')
            chain_info = parts[2].strip()
            iden.append(chain_info)
            #chains.append(chain_info)
    return iden
position=[]
cookies=[]


def transform_chain_auth(chain_auth_str):
    temp=[]
    chain_part=chain_auth_str.split('[')[0].split()[-1]    
    if "[auth" in chain_auth_str:        
        auth_part = chain_auth_str.split('[')[1].replace('auth ', '').replace(']', '')
        transformed_str = f"{chain_part}/{auth_part}"
        chain_id.append(transformed_str)
    else:
        auth_part = chain_auth_str.split(",")
        cleaned_chains_list = [item.replace('Chains', '').strip() for item in auth_part]
        for i in range(len(cleaned_chains_list)):
            transformed_str=f"{cleaned_chains_list[i]}/"
            chain_id.append(transformed_str)
    return transformed_str
    
chain_id=[]
def combiner(chain, iden, pos):
    pie=[]
    pizza=[]
    count=0
    for i in range(len(chain)):
        chain_auth_str = chain[i]
        #chain_auth_str = chain_auth_str.split('[')[0].split()[-1]
        #print(len(chain_auth_str))
        result=transform_chain_auth(chain_auth_str)
        if len(iden)>1:
            for j in range(len(iden)):
                pizza.append(iden[j])
                pie.append(pos[i])
                count+=1
        else:
            pizza.append(result)
            pie.append(pos[i])
            count+=1
        chain_id.clear()
        
        with open("currSeperation.txt", "w+") as x:
            x.truncate(0)
            for i in range(count):
                x.write(f"{pizza[i]} {pie[i]} \n")
                print(i)
                
file_path = f"/Users/jasonhwang/Documents/pdb_trim/{PBD_code}_fasta.txt"
chains = parse_chains(file_path)
identify=identify(file_path)

for chain in chains:
    df=pd.DataFrame({"id":identify, "chains":chains})             
                
print(df)    
for i in range(len(df)):
    if target in df.iloc[i, 0]:
        if "light" in df.iloc[i, 0]:
            position.append("C")
            cookies.append(df.iloc[i, 1])
        if "heavy" in df.iloc[i, 0]:
            position.append("B")
            cookies.append(df.iloc[i, 1])
    if "Spike" in df.iloc[i, 0]:
        position.append("A")
        cookies.append(df.iloc[i, 1])
            
combiner(cookies, chain_id, position)



