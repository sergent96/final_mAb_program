#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 13:21:23 2024

@author: jasonhwang
"""
import pandas as pd

x=pd.read_csv("/Users/jasonhwang/Documents/pdb_trim/PDB Files - Sheet1.csv", header=None)
x.columns=["Antibody Name", "PDB", "IGHV Gene Usage", "IGK/LV Gene Usage", "Target", "Note 1", "Note 2", "D614G", "BA.1", "BA.2", "BA.2.75", 	"BA.5", "BQ.1.1", "XBB"]
#note: PDB=x.iloc[i, 1]
#name of abs=x.iloc[i, 0]
 
#webbrowser.open(query)
PBD_code="7nd6"
selection="vh"
target="COVOX-158"
seq_SHM="heavy"

with open("currData_selection.txt", "w+") as vader:
    vader.truncate(0)
    vader.write(f"{PBD_code} {selection} {target} {seq_SHM}")
print("a")

import urllib.request, urllib.error, urllib.parse

url = f"https://www.rcsb.org/fasta/entry/{PBD_code}/display"

response = urllib.request.urlopen(url)
webContent = response.read().decode('UTF-8')

f = open(f'{PBD_code}_fasta.txt', 'w')
f.write(webContent)
f.close


