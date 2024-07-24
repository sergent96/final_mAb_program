#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 17:21:43 2024

@author: jasonhwang
"""
import csv
import pandas as pd

g=open(f"currData.txt").read().split()
print(g)
chain=str(g[1])
name_good="fold_bd57_0394"


with open(f"a_{name_good}_wt.csv", "w") as f:
    writer=csv.writer(f)
    for i in range(5):
        name=f"{name_good}_wt_model_{i}"
        k=open(f"final_{name}_{chain}.csv")
        reader=csv.reader(k)
        next(reader)
        for row in reader:
            writer.writerow(row)
        
print("combiner completed")

df=pd.read_csv(f"a_{name_good}_wt.csv")

df.columns=["count", "dist", "vh_res", "vh_pos", "vh_aa", "vh_code", "res_loc", "rbd_res", "rbd_pos", "rbd_aa", "rbd_code", "rbd_loc", "rbm", "con", "ACE", "v30v4", "p5a", "p22a"]
df["vh_coord"]=df["vh_pos"].astype(str)+df["vh_aa"]
df["rbd_coord"]=df["rbd_pos"].astype(str)+df["rbd_aa"]
df["contact"]=df["vh_coord"]+" "+df["rbd_coord"]

uniq=set(df["contact"])
uniq_list=list(uniq)
length_diff=len(df)-len(uniq_list)

if length_diff > 0:
    # Option 1: Repeat elements to extend the list
    extended_list =  uniq_list + [None] * length_diff
df["uniq_val"]=extended_list

contact_counts = df["contact"].value_counts()
df["uniq_count"] = (df["uniq_val"].map(contact_counts).fillna(0).astype(int))
df["uniq_count_normalized"] = df["uniq_count"] / 5
count_08 = (df["uniq_count_normalized"] == 0.8).sum()
count_1 = (df["uniq_count_normalized"] == 1).sum()

total_count=count_08+count_1
#df = df[df["uniq_count_normalized"] >= 0.8]

with open("data.csv", "a") as l:
    l.write(f"{name_good} | {total_count} \n")

print(df.head())

df.to_csv(f"0_final_{name_good}.csv")