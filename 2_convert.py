#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:36:28 2024

@author: jasonhwang
"""

import os

g=open(f"currData.txt").read().split()
print(g)
output=open(f"{g[2]}_{g[1]}.csv", "w")

#open file, read, split
text_list=open(f"/Users/jasonhwang/Documents/MAb_Af3/{g[2]}_{g[1]}.txt").read().split()
#group into 3, conditions for split
grouped_text_list=list(zip(*(iter(text_list),)*3))[1:] 
for ele1, ele2, ele3 in grouped_text_list: 
    if float(ele1)  <= 5.0:  #if less than 5, area of interest
        output.write(f"{ele1}|{ele2}|{ele3}\n")   #write data into file
print("completed")