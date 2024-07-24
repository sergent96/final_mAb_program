#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:52:06 2024

@author: jasonhwang
"""

import pandas as pd
from Bio import Align
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans



ic50=pd.read_csv("/Users/jasonhwang/Documents/heat_map/ic50_cov_Val.csv")
abs_name=ic50["Antibody  Name"]
gene_vl=ic50["Light chain V gene"]
gene_vh=ic50["Heavy chain V gene"]
sources=ic50["source"]
print(len(sources))
sorc=[]
for i in range(len(sources)):
    sorc.append("6")
for i in range(len(sources)):
    if "BA.1 convalescents" in sources[i]:
        sorc[i]="0"
    if "BA.2 convalescents" in sources[i]:
        sorc[i]="1"
    if "BA.5 convalescents" in sources[i]:
        sorc[i]="2"
    if "WT convalescents" in sources[i]:
        sorc[i]="3"
    if "SARS convalescents" in sources[i]:
        sorc[i]="4"
    if "WT vaccinees" in sources[i]:
        sorc[i]="5"

print(len(sorc))
print(sorc)

df=pd.DataFrame({"gene":gene_vh, "data":sorc})
print(df)
x=sorc



name_ascii=[]
for i in range(len(abs_name)):
    count=0
    for j in range(len(abs_name[i])):
        z=ord(abs_name[i][j])
        count+=z
    name_ascii.append(count)
name_ascii_max=max(name_ascii)
for i in range(len(name_ascii)):
    name_ascii[i]=name_ascii[i]/name_ascii_max*10

from sklearn.cluster import KMeans

name_ascii=[]
for i in range(len(abs_name)):
    count=0
    for j in range(len(abs_name[i])):
        z=ord(abs_name[i][j])
        count+=z
    name_ascii.append(count)
name_ascii_max=max(name_ascii)
for i in range(len(name_ascii)):
    name_ascii[i]=name_ascii[i]/name_ascii_max*10

#plt.title('Elbow method')
#plt.xlabel('Number of clusters')
#plt.ylabel('Inertia')

#kmeans = KMeans(n_clusters=5)
#kmeans.fit(data)
#plt.scatter(x, y, c=kmeans.labels_)
#plt.show()

#cluster_label=kmeans.labels_

BG=ic50["D614G"]
for i in range(len(BG)):
    if BG[i]==">10":
        BG[i]="10"
BA1=ic50["BA.1"]
for i in range(len(BA1)):
    if BA1[i]==">10":
        BA1[i]="10"
BA2=ic50["BA.2"]
for i in range(len(BA2)):
    if BA2[i]==">10":
        BA2[i]="10"
BA5=ic50["BA.5"]
for i in range(len(BA5)):
    if BA5[i]==">10":
        BA5[i]="10"
BA275=ic50["BA.2.75"]
for i in range(len(BA275)):
    if BA275[i]==">10":
        BA275[i]="10"
BA2752=ic50["BA.2.75.2"]
for i in range(len(BA2752)):
    if BA2752[i]==">10":
        BA2752[i]="10"
CA1=ic50["CA.1"]
for i in range(len(CA1)):
    if CA1[i]==">10":
        CA1[i]="10"
BQ11=ic50["BQ.1.1"]
for i in range(len(BQ11)):
    if BQ11[i]==">10":
        BQ11[i]="10"
BR2=ic50["BR.2"]
for i in range(len(BR2)):
    if BR2[i]==">10":
        BR2[i]="10"
BM=ic50["BM.1.1.1"]
for i in range(len(BM)):
    if BM[i]==">10":
        BM[i]="10"
XBB=ic50["XBB"]
for i in range(len(XBB)):
    if XBB[i]==">10":
        XBB[i]="10"
cluster_data=pd.DataFrame({
    "name":abs_name, 
    "name_ascii code":name_ascii, 
    "vh seq":gene_vh, 
    "cluster":x, 
    "D614G":BG, 
    "BA.1":BA1, 
    "BA.2":BA2, 
    "BA.5":BA5, 
    "BA.2.75":BA275, 
    "BA.2.75.2":BA2752, 
    "CA.1":CA1, 
    "BQ.1.1":BQ11, 
    "BR.2":BR2, 
    "BM.1.1.1":BM, 
    "XBB":XBB
    })

sorted_cluster=cluster_data.sort_values(by="cluster", ascending=True)
sorted_cluster.to_csv("ic50_clustered_data_vh.csv")


