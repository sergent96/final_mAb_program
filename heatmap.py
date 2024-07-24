#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 14:41:01 2024

@author: jasonhwang
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

file_path = '/Users/jasonhwang/Documents/heat_map/ic50_clustered_data_vh.csv'
data = pd.read_csv(file_path)
data_sor=data
data_col=data

# Set 'Antibody Name' as the index
#data_sor.set_index('source', inplace=True)
#heatmap_data_sorc = data_sor[["D614G", "BA.1", "BA.2", "BA.2.75", "BA.5", "BQ.1.1", "XBB"]]

data_col.set_index("cluster", inplace=True)
heatmap_data_col=data_col[["D614G", "BA.1", "BA.2", "BA.2.75", "BA.5", "BQ.1.1", "XBB"]]

#plt.figure(figsize=(10, 8))
#sns.heatmap(heatmap_data_sorc, cmap="Reds", annot=True, fmt=".2f", cbar=True)
#plt.title("Heatmap of IC50 Values: source")
#plt.show()

plt.figure(figsize=(10, 10))
sns.heatmap(heatmap_data_col, cmap="Reds", annot=True, fmt=".2f", cbar=True)
plt.title("Heatmap of IC50 Values: columns")
plt.show()