#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 14:36:52 2024

@author: jasonhwang
"""

import pandas as pd
from Bio import Align

# Load the CSV file
comb = pd.read_csv("/Users/jasonhwang/Documents/SHM_group_breadth.csv")

# Extract sequence IDs and sequences
seq_id = comb["Antibody Name"]
seq_comb = comb["HCDR3"]
print(seq_id)
# Create a dictionary of sequences
dict_seq = {seq_id[i]: seq_comb[i] for i in range(len(seq_id))}
seq_tracker = [seq_comb[i] for i in range(len(seq_id))]

print(len(dict_seq))
# Initialize the aligner
aligner = Align.PairwiseAligner(match_score=1.0, mode="local")
aligner.mode = 'local'

# Reference sequence
seq1 = "EVQLVESGGGLIQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
seq_det = []

for seq2 in seq_tracker:
    #print(seq1)

    # Align seq2 with the reference sequence seq1
    alignments = aligner.align(seq2, seq1)
    
    # Get the best alignment (first alignment in the list)
    best_alignment = alignments[0]
    
    # Store the sequence and its alignment
    seq_det.append(seq2)
    seq_det.append(best_alignment)

# Write the results to a file
with open("seq_aligner.txt", "w") as f:
    for i in range(len(seq_det)):
        if i % 2 == 0:
            f.write(f"{seq_id[i/2]} \n")
            f.write(f"Sequence: {seq_det[i]} \n")
        else:
            f.write(f"Best Alignment: \n \n{seq_det[i]} \n")

print("Alignment process completed.")