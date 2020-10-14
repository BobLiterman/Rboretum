#!/usr/bin/env python3
"""
Arguments:
(1) path_to_align: Absolute path to alignment file
(2) align_name: Label for dataset/subset
"""
import sys
import os
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
import align_func
    
def fetchSpeciesComposition(path_to_align,spp_info,align_name):
    
    # Set path to alignment
    global alignment_path
    alignment_path = str(path_to_align)
    
    # Get species list from semicolon-separated string
    global spp_list
    spp_list = sorted(str(spp_info).split(";"))

    # Set alignment name
    global alignment_name
    alignment_name = str(align_name)

    # Read in alignment and prune to desired species if requested
    if not Alignment_Reader.getPrunedAlignment():
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path)+" with provided species list.")
        
    # Get alignment length
    alignment_length = int(pruned_alignment.get_alignment_length())
    
    # Get counts
    a_totals = []
    c_totals = []
    t_totals = []
    g_totals = []
    n_totals = []
    gap_totals = []
    
    for seq in pruned_alignment:
        a_totals.append(seq.seq.count('a')+seq.seq.count('A'))
        c_totals.append(seq.seq.count('c')+seq.seq.count('C'))
        t_totals.append(seq.seq.count('t')+seq.seq.count('T'))
        g_totals.append(seq.seq.count('g')+seq.seq.count('G'))
        n_totals.append(seq.seq.count('n')+seq.seq.count('N'))
        gap_totals.append(seq.seq.count('-'))

    # Get total base counts and gc counts
    total_base_counts = []
    gc_counts = []
    base_zip = zip(a_totals, c_totals,t_totals,g_totals)
    for a_i,c_i,t_i,g_i in base_zip:
        total_base_counts.append(a_i+c_i+t_i+g_i)
        gc_counts.append(g_i+c_i)
    
    # Get percent GC, N, and gap
    percent_gc = []
    percent_gc_zip = zip(gc_counts,total_base_counts)
    for gc_i,total_i in percent_gc_zip:
        if total_i == 0:
            percent_gc.append(np.nan)
        else:
            percent_gc.append(gc_i/float(total_i))
    
    # Get Percent N
    percent_n = []
    for n_i in n_totals:
        percent_n.append(float(n_i)/alignment_length)
    
    # Get Percent Gap
    percent_gap = []
    for gap_i in gap_totals:
        percent_gap.append(float(gap_i)/alignment_length)
    
    df = pd.DataFrame(list(zip(spp_list,total_base_counts,n_totals,gap_totals,percent_gc,percent_n,percent_gap)), 
               columns =['Taxon','Total_Bases','Total_N','Total_Gaps','Percent_GC','Percent_N','Percent_Gap'])
    df['Alignment_Name'] = alignment_name
    return df