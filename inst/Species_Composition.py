#!/usr/bin/env python3
"""
Arguments:
(1) path_to_align: Absolute path to alignment file
(2) spp_info: Semicolon-separated list of taxa to subset from alignment
(3) align_name: Label for dataset/subset
"""
import sys
import os
import inspect

# Fix path for reticulate
src_file_path = inspect.getfile(lambda: None)

if src_file_path not in sys.path:
    sys.path.append(src_file_path)
    
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from readAlignment import getPrunedAlignment
    
def fetchSpeciesComposition(path_to_align,spp_info,align_name):
    
    # Set path to alignment
    alignment_path = str(path_to_align)
    
    # Get species list from semicolon-separated string
    spp_list = sorted(str(spp_info).split(";"))

    # Set alignment name
    alignment_name = str(align_name)

    # Read in alignment and prune to desired species if requested
    try:
        pruned_alignment = getPrunedAlignment(alignment_path,spp_list)        
    except:
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

if __name__ == "__main__":
    fetchSpeciesComposition(sys.argv[1],sys.argv[2],sys.argv[3])