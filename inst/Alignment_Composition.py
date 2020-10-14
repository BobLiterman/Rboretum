#!/usr/bin/env python3
"""
Arguments:
(1) path_to_align: Absolute path to alignment file
(2) spp_info: Semicolon-separated list of taxa to subset from alignment
(3) align_name: Label for dataset/subset
"""
import sys
import os
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
try:
    from getPrunedAlignment import getPrunedAlignment
except ImportError:
    from .getPrunedAlignment import getPrunedAlignment

def countAs():
    # countAs returns the count of A/a in the alignment
    global pruned_alignment
    a_total = []
    for seq in pruned_alignment:
        a_total.append(seq.seq.count('a')+seq.seq.count('A'))

    return int(sum(a_total))

def countCs():
    # countCs returns the count of C/c in the alignment
    global pruned_alignment
    c_total = []
    for seq in pruned_alignment:
        c_total.append(seq.seq.count('c')+seq.seq.count('C'))
    
    return int(sum(c_total))

def countGs():
    # countGs returns the count of G/g in the alignment
    global pruned_alignment
    g_total = []
    for seq in pruned_alignment:
        g_total.append(seq.seq.count('g')+seq.seq.count('G'))

    return int(sum(g_total))
        
def countTs():
    # countTs returns the count of T/t in the alignment
    global pruned_alignment
    t_total = []
    for seq in pruned_alignment:
        t_total.append(seq.seq.count('t')+seq.seq.count('T'))

    return int(sum(t_total))

def countNs():
    # countNs returns the count of N/n in the alignment
    global pruned_alignment
    n_total = []
    for seq in pruned_alignment:
        n_total.append(seq.seq.count('n')+seq.seq.count('N'))

    return int(sum(n_total))

def countGaps():
    # countGaps returns the count of "-" in the alignment
    global pruned_alignment
    gap_total = []
    for seq in pruned_alignment:
        gap_total.append(seq.seq.count('-'))

    return int(sum(gap_total))

def fetchAlignmentComposition(path_to_align,spp_info,align_name):
    
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
    if not getPrunedAlignment.getPrunedAlignment():
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path)+" with provided species list.")
    
    # Get alignment length
    alignment_length = int(pruned_alignment.get_alignment_length())
    
    # Count A, C, T, G, N, -
    a_count = countAs()
    c_count = countCs()
    g_count = countGs()
    t_count = countTs()
    n_count = countNs()
    gap_count = countGaps()

    # Sum ACTG + GC bases and get percent GC
    all_base_count = sum([a_count,c_count,g_count,t_count])
    gc_count = float(sum([c_count,g_count]))
    
    # Calculate Percent GC
    if all_base_count == 0:
        percent_gc = np.nan
    else:
        percent_gc = gc_count/all_base_count
    
    # Calculate Percent N + Percent Gap
    percent_n = float(n_count)/(alignment_length*len(spp_list))
    percent_gap = float(gap_count)/(alignment_length*len(spp_list))
    
    return(pd.DataFrame([[alignment_name,alignment_length,percent_gc,percent_n,percent_gap]],columns=['Alignment_Name','Alignment_Length','Percent_GC','Percent_N','Percent_Gap']))