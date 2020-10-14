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
import align_func

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
    if not Alignment_Reader.getPrunedAlignment():
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path)+" with provided species list.")
    
    # Get alignment length
    alignment_length = int(pruned_alignment.get_alignment_length())
    
    # Count A, C, T, G, N, -
    a_count = Site_Counter.countAs()
    c_count = Site_Counter.countCs()
    g_count = Site_Counter.countGs()
    t_count = Site_Counter.countTs()
    n_count = Site_Counter.countNs()
    gap_count = Site_Counter.countGaps()

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