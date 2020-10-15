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
from siteCounter import countAs,countCs,countGs,countTs,countNs,countGaps
from readAlignment import getPrunedAlignment

def fetchAlignmentComposition(path_to_align,spp_info,align_name):
    
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
    
    # Count A, C, T, G, N, -
    a_count = countAs(pruned_alignment)
    c_count = countCs(pruned_alignment)
    g_count = countGs(pruned_alignment)
    t_count = countTs(pruned_alignment)
    n_count = countNs(pruned_alignment)
    gap_count = countGaps(pruned_alignment)

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
    
    return_df = pd.DataFrame([[alignment_name,alignment_length,percent_gc,percent_n,percent_gap]],columns=['Alignment_Name','Alignment_Length','Percent_GC','Percent_N','Percent_Gap'])
    return(return_df)

if __name__ == "__main__":
    fetchAlignmentComposition(sys.argv[1],sys.argv[2],sys.argv[3])