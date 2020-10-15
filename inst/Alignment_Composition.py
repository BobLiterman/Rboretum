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

# countAs returns the count of A/a in the alignment
def countAs(pruned_alignment):
    a_total = []
    for seq in pruned_alignment:
        a_total.append(seq.seq.count('a')+seq.seq.count('A'))

    return int(sum(a_total))

# countCs returns the count of C/c in the alignment
def countCs(pruned_alignment):
    c_total = []
    for seq in pruned_alignment:
        c_total.append(seq.seq.count('c')+seq.seq.count('C'))
    
    return int(sum(c_total))

# countGs returns the count of G/g in the alignment
def countGs(pruned_alignment):
    g_total = []
    for seq in pruned_alignment:
        g_total.append(seq.seq.count('g')+seq.seq.count('G'))

    return int(sum(g_total))

# countTs returns the count of T/t in the alignment        
def countTs(pruned_alignment):
    t_total = []
    for seq in pruned_alignment:
        t_total.append(seq.seq.count('t')+seq.seq.count('T'))

    return int(sum(t_total))

# countNs returns the count of N/n in the alignment
def countNs(pruned_alignment):
    n_total = []
    for seq in pruned_alignment:
        n_total.append(seq.seq.count('n')+seq.seq.count('N'))

    return int(sum(n_total))

# countGaps returns the count of "-" in the alignment
def countGaps(pruned_alignment):
    gap_total = []
    for seq in pruned_alignment:
        gap_total.append(seq.seq.count('-'))

    return int(sum(gap_total))

def fetchAlignmentComposition(path_to_align,spp_info,align_name):
    
    # Set path to alignment
    alignment_path = str(path_to_align)
    
    # Get species list from semicolon-separated string
    spp_list = sorted(str(spp_info).split(";"))

    # Set alignment name
    alignment_name = str(align_name)

    # Read in alignment and prune to desired species if requested
    try:
        formats = {'nex': 'nexus', 'nexus': 'nexus',
                   'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
                   'fa': 'fasta', 'fasta': 'fasta'}
        
        fformat = formats[alignment_path.split('.')[-1]]
        raw_alignment = AlignIO.read(alignment_path, fformat)

    # If alignment cannot be read in, raise exception
    except:
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path))

    # Get species from raw alignment
    raw_spp = list()
    for seq_record in raw_alignment:
        raw_spp.append(str(seq_record.id))
    
    # Convert numeric IDs to string prior to sorting
    raw_spp = [str(i) for i in raw_spp] 
    raw_spp.sort()
    
    # If fewer than three species exist in alignment or species list, raise exception
    if (len(list(set(spp_list))) < 3) or (len(list(set(raw_spp))) < 3):
        sys.exit("ERROR: Cannot process fewer than 3 species...")
        
    # If requested species are not in alignment, raise exception
    spp_diff = list(set(spp_list) - set(raw_spp))

    if len(spp_diff) > 0:
        sys.exit("ERROR: Requested species not found in "+os.path.basename(alignment_path)+"...")
    
    # Re-order alignment to sorted order of species list
    else:
        # Create dummy alignment        
        pruned_alignment = raw_alignment[0:0]
        
        # Populate alignment by adding taxa sorted by taxon ID
        for i in spp_list:
            pruned_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))
        
        # If resulting alignment is empty, raise exception
        if int(pruned_alignment.get_alignment_length()) == 0:
            sys.exit("ERROR: Alignment processed, but appears to have no bases...")

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

def main(path_to_align,spp_info,align_name):
    return(fetchAlignmentComposition(path_to_align,spp_info,align_name))
    
if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])