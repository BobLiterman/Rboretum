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
import pickle
import tempfile

# countAs returns the count of A/a in the alignment
def countAs():
    global pruned_alignment
    a_total = []
    for seq in pruned_alignment:
        a_total.append(seq.seq.count('a')+seq.seq.count('A'))

    return int(sum(a_total))

# countCs returns the count of C/c in the alignment
def countCs():
    global pruned_alignment
    c_total = []
    for seq in pruned_alignment:
        c_total.append(seq.seq.count('c')+seq.seq.count('C'))
    
    return int(sum(c_total))

# countGs returns the count of G/g in the alignment
def countGs():
    global pruned_alignment
    g_total = []
    for seq in pruned_alignment:
        g_total.append(seq.seq.count('g')+seq.seq.count('G'))

    return int(sum(g_total))

# countTs returns the count of T/t in the alignment        
def countTs():
    global pruned_alignment
    t_total = []
    for seq in pruned_alignment:
        t_total.append(seq.seq.count('t')+seq.seq.count('T'))

    return int(sum(t_total))

# countDegen returns the count of degenerate bases in the alignment
def countDegen():
    global pruned_alignment
    degen_total = []
    for seq in pruned_alignment:
        degen_total.append(seq.seq.count('w')+seq.seq.count('W')+seq.seq.count('s')+seq.seq.count('S')+seq.seq.count('m')+seq.seq.count('M')+seq.seq.count('k')+seq.seq.count('K')+seq.seq.count('r')+seq.seq.count('R')+seq.seq.count('y')+seq.seq.count('Y')+seq.seq.count('b')+seq.seq.count('B')+seq.seq.count('d')+seq.seq.count('D')+seq.seq.count('h')+seq.seq.count('H')+seq.seq.count('v')+seq.seq.count('V'))

    return int(sum(degen_total))

# countNs returns the count of N/n in the alignment
def countNs():
    global pruned_alignment
    n_total = []
    for seq in pruned_alignment:
        n_total.append(seq.seq.count('n')+seq.seq.count('N'))

    return int(sum(n_total))

# countGaps returns the count of "-" in the alignment
def countGaps():
    global pruned_alignment
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
    raw_spp = []
    for seq_record in raw_alignment:
        raw_spp.append(str(seq_record.id))

    # If fewer than three species exist in alignment or species list, raise exception
    if (len(spp_list) < 3) or (len(raw_spp) < 3):
        sys.exit("ERROR: Cannot process fewer than 3 species...")
        
    # If requested species are not in alignment, raise exception
    if not all(elem in spp_list for elem in raw_spp):
        sys.exit("ERROR: Requested species not found in "+os.path.basename(alignment_path)+"...")
    
    # Re-order alignment to sorted order of species list
    else:
        # Create dummy alignment        
        global pruned_alignment
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
    a_count = countAs()
    c_count = countCs()
    g_count = countGs()
    t_count = countTs()
    degen_count = countDegen()
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
    percent_degen = float(degen_count)/(alignment_length*len(spp_list))
    
    return_df = pd.DataFrame([[alignment_name,alignment_length,percent_gc,percent_degen,percent_n,percent_gap]],columns=['Alignment_Name','Alignment_Length','Percent_GC','Percent_Degenerate','Percent_N','Percent_Gap'])
    
    temp_df = tempfile.NamedTemporaryFile()
    temp_df_name = temp_df.name
    temp_df.close()
    with open(temp_df_name, "wb") as f:
        pickle.dump(return_df, f)
    f.close()
    
    print(temp_df_name)
    
def main(path_to_align,spp_info,align_name):
    fetchAlignmentComposition(path_to_align,spp_info,align_name)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])