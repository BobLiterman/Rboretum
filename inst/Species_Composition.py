#!/usr/bin/env python3
"""
Arguments:
(1) path_to_align: Absolute path to alignment file
(2) spp_info: Semicolon-separated list of taxa to subset from alignment
(3) align_name: Label for dataset/subset
"""
import sys
import os

if "." not in sys.path:
    sys.path.append(".")
    
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from .readAlignment import getPrunedAlignment
    
def fetchSpeciesComposition(path_to_align,spp_info,align_name):
        
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

def main(path_to_align,spp_info,align_name):
    return(fetchSpeciesComposition(path_to_align,spp_info,align_name))
    
if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])