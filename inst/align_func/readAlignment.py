#!/usr/bin/env python3
import sys
import os
import numpy as np
from Bio import AlignIO, SeqIO

def getPrunedAlignment():
    # getPrunedAlignment returns an alignment if:
        # (1) Three or more species are present in the alignment and spp_list [if specified]
        # (2) All species provided by spp_list are present in the alignment [if specified]
        # (3) Alignment can be found and read into Python
    
    # Read in alignment
    try:
        formats = {'nex': 'nexus', 'nexus': 'nexus',
                   'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
                   'fa': 'fasta', 'fasta': 'fasta'}
        
        fformat = formats[alignment_path.split('.')[-1]]
        raw_alignment = AlignIO.read(alignment_path, fformat)

    # If alignment cannot be read in, raise exception
    except:
        raise

    # Get species from raw alignment
    raw_spp = list()
    for seq_record in raw_alignment:
        raw_spp.append(str(seq_record.id))
    
    # Convert numeric IDs to string prior to sorting
    raw_spp = [str(i) for i in raw_spp] 
    raw_spp.sort()
    
    # If fewer than three species exist in alignment or species list, raise exception
    if (len(list(set(spp_list))) < 3) or (len(list(set(raw_spp))) < 3):
        raise
        
    # If requested species are not in alignment, raise exception
    spp_diff = list(set(spp_list) - set(raw_spp))

    if len(spp_diff) > 0:
        raise
    
    # Re-order alignment to sorted order of species list
    else:
        # Create dummy alignment        
        pruned_alignment = raw_alignment[0:0]
        
        # Populate alignment by adding taxa sorted by taxon ID
        for i in spp_list:
            pruned_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))
        
        # If resulting alignment is empty, raise exception
        if int(pruned_alignment.get_alignment_length()) == 0:
            raise
        else:
            # Return True to indicate a valid alignment was processed
            return pruned_alignment