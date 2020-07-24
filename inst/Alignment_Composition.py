#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and:

    1) Prune down to a supplied subset of taxa
    2) Breakdown the signal at each site
"""

import sys
import os
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from operator import itemgetter
from collections import Counter
import multiprocessing as mp
from itertools import chain
from getAlignmentSpecies import getAlignmentSpecies

def getAlignmentComposition(alignment_path1,alignment_name1,spp_info):
    #ALIGNMENT_NAME########
    global alignment_path
    global alignment_name
    global info_gap
    global spp_list
    global bases

    # Prepare arguments
    if not alignment_path1:
        sys.exit("ERROR: No alignment path specified.")
    else:
        alignment_path = str(alignment_path1)
    
    # Get species list, or process all species
    if not spp_info:
        spp_list = sorted(getAlignmentSpecies(alignment_path))
    else:
        if len(spp_info) == 1 and ";" in spp_info:
            spp_list = sorted(str(spp_string).split(";"))
        else:
            spp_list = sorted(spp_info)

    # If alignment_filename contains all species from spp_list (>= 3 species), continue
    if not getPrunedAlignment():
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path)+" with provided species list.")
    
    # Set alignment name
    if not alignment_name1:
        alignment_name = str(os.path.basename(alignment_path))
    else:
        try:
            alignment_name = str(alignment_name1)
        except:
            alignment_name = str(os.path.basename(alignment_path))
    
    alignment_length = int(pruned_alignment.get_alignment_length())
    
    a_count = countAs()
    c_count = countCs()
    g_count = countGs()
    t_count = countTs()

    all_base_count = sum(a_count,c_count,g_count,t_count)
    
    gc_count = sum(c_count,g_count)
    percent_gc = gc_count/all_base_count

    n_count = countNs()
    percent_n = n_count/alignment_length

    gap_count = countGaps()
    percent_gap = gap_count/alignment_length
    
    return(pd.DataFrame([[alignment_name,alignment_length,a_count,c_count,g_count,t_count,percent_gc,percent_n,percent_gap]],columns=['Alignment','Alignment_Length','A_Count','C_Count','G_Count','T_Count','Percent_GC','Percent_N','Percent_Gap']))

def getPrunedAlignment():
    # getPrunedAlignment returns True if:
        # (1) All species in tree are in alignment
        # (2) Three or more species are present.
    # If True,sets global pruned_alignment is pruned and alphabetized match tree species

    formats = {'nex': 'nexus', 'nexus': 'nexus',
               'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
               'fa': 'fasta', 'fasta': 'fasta'}
    
    fformat = formats[alignment_path.split('.')[-1]]
    
    # Read in alignment
    try:
        raw_alignment = AlignIO.read(alignment_path, fformat)
    except:
        return False

    # Get alignment species
    raw_spp = list()
    for seq_record in raw_alignment:
        raw_spp.append(str(seq_record.id))
    
    # Get sorted species list
    global sort_species
    sort_species = sorted(raw_spp)
    
    # If requested species are not in alignment, or species list/alignment contains fewer than 3 species, return False
    if((len(list(set(spp_list)-set(raw_spp))) > 0) or (len(list(set(spp_list))) < 3) or (len(list(set(raw_spp))) < 3)):
        return False
    
    # Re-order alignment to sorted order of species list
    else:
        # Create dummy alignment        
        global pruned_alignment
        pruned_alignment = raw_alignment[0:0]
        
        # Populate alignment by adding taxa sorted by taxon ID
        for i in sort_species:
            pruned_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))
        
        # Return True to indicate a valid alignment was processed
        return True

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

    return sum(n_total)

def countGaps():
    # countGaps returns the count of "-" in the alignment
    global pruned_alignment
    gap_total = []
    for seq in pruned_alignment:
        gap_total.append(seq.seq.count('-'))

    return sum(gap_total)