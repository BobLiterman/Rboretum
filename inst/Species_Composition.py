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
from operator import itemgetter
from collections import Counter
import multiprocessing as mp
from itertools import chain

def fetchSpeciesComposition(path_to_align,align_name):
    
    # Set path to alignment
    global alignment_path
    alignment_path = str(path_to_align)
    
    # Set alignment name
    global alignment_name
    alignment_name = str(align_name)

    # If alignment at path_to_align contains all species from spp_list (>= 3 species), continue
    if not getPrunedAlignment():
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path))
    
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
        percent_gc.append(gc_i/float(total_i))

    percent_n = []
    for n_i in n_totals:
        percent_n.append(float(n_i)/alignment_length)
    
    percent_gap = []
    for gap_i in gap_totals:
        percent_gap.append(float(gap_i)/alignment_length)
    
    df = pd.DataFrame(list(zip(raw_spp,total_base_counts,n_totals,gap_totals,percent_gc,percent_n,percent_gap)), 
               columns =['Taxon','Total_Bases','Total_N','Total_Gaps','Percent_GC','Percent_N','Percent_Gap'])
    df['Alignment_Name'] = alignment_name
    return df

def getPrunedAlignment():
    # getPrunedAlignment returns True if:
        # (1) All species in tree are in alignment
        # (2) Three or more species are present.
    # If True,sets global pruned_alignment is pruned and alphabetized match tree species

    # Read in alignment
    try:
        formats = {'nex': 'nexus', 'nexus': 'nexus',
                   'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
                   'fa': 'fasta', 'fasta': 'fasta'}
        
        fformat = formats[alignment_path.split('.')[-1]]
        raw_alignment = AlignIO.read(alignment_path, fformat)
        
        # Get alignment species
        global raw_spp
        raw_spp = list()
        
        for seq_record in raw_alignment:
            raw_spp.append(str(seq_record.id))
        sort_species = sorted(raw_spp)
                
        global pruned_alignment
        pruned_alignment = raw_alignment[0:0]
        
        # Populate alignment by adding taxa sorted by taxon ID
        for i in sort_species:
            pruned_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))
        
        # Return True to indicate a valid alignment was processed
        return True
    except:
        return False