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

def alignmentPercents(alignment_path1,spp_string):

    global alignment_path
    global info_gap
    global spp_list
    global bases

    # Prepare arguments
    alignment_path = str(alignment_path1)
    spp_list = sorted(str(spp_string).split(";"))

    # If alignment_filename contains all species from spp_list (>= 3 species), continue
    if not getPrunedAlignment():
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path)+" with provided species list.")
    
    return(pd.DataFrame([[os.path.basename(alignment_path),int(pruned_alignment.get_alignment_length()),getPercentGC(),getPercentN(),getPercentGap()]], columns=['Alignment','Alignment_Length','Percent_GC','Percent_N','Percent_Gap']))

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

def getPercentGC():
    # getPercentGC returns the percent GC content for the entire alignment
    # Bases that are not A/a, T/t, C/c, or G/g are not considered as part of the total base count (i.e. Ns and gaps are ignored)
    global pruned_alignment
    all_total = []
    gc_total = []
    for seq in pruned_alignment:
        all_total.append(seq.seq.count('a')+seq.seq.count('t')+seq.seq.count('A')+seq.seq.count('T')+seq.seq.count('g')+seq.seq.count('c')+seq.seq.count('G')+seq.seq.count('C'))
        gc_total.append(seq.seq.count('g')+seq.seq.count('c')+seq.seq.count('G')+seq.seq.count('C'))

    return(sum(gc_total)/sum(all_total))

def getPercentN():
    # getPercentN returns percentage of the entire alignment that has N base calls
    global pruned_alignment
    all_total = []
    n_total = []
    for seq in pruned_alignment:
        all_total.append(seq.seq.count('-')+seq.seq.count('n')+seq.seq.count('N')+seq.seq.count('a')+seq.seq.count('t')+seq.seq.count('A')+seq.seq.count('T')+seq.seq.count('g')+seq.seq.count('c')+seq.seq.count('G')+seq.seq.count('C'))
        n_total.append(seq.seq.count('n')+seq.seq.count('N'))

    return(sum(n_total)/sum(all_total))

def getPercentGap():
    # getPercentGap returns percentage of the entire alignment with indels ('-')
    global pruned_alignment
    all_total = []
    gap_total = []
    for seq in pruned_alignment:
        all_total.append(seq.seq.count('-')+seq.seq.count('n')+seq.seq.count('N')+seq.seq.count('a')+seq.seq.count('t')+seq.seq.count('A')+seq.seq.count('T')+seq.seq.count('g')+seq.seq.count('c')+seq.seq.count('G')+seq.seq.count('C'))
        gap_total.append(seq.seq.count('-'))

    return(sum(gap_total)/sum(all_total))