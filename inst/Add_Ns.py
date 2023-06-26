#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and:

    1) Prune down to a supplied subset of taxa
    2) Breakdown the signal at each site
"""

import sys
import os
import numpy as np
import math
import pandas as pd
from Bio import AlignIO, SeqIO
from operator import itemgetter
from collections import Counter
import multiprocessing as mp
from itertools import chain
import copy
import pickle
import tempfile
from functools import partial
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Convert 'percent_n' of 'seq' to N
def add_n(seq,percent_n):
    
    # Get sequence + length
    seq = list(seq)
    seq_length = len(seq)
    
    # Get number of Ns to add
    num_n = math.floor(percent_n * seq_length)
    
    # Generate num_n positions
    pos = np.random.randint(1,seq_length,num_n)
    for i in pos:
        seq[i] = 'N'
    
    return ''.join(seq)
    
def simulate_missing(path_to_align,spp_info,prop_n):

    # Set path to alignment
    alignment_path = str(path_to_align)
    full_path = os.path.abspath(alignment_path)
    base_name = os.path.basename(full_path)
    dir_name = os.path.dirname(full_path)

    # Get species list
    spp_list = str(spp_info).split(",")
    
    # Get N proportions
    n_percents = [float(x) for x in str(prop_n).split(",")]

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

    if all(elem in spp_list  for elem in raw_spp):
        # Create dummy alignment  
        global pruned_alignment
        pruned_alignment = raw_alignment[0:0]
        
        # Populate alignment by adding taxa sorted by taxon ID
        for i in range(0,len(spp_list)):
            
            spp_to_add = spp_list[i]
            spp_n_percent = n_percents[i]
            
            raw_index = raw_spp.index(spp_to_add)
            raw_id = raw_alignment[raw_index].id
            raw_seq = raw_alignment[raw_index].seq
            
            new_seq = add_n(raw_seq,spp_n_percent)
            #pruned_alignment.add_sequence(str(raw_id), new_seq)
            pruned_alignment.append(SeqRecord(Seq(new_seq), id=str(raw_id)))
        
        # If resulting alignment is empty, raise exception
        if int(pruned_alignment.get_alignment_length()) == 0:
            sys.exit("ERROR: Alignment processed, but appears to have no bases...")
        else:
            with open(dir_name+'\\'+base_name.replace("."+alignment_path.split('.')[-1],"_SimN."+alignment_path.split('.')[-1]), "w") as handle:
                SeqIO.write(pruned_alignment, handle, "phylip")
                
    else:
        sys.exit("ERROR: Requested species not found in "+os.path.basename(alignment_path)+"...")

def main(path_to_align,spp_info,prop_n):
    simulate_missing(path_to_align,spp_info,prop_n)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])
