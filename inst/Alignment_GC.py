#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and return the percent GC.
"""

import sys
import os
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqUtils import GC


def getGC(alignment_path):

    formats = {'nex': 'nexus', 'nexus': 'nexus',
               'phy': 'phylip-relaxed', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip-relaxed',
               'fa': 'fasta', 'fasta': 'fasta'}
    fformat = formats[alignment_path.split('.')[-1]]

    pruned_alignment = AlignIO.read(alignment_path, fformat)

    all_total = []
    gc_total = []
    for seq in pruned_alignment:
        all_total.append(seq.seq.count('a')+seq.seq.count('t')+seq.seq.count('A')+seq.seq.count('T')+seq.seq.count('g')+seq.seq.count('c')+seq.seq.count('G')+seq.seq.count('C'))
        gc_total.append(seq.seq.count('g')+seq.seq.count('c')+seq.seq.count('G')+seq.seq.count('C'))

    return(sum(gc_total)/sum(all_total))
