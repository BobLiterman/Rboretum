#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and return the sequence names
"""
from Bio import AlignIO, SeqIO

def getAlignmentSpecies(alignment_path):
    
    formats = {'nex': 'nexus', 'nexus': 'nexus',
               'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
               'fa': 'fasta', 'fasta': 'fasta'}

    fformat = formats[alignment_path.split('.')[-1]]
    
    try:
        align = AlignIO.read(alignment_path, fformat)

        # Get alignment species
        align_species = list()
        for seq_record in align:
            align_species.append(seq_record.id)

        sort_species = [str(i) for i in align_species]
        sort_species.sort()
        return sort_species
    
    except:
        return []