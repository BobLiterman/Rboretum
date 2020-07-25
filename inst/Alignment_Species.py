#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and return the sequence names as a semicolon-separated string
"""
from Bio import AlignIO, SeqIO

def fetchAlignmentSpecies(alignment_path):
    
    formats = {'nex': 'nexus', 'nexus': 'nexus',
               'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
               'fa': 'fasta', 'fasta': 'fasta'}

    fformat = formats[alignment_path.split('.')[-1]]
    
    try:
        align = AlignIO.read(alignment_path, fformat)

        # Get alignment species
        align_species = list()
        for seq_record in align:
            align_species.append(str(seq_record.id))
        
        return align_species
    
    except:
        sys.exit("ERROR: Cannot process file at "+alignment_path)