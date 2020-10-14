#!/usr/bin/env python3
from Bio import AlignIO, SeqIO

# countAs returns the count of A/a in the alignment
def countAs(pruned_alignment):
    a_total = []
    for seq in pruned_alignment:
        a_total.append(seq.seq.count('a')+seq.seq.count('A'))

    return int(sum(a_total))

# countCs returns the count of C/c in the alignment
def countCs(pruned_alignment):
    c_total = []
    for seq in pruned_alignment:
        c_total.append(seq.seq.count('c')+seq.seq.count('C'))
    
    return int(sum(c_total))

# countGs returns the count of G/g in the alignment
def countGs(pruned_alignment):
    g_total = []
    for seq in pruned_alignment:
        g_total.append(seq.seq.count('g')+seq.seq.count('G'))

    return int(sum(g_total))

# countTs returns the count of T/t in the alignment        
def countTs(pruned_alignment):
    t_total = []
    for seq in pruned_alignment:
        t_total.append(seq.seq.count('t')+seq.seq.count('T'))

    return int(sum(t_total))

# countNs returns the count of N/n in the alignment
def countNs(pruned_alignment):
    n_total = []
    for seq in pruned_alignment:
        n_total.append(seq.seq.count('n')+seq.seq.count('N'))

    return int(sum(n_total))

# countGaps returns the count of "-" in the alignment
def countGaps(pruned_alignment):
    gap_total = []
    for seq in pruned_alignment:
        gap_total.append(seq.seq.count('-'))

    return int(sum(gap_total))