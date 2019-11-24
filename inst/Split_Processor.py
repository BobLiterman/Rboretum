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

# Returns indexes  of occurrences of list items in a list or string
# https://stackoverflow.com/questions/13009675/find-all-the-occurrences-of-a-character-in-a-string
def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter in ch]

# Returns True if:
# (1) All species in tree are in alignment
# (2) Three or more species are present.
# If True,sets global pruned_alignment is pruned and alphabetized match tree species
def subset_alignment():

    global pruned_alignment

    formats = {'nex': 'nexus', 'nexus': 'nexus',
               'phy': 'phylip', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip',
               'fa': 'fasta', 'fasta': 'fasta'}
    fformat = formats[alignment_path.split('.')[-1]]

    pruned_alignment = AlignIO.read(alignment_path, fformat)

    # Get alignment species
    raw_spp = list()
    for seq_record in pruned_alignment:
        raw_spp.append(seq_record.id)

    global align_list
    align_list = sorted(raw_spp)

    # (1) If requested speces are not in alignment, or species list/alignment contains fewer than 3 species, return False

    if((len(list(set(spp_list)-set(raw_spp))) > 0) or (len(list(set(spp_list))) < 3) or (len(list(set(raw_spp))) < 3)):
        return False

    else:
        # Re-order alignment to sorted order of species list
        # Create dummy alignment
        empty_align = pruned_alignment[0:0]

        for i in sorted(spp_list):
            empty_align.add_sequence(str(pruned_alignment[raw_spp.index(i)].id), str(
                pruned_alignment[raw_spp.index(i)].seq))

        pruned_alignment = empty_align
        return True

# Takes integer argument and describes pattern of variation at that position in pruned_alignment
def create_site_dict(pos):

    global pruned_alignment
    temp_dict = dict()

    # Pull alignment column
    seq_string = pruned_alignment[:, pos]

    # Establish valid bases
    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c']

    # Get column sequence as string
    site_set = Counter(seq_string)

    # Get different alleles as list
    allele_list = list(site_set)

    # If column contains non-base value (ACTG,actg), return "a"
    if(not set(allele_list).issubset(set(bases))):
        temp_dict[pos] = "non_base"

    # If column contains no 'non-bases'...
    else:
        allele_count = len(allele_list)
        site_counts = list(site_set.values())

        # If one allele detected, the return is invariant ("b")
        if(allele_count == 1):
            temp_dict[pos] = "invariant"

        # If more than one allele is detected...
        else:

            # If any variation is singleton in nature, return "c"
            if(1 in site_counts):
                temp_dict[pos] = "singleton"

            # If there are no singletons, and no non-bases, sites are either biallelic ("d"), triallelic ("e"), or quadallelic ("f")
            else:
                if(allele_count == 2):
                    temp_dict[pos] = "biallelic"
                if(allele_count == 3):
                    temp_dict[pos] = "triallelic"
                if(allele_count == 4):
                    temp_dict[pos] = "quadallelic"

    return temp_dict

# IF GAPS ARE ENABLED:
def create_site_dict_gap(pos):

    global pruned_alignment
    temp_dict = dict()

    # Pull alignment column
    seq_string = pruned_alignment[:, pos]

    # Establish valid bases
    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c', '-']

    # Get column sequence as string
    site_set = Counter(seq_string)

    # Get different alleles as list
    allele_list = list(site_set)

    # If column contains non-base value (ACTG-,actg-), return "a"
    if(not set(allele_list).issubset(set(bases))):
        temp_dict[pos] = "non_base"

    # If column contains no 'non-bases'...
    else:
        allele_count = len(allele_list)
        site_counts = list(site_set.values())

        # If column includes gaps...
        if("-" in allele_list):

            # If one allele detected, the return is invariant ("g")
            if(allele_count == 1):
                temp_dict[pos] = "gap_invariant"

            # If more than one allele is detected...
            else:

                # If any variation is singleton in nature, return "h"
                if(1 in site_counts):
                    temp_dict[pos] = "gap_singleton"

                # If there are no singletons, and no non-bases, sites are either biallelic ("i"), triallelic ("j"), quadallelic ("k"), or pentallelic ("l")
                else:
                    if(allele_count == 2):
                        temp_dict[pos] = "gap_biallelic"
                    if(allele_count == 3):
                        temp_dict[pos] = "gap_triallelic"
                    if(allele_count == 4):
                        temp_dict[pos] = "gap_quadallelic"
                    if(allele_count == 5):
                        temp_dict[pos] = "pentallelic"

        # If column lacks gaps, key values are same as non-gap
        else:
            # If one allele detected, the return is invariant ("b")
            if(allele_count == 1):
                temp_dict[pos] = "invariant"

            # If more than one allele is detected...
            else:

                # If any variation is singleton in nature, return "c"
                if(1 in site_counts):
                    temp_dict[pos] = "singleton"

                # If there are no singletons, and no non-bases, sites are either biallelic ("d"), triallelic ("e"), or quadallelic ("f")
                else:
                    if(allele_count == 2):
                        temp_dict[pos] = "biallelic"
                    if(allele_count == 3):
                        temp_dict[pos] = "triallelic"
                    if(allele_count == 4):
                        temp_dict[pos] = "quadallelic"

    return temp_dict

def process_nonbase(pos):

    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c']

    seq_string = pruned_alignment[:, pos]

    non_base_list = list(set(seq_string) - set(bases))
    non_base_index = findOccurrences(seq_string, non_base_list)
    non_base_taxa = itemgetter(*non_base_index)(sorted(spp_list))
    non_base_count = int(len(non_base_index))

    if len(non_base_index) == 1:
        non_base_taxa_string = non_base_taxa
    if len(non_base_index) > 1:
        non_base_taxa_string = ";".join(sorted(non_base_taxa))

    if len(spp_list) - len(non_base_taxa) >= 3:
        
        # Get species list and sequence data from taxa with appropriate bases
        new_seq = [v for k,v in enumerate(seq_string) if k not in non_base_index]
        new_spp = [v for k,v in enumerate(spp_list) if k not in non_base_index]

        # Get different alleles as list
        site_set = Counter(new_seq)
        allele_list = list(site_set)
        allele_count = len(allele_list)
        site_counts = list(site_set.values())
    
        if len(site_counts)==1:

            return(pd.DataFrame([[pos,'invariant',non_base_taxa_string,non_base_count,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif 1 in site_counts:
            
            singleton_list = [base for base, count in site_set.items() if count == 1]
            singelton_index = findOccurrences(new_seq, singleton_list)
            singleton_taxa = itemgetter(*singelton_index)(sorted(new_spp))
            
            if len(singelton_index) == 1:
                taxa =  singleton_taxa
            elif len(singelton_index) > 1:
                taxa = ";".join(sorted(singleton_taxa))

            if len(new_spp) - len(singleton_taxa) >= 3:
                    
                # Get species list and sequence data from taxa with appropriate bases
                new_seq2 = [v for k,v in enumerate(new_seq) if k not in singelton_index]
                new_spp2 = [v for k,v in enumerate(new_spp) if k not in singelton_index]
        
                # Get different alleles as list
                site_set2 = Counter(new_seq2)
                allele_list2 = list(site_set2)
                allele_count2 = len(allele_list2)
                site_counts2 = list(site_set2.values())
            
                if len(site_counts2)==1:
                    return(pd.DataFrame([[pos,'invariant',non_base_taxa_string,non_base_count,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
                elif len(site_counts2)==2:
                    
                    split_1_base = sorted(list(set(new_seq2)))[0]
                    base_1_index = findOccurrences(new_seq2, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp2)))
                    
                    split_2_base = sorted(list(set(new_seq2)))[1]
                    base_2_index = findOccurrences(new_seq2, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp2)))
        
                    return(pd.DataFrame([[pos,'biallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
                elif len(site_counts2)==3:
        
                    split_1_base = sorted(list(set(new_seq2)))[0]
                    base_1_index = findOccurrences(new_seq2, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp2)))
                    
                    split_2_base = sorted(list(set(new_seq2)))[1]
                    base_2_index = findOccurrences(new_seq2, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp2)))
        
                    split_3_base = sorted(list(set(new_seq2)))[2]
                    base_3_index = findOccurrences(new_seq2, split_3_base)
                    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp2)))          
        
                    return(pd.DataFrame([[pos,'triallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
                elif len(site_counts2)==4:
        
                    split_1_base = sorted(list(set(new_seq2)))[0]
                    base_1_index = findOccurrences(new_seq2, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp2)))
                    
                    split_2_base = sorted(list(set(new_seq2)))[1]
                    base_2_index = findOccurrences(new_seq2, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp2)))
        
                    split_3_base = sorted(list(set(new_seq2)))[2]
                    base_3_index = findOccurrences(new_seq2, split_3_base)
                    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp2)))
                    
                    split_4_base = sorted(list(set(new_seq2)))[3]
                    base_4_index = findOccurrences(new_seq2, split_4_base)
                    base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp2)))  

                    return(pd.DataFrame([[pos,'quadallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
            else:
                return(pd.DataFrame([[pos,'singleton',non_base_taxa_string,non_base_count,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            
        elif len(site_counts)==2:
            
            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))
            
            return(pd.DataFrame([[pos,'biallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
        elif len(site_counts)==3:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))          
            
            return(pd.DataFrame([[pos,'triallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
        elif len(site_counts)==4:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))
            
            split_4_base = sorted(list(set(new_seq)))[3]
            base_4_index = findOccurrences(new_seq, split_4_base)
            base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp)))  
            
            return(pd.DataFrame([[pos,'quadallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    else:
        return(pd.DataFrame([[pos,'non_base',non_base_taxa_string,non_base_count,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def process_nonbase_gap(pos):

    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c','-']

    seq_string = pruned_alignment[:, pos]

    non_base_list = list(set(seq_string) - set(bases))
    non_base_index = findOccurrences(seq_string, non_base_list)
    non_base_taxa = itemgetter(*non_base_index)(sorted(spp_list))
    non_base_count = int(len(non_base_index))

    if len(non_base_index) == 1:
        non_base_taxa_string = non_base_taxa
    if len(non_base_index) > 1:
        non_base_taxa_string = ";".join(sorted(non_base_taxa))

    if len(spp_list) - len(non_base_taxa) >= 3:
        
        # Get species list and sequence data from taxa with appropriate bases
        new_seq = [v for k,v in enumerate(seq_string) if k not in non_base_index]
        new_spp = [v for k,v in enumerate(spp_list) if k not in non_base_index]

        # Get different alleles as list
        site_set = Counter(new_seq)
        allele_list = list(site_set)
        allele_count = len(allele_list)
        site_counts = list(site_set.values())
    
        if len(site_counts)==1:
            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_invariant',non_base_taxa_string,non_base_count,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'invariant',non_base_taxa_string,non_base_count,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif 1 in site_counts:
            
            singleton_list = [base for base, count in site_set.items() if count == 1]
            singelton_index = findOccurrences(new_seq, singleton_list)
            singleton_taxa = itemgetter(*singelton_index)(sorted(new_spp))
            
            if len(singelton_index) == 1:
                taxa =  singleton_taxa
            elif len(singelton_index) > 1:
                taxa = ";".join(sorted(singleton_taxa))

            if len(new_spp) - len(singleton_taxa) >= 3:
                    
                # Get species list and sequence data from taxa with appropriate bases
                new_seq2 = [v for k,v in enumerate(new_seq) if k not in singelton_index]
                new_spp2 = [v for k,v in enumerate(new_spp) if k not in singelton_index]
        
                # Get different alleles as list
                site_set2 = Counter(new_seq2)
                allele_list2 = list(site_set2)
                allele_count2 = len(allele_list2)
                site_counts2 = list(site_set2.values())
            
                if len(site_counts2)==1:
                    if '-' in new_seq:
                        return(pd.DataFrame([[pos,'gap_invariant',non_base_taxa_string,non_base_count,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
                    else:
                        return(pd.DataFrame([[pos,'invariant',non_base_taxa_string,non_base_count,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
                elif len(site_counts2)==2:
                    
                    split_1_base = sorted(list(set(new_seq2)))[0]
                    base_1_index = findOccurrences(new_seq2, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp2)))
                    
                    split_2_base = sorted(list(set(new_seq2)))[1]
                    base_2_index = findOccurrences(new_seq2, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp2)))
        
                    if '-' in new_seq:
                        return(pd.DataFrame([[pos,'gap_biallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
                    else:
                        return(pd.DataFrame([[pos,'biallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
                elif len(site_counts2)==3:
        
                    split_1_base = sorted(list(set(new_seq2)))[0]
                    base_1_index = findOccurrences(new_seq2, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp2)))
                    
                    split_2_base = sorted(list(set(new_seq2)))[1]
                    base_2_index = findOccurrences(new_seq2, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp2)))
        
                    split_3_base = sorted(list(set(new_seq2)))[2]
                    base_3_index = findOccurrences(new_seq2, split_3_base)
                    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp2)))          
        
                    if '-' in new_seq:
                        return(pd.DataFrame([[pos,'gap_triallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
                    else:
                        return(pd.DataFrame([[pos,'triallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
                elif len(site_counts2)==4:
        
                    split_1_base = sorted(list(set(new_seq2)))[0]
                    base_1_index = findOccurrences(new_seq2, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp2)))
                    
                    split_2_base = sorted(list(set(new_seq2)))[1]
                    base_2_index = findOccurrences(new_seq2, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp2)))
        
                    split_3_base = sorted(list(set(new_seq2)))[2]
                    base_3_index = findOccurrences(new_seq2, split_3_base)
                    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp2)))
                    
                    split_4_base = sorted(list(set(new_seq2)))[3]
                    base_4_index = findOccurrences(new_seq2, split_4_base)
                    base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp2)))  
        
                    if '-' in new_seq:
                        return(pd.DataFrame([[pos,'gap_quadallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
                    else:
                        return(pd.DataFrame([[pos,'quadallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
                elif len(site_counts2)==5:
        
                    split_1_base = sorted(list(set(new_seq2)))[0]
                    base_1_index = findOccurrences(new_seq2, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp2)))
                    
                    split_2_base = sorted(list(set(new_seq2)))[1]
                    base_2_index = findOccurrences(new_seq2, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp2)))
        
                    split_3_base = sorted(list(set(new_seq2)))[2]
                    base_3_index = findOccurrences(new_seq2, split_3_base)
                    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp2)))
                    
                    split_4_base = sorted(list(set(new_seq2)))[3]
                    base_4_index = findOccurrences(new_seq2, split_4_base)
                    base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp2)))  
        
                    split_5_base = sorted(list(set(new_seq2)))[4]
                    base_5_index = findOccurrences(new_seq2, split_5_base)
                    base_5_taxa = ";".join(itemgetter(*base_5_index)(sorted(new_spp2)))  
        
                    return(pd.DataFrame([[pos,'pentallelic',non_base_taxa_string,non_base_count,taxa,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,base_5_taxa]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
            else:
                if '-' in new_seq:
                    return(pd.DataFrame([[pos,'gap_singleton',non_base_taxa_string,non_base_count,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
                else:
                    return(pd.DataFrame([[pos,'singleton',non_base_taxa_string,non_base_count,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==2:
            
            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_biallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'biallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==3:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))          

            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_triallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'triallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==4:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))
            
            split_4_base = sorted(list(set(new_seq)))[3]
            base_4_index = findOccurrences(new_seq, split_4_base)
            base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp)))  

            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_quadallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'quadallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==5:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))
            
            split_4_base = sorted(list(set(new_seq)))[3]
            base_4_index = findOccurrences(new_seq, split_4_base)
            base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp)))  

            split_5_base = sorted(list(set(new_seq)))[4]
            base_5_index = findOccurrences(new_seq, split_5_base)
            base_5_taxa = ";".join(itemgetter(*base_5_index)(sorted(new_spp)))  

            return(pd.DataFrame([[pos,'pentallelic',non_base_taxa_string,non_base_count,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,base_5_taxa]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    else:
        return(pd.DataFrame([[pos,'non_base',non_base_taxa_string,non_base_count,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def process_invariant(pos):
    seq_string = pruned_alignment[:, pos]
    
    if '-' in seq_string:
        return(pd.DataFrame([[pos,'gap_invariant',np.nan,0,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
    else:
        return(pd.DataFrame([[pos,'invariant',np.nan,0,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def process_singletons(pos):

    seq_string = pruned_alignment[:, pos]

    # Get column sequence as string
    site_set = dict(Counter(seq_string))

    # https://stackoverflow.com/questions/8023306/get-key-by-value-in-dictionary
    singleton_list = [base for base, count in site_set.items() if count == 1]
    singelton_index = findOccurrences(seq_string, singleton_list)
    singleton_taxa = itemgetter(*singelton_index)(sorted(spp_list))

    if len(singelton_index) == 1:
        taxa =  singleton_taxa
    if len(singelton_index) > 1:
        taxa = ";".join(sorted(singleton_taxa))

    if len(spp_list) - len(singleton_taxa) >= 3:
        
        # Get species list and sequence data from taxa with appropriate bases
        new_seq = [v for k,v in enumerate(seq_string) if k not in singelton_index]
        new_spp = [v for k,v in enumerate(spp_list) if k not in singelton_index]

        # Get different alleles as list
        site_set = Counter(new_seq)
        allele_list = list(site_set)
        allele_count = len(allele_list)
        site_counts = list(site_set.values())
    
        if len(site_counts)==1:
            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_invariant',np.nan,0,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'invariant',np.nan,0,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==2:
            
            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_biallelic',np.nan,0,taxa,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'biallelic',np.nan,0,taxa,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==3:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))          

            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_triallelic',np.nan,0,taxa,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'triallelic',np.nan,0,taxa,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==4:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))
            
            split_4_base = sorted(list(set(new_seq)))[3]
            base_4_index = findOccurrences(new_seq, split_4_base)
            base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp)))  

            if '-' in new_seq:
                return(pd.DataFrame([[pos,'gap_quadallelic',np.nan,0,taxa,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
            else:
                return(pd.DataFrame([[pos,'quadallelic',np.nan,0,taxa,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

        elif len(site_counts)==5:

            split_1_base = sorted(list(set(new_seq)))[0]
            base_1_index = findOccurrences(new_seq, split_1_base)
            base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
            
            split_2_base = sorted(list(set(new_seq)))[1]
            base_2_index = findOccurrences(new_seq, split_2_base)
            base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))

            split_3_base = sorted(list(set(new_seq)))[2]
            base_3_index = findOccurrences(new_seq, split_3_base)
            base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(new_spp)))
            
            split_4_base = sorted(list(set(new_seq)))[3]
            base_4_index = findOccurrences(new_seq, split_4_base)
            base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(new_spp)))  

            split_5_base = sorted(list(set(new_seq)))[4]
            base_5_index = findOccurrences(new_seq, split_5_base)
            base_5_taxa = ";".join(itemgetter(*base_5_index)(sorted(new_spp)))  

            return(pd.DataFrame([[pos,'pentallelic',np.nan,0,taxa,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,base_5_taxa]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    else:
        return(pd.DataFrame([[pos,'singleton',np.nan,0,taxa,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))


def process_biallelic(pos):

    seq_string = pruned_alignment[:, pos]

    # Get column sequence as string
    split_1_base = sorted(list(set(seq_string)))[0]
    split_2_base = sorted(list(set(seq_string)))[1]

    base_1_index = findOccurrences(seq_string, split_1_base)
    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(spp_list)))

    base_2_index = findOccurrences(seq_string, split_2_base)
    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(spp_list)))
    
    if '-' in seq_string:
        return(pd.DataFrame([[pos,'gap_biallelic',np.nan,0,np.nan,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
    else:
        return(pd.DataFrame([[pos,'biallelic',np.nan,0,np.nan,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
      
def process_triallelic(pos):

    seq_string = pruned_alignment[:, pos]

    # Get column sequence as string
    split_1_base = sorted(list(set(seq_string)))[0]
    split_2_base = sorted(list(set(seq_string)))[1]
    split_3_base = sorted(list(set(seq_string)))[2]

    base_1_index = findOccurrences(seq_string, split_1_base)
    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(spp_list)))

    base_2_index = findOccurrences(seq_string, split_2_base)
    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(spp_list)))

    base_3_index = findOccurrences(seq_string, split_3_base)
    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(spp_list)))

    if '-' in seq_string:
        return(pd.DataFrame([[pos,'gap_triallelic',np.nan,0,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
    else:
        return(pd.DataFrame([[pos,'triallelic',np.nan,0,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def process_quadallelic(pos):

    seq_string = pruned_alignment[:, pos]

    # Get column sequence as string
    split_1_base = sorted(list(set(seq_string)))[0]
    split_2_base = sorted(list(set(seq_string)))[1]
    split_3_base = sorted(list(set(seq_string)))[2]
    split_4_base = sorted(list(set(seq_string)))[3]

    base_1_index = findOccurrences(seq_string, split_1_base)
    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(spp_list)))

    base_2_index = findOccurrences(seq_string, split_2_base)
    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(spp_list)))

    base_3_index = findOccurrences(seq_string, split_3_base)
    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(spp_list)))

    base_4_index = findOccurrences(seq_string, split_4_base)
    base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(spp_list)))

    if '-' in seq_string:
        return(pd.DataFrame([[pos,'gap_quadallelic',np.nan,0,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
    else:
        return(pd.DataFrame([[pos,'quadallelic',np.nan,0,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def process_pentallelic(pos):

    seq_string = pruned_alignment[:, pos]

    # Get column sequence as string
    split_1_base = sorted(list(set(seq_string)))[0]
    split_2_base = sorted(list(set(seq_string)))[1]
    split_3_base = sorted(list(set(seq_string)))[2]
    split_4_base = sorted(list(set(seq_string)))[3]
    split_5_base = sorted(list(set(seq_string)))[4]

    base_1_index = findOccurrences(seq_string, split_1_base)
    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(spp_list)))

    base_2_index = findOccurrences(seq_string, split_2_base)
    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(spp_list)))

    base_3_index = findOccurrences(seq_string, split_3_base)
    base_3_taxa = ";".join(itemgetter(*base_3_index)(sorted(spp_list)))

    base_4_index = findOccurrences(seq_string, split_4_base)
    base_4_taxa = ";".join(itemgetter(*base_4_index)(sorted(spp_list)))

    base_5_index = findOccurrences(seq_string, split_5_base)
    base_5_taxa = ";".join(itemgetter(*base_5_index)(sorted(spp_list)))

    return(pd.DataFrame([[pos,'pentallelic',np.nan,0,np.nan,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,base_5_taxa]], columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def sitePasser(pos_pattern):
    pos = pos_pattern[0]
    pattern = pos_pattern[1]
    
    if pattern in ['invariant','gap_invariant']:
        return process_invariant(pos)

    elif pattern == 'non_base':
        if info_gap == "0":
            return process_nonbase(pos)
        elif info_gap == "1":
            return process_nonbase_gap(pos)
    
    elif pattern in ['singleton','gap_singleton']:
        return process_singletons(pos)
        
    elif pattern in ['biallelic','gap_biallelic']:
        return process_biallelic(pos)

    elif pattern in ['triallelic','gap_triallelic']:
        return process_triallelic(pos)

    elif pattern in ['quadallelic','gap_quadallelic']:
        return process_quadallelic(pos)

    elif pattern == 'pentallelic':
        return process_pentallelic(pos)
      

def getSplitSupport(alignment_path1,info_gap1,spp_string):

    global alignment_path
    global info_gap
    global spp_list

    alignment_path = str(alignment_path1)
    info_gap = str(info_gap1)
    spp_list = sorted(str(spp_string).split(";"))
    
    # Get site split info, or False if alignment cannot be processed for given species set
    all_sites_df = splitMain()

    if all_sites_df.shape[0] == 0:
        return -1
    else:
        return all_sites_df

def splitMain():

    # Create dummy return for errors
    empty = pd.DataFrame()

    # If alignment_filename contains all species from spp_list (>= 3 species), continue
    if subset_alignment():

        # Define each position in the alignment by its 0-based position
        alignment_positions = range(0, pruned_alignment.get_alignment_length())
        # Detect CPUs and create a pool of workers
        if mp.cpu_count() == 1:
            pool_cpu = 1
        else:
            pool_cpu = mp.cpu_count() - 1

        # For each position, grab bases and describe variation pattern
        global site_dict

        if info_gap == "0":
            with mp.Pool(pool_cpu) as pool:
                site_dict = pool.map(create_site_dict, alignment_positions)
                print("Processed site patterns with gaps interpreted as missing data...assessing signal...")

        elif info_gap == "1":
            with mp.Pool(pool_cpu) as pool:
                site_dict = pool.map(create_site_dict_gap, alignment_positions)
                print("Processed site patterns with gaps interpreted as informative indels...assessing signal...")

        # Process site patterns
        pattern_list = []

        for i in range(0, len(site_dict)):
            pattern_list.append([i, site_dict[i][i]])

        with mp.Pool(pool_cpu) as pool:
            results = pool.map(sitePasser,pattern_list,1000)
        
        print("Signal processed...compiling final dataframe...")
        results = pd.concat(results).sort_values(by=['Zeroed_Site_Position']).reset_index()

        return results
    else:
        return empty
