#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and:

    1) Prune down to a supplied subset of taxa
    2) Remove sites that have indels or missing data (Any column != [A,C,T,G,a,c,t,g])

    THEN

    3) Isolate singelton sites (alignment_sing.nex)
    4) Isolate all parsimony informative sites (alignment_pi.nex)
    5) Isolate biallelic PI sites (alignment_bi.nex)
    6) Isolate triallelic PI sites (alignmment_tri.nex)
    7) Isolate quadallelic PI sites (alignment_quad.nex)
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
def subset_alignment(alignment_path, spp_list):

    global pruned_alignment

    formats = {'nex': 'nexus', 'nexus': 'nexus',
               'phy': 'phylip-relaxed', 'phylip-relaxed': 'phylip-relaxed', 'phylip': 'phylip-relaxed',
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

    return(temp_dict)

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

    return(temp_dict)

def process_nonbase(pos):

    global pruned_alignment
    global spp_list

    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c']

    seq_string = pruned_alignment[:, pos]

    non_base_list = list(set(seq_string) - set(bases))
    non_base_index = findOccurrences(seq_string, non_base_list)
    non_base_taxa = itemgetter(*non_base_index)(sorted(spp_list))

    if len(non_base_index) == 1:
        return(non_base_taxa)
    if len(non_base_index) > 1:
        return(";".join(sorted(non_base_taxa)))

def process_nonbase_gap(pos):

    global pruned_alignment
    global spp_list

    bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c', '-']

    seq_string = pruned_alignment[:, pos]

    non_base_list = list(set(seq_string) - set(bases))
    non_base_index = findOccurrences(seq_string, non_base_list)
    non_base_taxa = itemgetter(*non_base_index)(sorted(spp_list))

    if len(non_base_index) == 1:
        return(non_base_taxa)
    if len(non_base_index) > 1:
        return(";".join(sorted(non_base_taxa)))

# Tries to rescue signal from columns with missing data for up to 'max_missing' taxa
def reprocess_nonbase(nonbase_df,max_missing):
    # Create empty dataframe
    converted_nonbase = nonbase_df[0:0]
    final_columns = ['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']

    # For each row with some missing data:
    for i in range(0,nonbase_df.shape[0]):
            
        temp_row = nonbase_df[i:i+1].reset_index().copy()
        temp_pos = int(temp_row['Zeroed_Site_Position'][0])
        non_base_taxa = sorted(str(temp_row['Non_Base_Taxa'][0]).split(sep=";"))

        # If more than max_missing are missing, replace row
        if len(non_base_taxa) > max_missing:
            converted_nonbase = pd.concat([converted_nonbase,temp_row],axis=0, sort=False, ignore_index=True)

        else:
            non_base_index = [j for j,x in enumerate(spp_list) if x in non_base_taxa]

            old_seq = pruned_alignment[:, temp_pos]

            # Get species list and sequence data from taxa with appropriate bases
            new_seq = [v for k,v in enumerate(old_seq) if k not in non_base_index]
            new_spp = [v for k,v in enumerate(spp_list) if k not in non_base_index]

            # Get different alleles as list
            site_set = Counter(new_seq)
            allele_list = list(site_set)
            allele_count = len(allele_list)
            site_counts = list(site_set.values())

            # Ignore singletons or invariant sites
            
            if len(site_counts)==1:
                temp_row.loc[0,'Site_Pattern'] = "non_base_invariant"
                converted_nonbase = pd.concat([converted_nonbase,temp_row],axis=0, sort=False, ignore_index=True)

            elif 1 in site_counts:
                temp_row.loc[0,'Site_Pattern'] = "non_base_singleton"
                converted_nonbase = pd.concat([converted_nonbase,temp_row],axis=0, sort=False, ignore_index=True)
                
            else:
                if len(site_counts)==2:
                    split_1_base = sorted(list(set(new_seq)))[0]
                    base_1_index = findOccurrences(new_seq, split_1_base)
                    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(new_spp)))
                    
                    split_2_base = sorted(list(set(new_seq)))[1]
                    base_2_index = findOccurrences(new_seq, split_2_base)
                    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(new_spp)))
                    
                    splits = sorted([base_1_taxa,base_2_taxa])
                    
                    temp_row.loc[0,'Split_1'] = splits[0]
                    temp_row.loc[0,'Split_2'] = splits[1]
                    
                    if "-" in allele_list:
                        temp_row.loc[0,'Site_Pattern'] = "non_base_gap_biallelic"
                    else:
                        temp_row.loc[0,'Site_Pattern'] = "non_base_biallelic"

                    converted_nonbase = pd.concat([converted_nonbase,temp_row],axis=0, sort=False, ignore_index=True)

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

                    splits = sorted([base_1_taxa,base_2_taxa,base_3_taxa])
                                        
                    temp_row.loc[0,'Split_1'] = splits[0]
                    temp_row.loc[0,'Split_2'] = splits[1]
                    temp_row.loc[0,'Split_3'] = splits[2]

                    if "-" in allele_list:
                        temp_row.loc[0,'Site_Pattern'] = "non_base_gap_triallelic"
                    else:
                        temp_row.loc[0,'Site_Pattern'] = "non_base_triallelic"
                    
                    converted_nonbase = pd.concat([converted_nonbase,temp_row],axis=0, sort=False, ignore_index=True)
                    
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
                                        
                    splits = sorted([base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa])
                    
                    temp_row.loc[0,'Split_1'] = splits[0]
                    temp_row.loc[0,'Split_2'] = splits[1]
                    temp_row.loc[0,'Split_3'] = splits[2]
                    temp_row.loc[0,'Split_4'] = splits[3]

                    if "-" in allele_list:
                        temp_row.loc[0,'Site_Pattern'] = "non_base_gap_quadallelic"
                    else:
                        temp_row.loc[0,'Site_Pattern'] = "non_base_quadallelic"

                    converted_nonbase = pd.concat([converted_nonbase,temp_row],axis=0, sort=False, ignore_index=True)

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
                    
                    splits = sorted([base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,base_5_taxa])

                    temp_row.loc[0,'Split_1'] = splits[0]
                    temp_row.loc[0,'Split_2'] = splits[1]
                    temp_row.loc[0,'Split_3'] = splits[2]
                    temp_row.loc[0,'Split_4'] = splits[3]
                    temp_row.loc[0,'Split_5'] = splits[4]
                    temp_row.loc[0,'Site_Pattern'] = "non_base_pentallelic"
                        
                    converted_nonbase = pd.concat([converted_nonbase,temp_row],axis=0, sort=False, ignore_index=True)

    
    converted_nonbase = converted_nonbase.sort_values(by=['Zeroed_Site_Position']).reindex(columns=final_columns).reset_index().drop('index',axis=1)
    
    return(converted_nonbase)

def process_singletons(pos):

    global pruned_alignment
    global spp_list

    seq_string = pruned_alignment[:, pos]

    # Get column sequence as string
    site_set = dict(Counter(seq_string))

    # https://stackoverflow.com/questions/8023306/get-key-by-value-in-dictionary
    singleton_list = [base for base, count in site_set.items() if count == 1]
    singelton_index = findOccurrences(seq_string, singleton_list)
    singleton_taxa = itemgetter(*singelton_index)(sorted(spp_list))

    if len(singelton_index) == 1:
        return(singleton_taxa)
    if len(singelton_index) > 1:
        return(";".join(sorted(singleton_taxa)))

def process_biallelic(pos):

    global pruned_alignment
    global spp_list

    seq_string = pruned_alignment[:, pos]

    # Get column sequence as string
    split_1_base = sorted(list(set(seq_string)))[0]
    split_2_base = sorted(list(set(seq_string)))[1]

    base_1_index = findOccurrences(seq_string, split_1_base)
    base_1_taxa = ";".join(itemgetter(*base_1_index)(sorted(spp_list)))

    base_2_index = findOccurrences(seq_string, split_2_base)
    base_2_taxa = ";".join(itemgetter(*base_2_index)(sorted(spp_list)))
    return(sorted([base_1_taxa, base_2_taxa]))

def process_triallelic(pos):

    global pruned_alignment
    global spp_list

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

    return(sorted([base_1_taxa, base_2_taxa, base_3_taxa]))

def process_quadallelic(pos):

    global pruned_alignment
    global spp_list

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

    return(sorted([base_1_taxa, base_2_taxa, base_3_taxa, base_4_taxa]))

def process_pentallelic(pos):

    global pruned_alignment
    global spp_list

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

    return(sorted([base_1_taxa, base_2_taxa, base_3_taxa, base_4_taxa, base_5_taxa]))

def getSplitSupport(alignment_path,info_gap,spp_string,max_missing):

    alignment_path = str(alignment_path)
    info_gap = str(info_gap)

    global spp_list
    spp_list = sorted(str(spp_string).split(";"))
    
    max_missing = int(max_missing)
    # Get site split info, or False if alignment cannot be processed for given species set
    all_sites_df = splitMain(alignment_path,info_gap,spp_list,max_missing)

    if all_sites_df.shape[0] == 0:
        return -1
    else:
        return(all_sites_df)

def splitMain(alignment_path,info_gap,spp_list,max_missing):

    # Create dummy return for errors
    empty = pd.DataFrame()

    # If alignment_filename contains all species from spp_list (>= 3 species), continue
    if(subset_alignment(alignment_path, spp_list)):

        # Define each position in the alignment by its 0-based position
        alignment_positions = range(0, pruned_alignment.get_alignment_length())

        # Detect CPUs and create a pool of workers
        if(mp.cpu_count() == 1 or mp.cpu_count() == 2):
            pool_cpu = 1
        else:
            pool_cpu = mp.cpu_count() - 1

        pool = mp.Pool(processes=pool_cpu)

        global site_dict

        # For each position, grab bases and describe variation pattern
        if info_gap == "0":
            site_dict = pool.map(create_site_dict, alignment_positions)
            pool.close()

        elif info_gap == "1":
            site_dict = pool.map(create_site_dict_gap, alignment_positions)
            pool.close()

        else:
            return(empty)

        # Process site patterns
        pattern_list = []

        for i in range(0, len(site_dict)):
            pattern_list.append(
                [i, site_dict[i][i]])

        # Create dataframe of each position (0-based) and it's variation pattern
        site_df = pd.DataFrame(pattern_list, columns=['Zeroed_Site_Position', 'Site_Pattern'])
        site_count = pruned_alignment.get_alignment_length()
        final_columns = ['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']

        invariant_count = 0
        non_base_count = 0
        singleton_count = 0
        biallelic_count = 0
        triallelic_count = 0
        quadallelic_count = 0
        pentallelic_count = 0

        # Process columns with invariant sites
        invariant_df = site_df[site_df['Site_Pattern'].isin(['invariant', 'gap_invariant'])].reset_index(drop=True)

        if invariant_df.shape[0] > 0:

            # Add dummy columns
            invariant_df['Non_Base_Taxa'] = np.nan
            invariant_df['Singleton_Taxa'] = np.nan
            invariant_df['Split_1'] = np.nan
            invariant_df['Split_2'] = np.nan
            invariant_df['Split_3'] = np.nan
            invariant_df['Split_4'] = np.nan
            invariant_df['Split_5'] = np.nan

        else:
            invariant_df = pd.DataFrame(columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'])

        invariant_count = invariant_df.shape[0]

        # Process columns for non-base sites
        non_base_df = site_df[site_df['Site_Pattern'].isin(['non_base'])].reset_index(drop=True)

        if non_base_df.shape[0] > 0:
            pool = mp.Pool(processes=pool_cpu)

            if info_gap == "0":
                non_base_taxa = pool.map(
                    process_nonbase, non_base_df['Zeroed_Site_Position'])
            if info_gap == "1":
                non_base_taxa = pool.map(
                    process_nonbase_gap, non_base_df['Zeroed_Site_Position'])

            pool.close()

            non_base_df['Non_Base_Taxa'] = non_base_taxa

            # Add dummy columns
            non_base_df['Singleton_Taxa'] = np.nan
            non_base_df['Split_1'] = np.nan
            non_base_df['Split_2'] = np.nan
            non_base_df['Split_3'] = np.nan
            non_base_df['Split_4'] = np.nan
            non_base_df['Split_5'] = np.nan
            
            # If missing taxa allowed, reprocess nonbase data
            if max_missing > 0:
                
                if (len(spp_list) - max_missing) < 4:
                    print("Can't allow so many missing taxa that total < 4. Running with max allowable missing")
                    non_base_df = reprocess_nonbase(non_base_df,len(spp_list)-4)
                
                else:
                    non_base_df = reprocess_nonbase(non_base_df,max_missing)


        else:
            non_base_df = pd.DataFrame(columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'])

        non_base_count = non_base_df.shape[0]

        # Process columns for singleton sites
        singleton_df = site_df[site_df['Site_Pattern'].isin(['singleton','gap_singleton'])].reset_index(drop=True)

        if singleton_df.shape[0] > 0:
            pool = mp.Pool(processes=pool_cpu)
            singleton_taxa = pool.map(process_singletons, singleton_df['Zeroed_Site_Position'])
            pool.close()
            singleton_df['Singleton_Taxa'] = singleton_taxa

            # Add dummy columns
            singleton_df['Non_Base_Taxa'] = np.nan
            singleton_df['Split_1'] = np.nan
            singleton_df['Split_2'] = np.nan
            singleton_df['Split_3'] = np.nan
            singleton_df['Split_4'] = np.nan
            singleton_df['Split_5'] = np.nan

        else:
            singleton_df = pd.DataFrame(columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'])

        singleton_count = singleton_df.shape[0]

        # Process columns for biallelic sites
        biallelic_df = site_df[site_df['Site_Pattern'].isin(['biallelic','gap_biallelic'])].reset_index(drop=True)

        if biallelic_df.shape[0] > 0:
            pool = mp.Pool(processes=pool_cpu)
            biallelic_sets = pool.map(process_biallelic, biallelic_df['Zeroed_Site_Position'])
            pool.close()
            biallelic_df = biallelic_df.join(pd.DataFrame(biallelic_sets,columns=['Split_1','Split_2']))

            # Add dummy columns
            biallelic_df['Non_Base_Taxa'] = np.nan
            biallelic_df['Singleton_Taxa'] = np.nan
            biallelic_df['Split_3'] = np.nan
            biallelic_df['Split_4'] = np.nan
            biallelic_df['Split_5'] = np.nan

        else:
            biallelic_df = pd.DataFrame(columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'])

        biallelic_count = biallelic_df.shape[0]

        # Process columns for triallelic sites
        triallelic_df = site_df[site_df['Site_Pattern'].isin(['triallelic','gap_triallelic'])].reset_index(drop=True)

        if triallelic_df.shape[0] > 0:
            pool = mp.Pool(processes=pool_cpu)
            triallelic_sets = pool.map(process_triallelic, triallelic_df['Zeroed_Site_Position'])
            pool.close()
            triallelic_df = triallelic_df.join(pd.DataFrame(triallelic_sets,columns=['Split_1','Split_2','Split_3']))

            # Add dummy columns
            triallelic_df['Non_Base_Taxa'] = np.nan
            triallelic_df['Singleton_Taxa'] = np.nan
            triallelic_df['Split_4'] = np.nan
            triallelic_df['Split_5'] = np.nan

        else:
            triallelic_df = pd.DataFrame(columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'])

        triallelic_count = triallelic_df.shape[0]

        # Process columns for quadallelic sites
        quadallelic_df = site_df[site_df['Site_Pattern'].isin(['quadallelic','gap_quadallelic'])].reset_index(drop=True)

        if quadallelic_df.shape[0] > 0:
            pool = mp.Pool(processes=pool_cpu)
            quadallelic_sets = pool.map(process_quadallelic, quadallelic_df['Zeroed_Site_Position'])
            pool.close()
            quadallelic_df = quadallelic_df.join(pd.DataFrame(quadallelic_sets,columns=['Split_1','Split_2','Split_3','Split_4']))

            # Add dummy columns
            quadallelic_df['Non_Base_Taxa'] = np.nan
            quadallelic_df['Singleton_Taxa'] = np.nan
            quadallelic_df['Split_5'] = np.nan

        else:
            quadallelic_df = pd.DataFrame(columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'])

        quadallelic_count = quadallelic_df.shape[0]

        # Process columns for pentallelic sites
        pentallelic_df = site_df[site_df['Site_Pattern'].isin(['pentallelic','gap_pentallelic'])].reset_index(drop=True)

        if pentallelic_df.shape[0] > 0:
            pool = mp.Pool(processes=pool_cpu)
            pentallelic_sets = pool.map(process_pentallelic, pentallelic_df['Zeroed_Site_Position'])
            pool.close()
            pentallelic_df = pentallelic_df.join(pd.DataFrame(pentallelic_sets,columns=['Split_1','Split_2','Split_3','Split_4','Split_5']))

            # Add dummy columns
            pentallelic_df['Non_Base_Taxa'] = np.nan
            pentallelic_df['Singleton_Taxa'] = np.nan

        else:
            pentallelic_df = pd.DataFrame(columns=['Zeroed_Site_Position','Site_Pattern','Non_Base_Taxa','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'])

        pentallelic_count = pentallelic_df.shape[0]

        # Combine dataframes

        if (invariant_count + non_base_count + singleton_count + biallelic_count + triallelic_count + quadallelic_count + pentallelic_count) > 1:
            all_sites_df = pd.concat([invariant_df,non_base_df,singleton_df,biallelic_df,triallelic_df,quadallelic_df,pentallelic_df,quadallelic_df], axis=0, sort=True).sort_values(by=['Zeroed_Site_Position'])
            all_sites_df = all_sites_df.reindex(columns=final_columns).reset_index()
            return(all_sites_df)

        else:
            return empty
    else:
        return empty
