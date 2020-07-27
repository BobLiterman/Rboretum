#!/usr/bin/env python3
"""
    This script will take an alignment file (NEXUS, FASTA, or Phylip-relaxed) and:

    1) Prune down to a supplied subset of taxa
    2) Breakdown the signal at each site
    
    NOT PULLING RIGHT TIP LABELS FOR SPLITS!
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
import copy

def patternProcessor(align_path,use_gaps,spp_info,align_name):

    global alignment_path
    global info_gap
    global spp_list
    global bases

    # Prepare arguments
    alignment_path = str(align_path)
    info_gap = str(use_gaps)
    
    # Set valid bases based on info_gap
    if use_gaps != "0" and use_gaps != "1":
        sys.exit("ERROR: info_gap must be '0' or '1'")
    
    if use_gaps == "0":
        bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c']
    
    else:
        bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c','-']
        
    spp_list = sorted((spp_info).split(";"))
    
    # Set alignment name
    global alignment_name
    alignment_name = str(align_name)
    
    # If alignment_filename contains all species from spp_list (>= 3 species), continue
    if not getPrunedAlignment():
        sys.exit("ERROR: Cannot process "+os.path.basename(alignment_path)+" with provided species list.")
    
    # For each position, grab bases and describe variation pattern
    if info_gap == "0":
        print("Processing site split patterns for " + os.path.basename(alignment_path) + ", with gaps interpreted as missing data...")
    if info_gap == "1":
        print("Processing site split patterns for " + os.path.basename(alignment_path) + ", with gaps interpreted as valid bases...")
        
    # Detect CPUs and create a pool of workers
    if mp.cpu_count() == 1:
        pool_cpu = 1
    else:
        pool_cpu = mp.cpu_count() - 1
    
    # Get a 0-based coordinate map for alignment
    alignment_positions = range(0, pruned_alignment.get_alignment_length())

    with mp.Pool(pool_cpu) as pool:
        pattern_results = pool.map(getSitePatterns, alignment_positions)
    
    print("Signal processed...compiling final dataframe...")
    pattern_results = pd.concat(pattern_results).sort_values(by=['Alignment_Position']).reset_index()
    
    # Add 1 to site position to convert from 0 to 1-based
    pattern_results.Alignment_Position = pattern_results.Alignment_Position + 1
    pattern_results['Alignment_Name'] = alignment_name

    return pattern_results

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
    except:
        return False

    # Get alignment species
    raw_spp = []
    for seq_record in raw_alignment:
        raw_spp.append(str(seq_record.id))
    
    # If species list/alignment contains fewer than 3 species, return False
    if len(spp_list) < 3 or len(raw_spp) < 3:
        return False
    
    # If any species from spp_list are not present in raw_spp, return False
    if not all(x in raw_spp for x in spp_list): 
        return False
    
    # Re-order alignment to sorted order of species list
    else:
        # Create dummy alignment        
        global pruned_alignment
        pruned_alignment = raw_alignment[0:0]
        
        # Populate alignment by adding taxa sorted by taxon ID
        for i in spp_list:
            pruned_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))
        
        # Return True to indicate a valid alignment was processed
        return True

def findOccurrences(s, ch):
    # findOccurrences returns indexes of occurrences of list items in a list or string
    return [i for i, letter in enumerate(s) if letter in ch]

def getSitePatterns(pos):
    
    global pruned_alignment
    global bases
    global spp_list
    
    temp_list = copy.deepcopy(spp_list)
    seq_string = pruned_alignment[:, pos]
            
    # Get alleles and their occurrences
    site_set = Counter(seq_string)
    allele_list = list(site_set)
    seq_bases = ";".join(allele_list)

    allele_count = len(allele_list)
    allele_counts = list(site_set.values())
    
    # If column contains any non-bases, remove them and reassess pattern unless fewer than three bases remain
    if(not set(allele_list).issubset(set(bases))):
        non_base_list = list(set(seq_string) - set(bases))
        non_base_index = findOccurrences(seq_string, non_base_list)
        non_base_count = int(len(non_base_index))

        # If fewer than three bases remain after removing non-bases, return 'non_base'result
        if len(temp_list) - non_base_count < 3:
            return(pd.DataFrame([[pos,'non_base',non_base_count,np.nan]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
        
        else:
            # Get species list and sequence data from taxa with appropriate bases
            seq_string = [v for k,v in enumerate(seq_string) if k not in non_base_index]
            temp_list = [v for k,v in enumerate(temp_list) if k not in non_base_index]
            
            # Re-get alleles and their occurrences
            site_set = Counter(seq_string)
            allele_list = list(site_set)

            allele_count = len(allele_list)
            allele_counts = list(site_set.values())

    else:
        non_base_count = 0

    # If column contains any singletons, remove them and reassess pattern unless fewer than three bases remain
    if(1 in allele_counts):
        singleton_list = [base for base, count in site_set.items() if count == 1]
        singelton_index = findOccurrences(seq_string, singleton_list)
        singleton_count = int(len(singelton_index))

        # If fewer than three bases remain after removing singletons, return 'non_base'result
        if len(temp_list) - singleton_count < 3:
            return(pd.DataFrame([[pos,'non_base',non_base_count,np.nan]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
        
        else:
            # Get species list and sequence data from taxa with appropriate bases
            seq_string = [v for k,v in enumerate(seq_string) if k not in singelton_index]
            temp_list = [v for k,v in enumerate(temp_list) if k not in singelton_index]
            
            # Re-get alleles and their occurrences
            site_set = Counter(seq_string)
            allele_list = list(site_set)

            allele_count = len(allele_list)
            allele_counts = list(site_set.values())

    else:
        singleton_count = 0

    # If only a single allele is present, return 'singleton' or 'invariant' result
    if(allele_count == 1):
        if singleton_count > 0:
            return(pd.DataFrame([[pos,'singleton',non_base_count,"No"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
        else:
            return(pd.DataFrame([[pos,'invariant',non_base_count,"No"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
    
    # Process biallelic sites
    if(allele_count == 2):
        if singleton_count > 0:
            return(pd.DataFrame([[pos,'biallelic',non_base_count,"No"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
        else:
            return(pd.DataFrame([[pos,'biallelic',non_base_count,"Yes"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))

    # Process triallelic sites
    if(allele_count == 3):
        if singleton_count > 0:
            return(pd.DataFrame([[pos,'triallelic',non_base_count,"No"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
        else:
            return(pd.DataFrame([[pos,'triallelic',non_base_count,"Yes"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))

    # Process quadallelic sites
    if(allele_count == 4):
        if singleton_count > 0:
            return(pd.DataFrame([[pos,'quadallelic',non_base_count,"No"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
        else:
            return(pd.DataFrame([[pos,'quadallelic',non_base_count,"Yes"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))

    # Process pentallelic sites
    if(allele_count == 5):
        if singleton_count > 0:
            return(pd.DataFrame([[pos,'pentallelic',non_base_count,"No"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))
        else:
            return(pd.DataFrame([[pos,'pentallelic',non_base_count,"Yes"]], columns=['Alignment_Position','Site_Pattern','Non_Base_Count','Parsimony_Informative']))

