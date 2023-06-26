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
import copy
import pickle
import tempfile
from functools import partial
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def findOccurrences(s, ch):
    # findOccurrences returns indexes of occurrences of list items in a list or string
    return [i for i, letter in enumerate(s) if letter in ch]

def windowsPool(pool_cpu,pos_list,alignment,spp,base):
    pool = mp.Pool(pool_cpu)
    var_getSiteSplits=partial(getSiteSplits_Windows, pruned_alignment = alignment,bases = base, spp_list = spp)
    split_results = pool.map(var_getSiteSplits, pos_list)
    pool.close()
    return(split_results)

def getSiteSplits_Windows(pos,pruned_alignment,bases,spp_list):
    
    temp_list = copy.deepcopy(spp_list)
    seq_string = pruned_alignment[:, pos]
            
    # Get alleles and their occurrences
    site_set = Counter(seq_string)
    allele_list = list(site_set)
    seq_bases = ";".join(allele_list)

    allele_count = len(allele_list)
    allele_counts = list(site_set.values())
    
    # Get gap taxa
    if "-" in allele_list:
        gap_index = findOccurrences(seq_string, "-")
        gap_taxa = itemgetter(*gap_index)(temp_list)
        gap_count = int(len(gap_index))

        if gap_count == 1:
            gap_taxa_string = gap_taxa
        if gap_count > 1:
            gap_taxa_string = ";".join(sorted(gap_taxa))
    
    else:
        gap_taxa_string = np.nan

    # If column contains any non-bases, remove them and reassess pattern unless fewer than three bases remain
    if(not set(allele_list).issubset(set(bases))):
        non_base_list = list(set(seq_string) - set(bases))
        non_base_index = findOccurrences(seq_string, non_base_list)
        non_base_taxa = itemgetter(*non_base_index)(temp_list)
        non_base_count = int(len(non_base_index))

        if non_base_count == 1:
            non_base_taxa_string = non_base_taxa
        if non_base_count > 1:
            non_base_taxa_string = ";".join(sorted(non_base_taxa))
        
        # If fewer than three bases remain after removing non-bases, return 'non_base'result
        if len(temp_list) - non_base_count < 3:
            return(pd.DataFrame([[pos,seq_bases,'non_base',non_base_taxa_string,non_base_count,np.nan,np.nan,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
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
        non_base_taxa_string = np.nan
        non_base_count = 0

    # If column contains any singletons, remove them and reassess pattern unless fewer than three bases remain
    if(1 in allele_counts):
        singleton_list = [base for base, count in site_set.items() if count == 1]
        singelton_index = findOccurrences(seq_string, singleton_list)
        singleton_taxa = itemgetter(*singelton_index)(temp_list)
        singleton_count = int(len(singelton_index))

        if singleton_count == 1:
            singleton_taxa_string = singleton_taxa
        if singleton_count > 1:
            singleton_taxa_string = ";".join(sorted(singleton_taxa))

        # If fewer than three bases remain after removing singletons, return 'non_base'result
        if len(temp_list) - singleton_count < 3:
            return(pd.DataFrame([[pos,seq_bases,'non_base',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
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
        singleton_taxa_string = np.nan
        singleton_count = 0

    # If only a single allele is present, return 'singleton' or 'invariant' result
    if(allele_count == 1):
        if singleton_count > 0:
            return(pd.DataFrame([[pos,seq_bases,'singleton',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        else:
            return(pd.DataFrame([[pos,seq_bases,'invariant',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
    
    # Process biallelic sites
    if(allele_count == 2):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'biallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    # Process triallelic sites
    if(allele_count == 3):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        split_3_base = sorted(allele_list)[2]
        base_3_index = findOccurrences(seq_string, split_3_base)
        base_3_taxa = ";".join(itemgetter(*base_3_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'triallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    # Process quadallelic sites
    if(allele_count == 4):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        split_3_base = sorted(allele_list)[2]
        base_3_index = findOccurrences(seq_string, split_3_base)
        base_3_taxa = ";".join(itemgetter(*base_3_index)(temp_list))
        
        split_4_base = sorted(allele_list)[3]
        base_4_index = findOccurrences(seq_string, split_4_base)
        base_4_taxa = ";".join(itemgetter(*base_4_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'quadallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    # Process pentallelic sites
    if(allele_count == 5):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        split_3_base = sorted(allele_list)[2]
        base_3_index = findOccurrences(seq_string, split_3_base)
        base_3_taxa = ";".join(itemgetter(*base_3_index)(temp_list))
        
        split_4_base = sorted(allele_list)[3]
        base_4_index = findOccurrences(seq_string, split_4_base)
        base_4_taxa = ";".join(itemgetter(*base_4_index)(temp_list))
        
        split_5_base = sorted(allele_list)[4]
        base_5_index = findOccurrences(seq_string, split_5_base)
        base_5_taxa = ";".join(itemgetter(*base_5_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'pentallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,base_5_taxa]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def getSiteSplits(pos):
    
    global pruned_alignment
    global spp_list
    global bases
    
    temp_list = copy.deepcopy(spp_list)
    seq_string = pruned_alignment[:, pos]
            
    # Get alleles and their occurrences
    site_set = Counter(seq_string)
    allele_list = list(site_set)
    seq_bases = ";".join(allele_list)

    allele_count = len(allele_list)
    allele_counts = list(site_set.values())
    
    # Get gap taxa
    if "-" in allele_list:
        gap_index = findOccurrences(seq_string, "-")
        gap_taxa = itemgetter(*gap_index)(temp_list)
        gap_count = int(len(gap_index))

        if gap_count == 1:
            gap_taxa_string = gap_taxa
        if gap_count > 1:
            gap_taxa_string = ";".join(sorted(gap_taxa))
    
    else:
        gap_taxa_string = np.nan

    # If column contains any non-bases, remove them and reassess pattern unless fewer than three bases remain
    if(not set(allele_list).issubset(set(bases))):
        non_base_list = list(set(seq_string) - set(bases))
        non_base_index = findOccurrences(seq_string, non_base_list)
        non_base_taxa = itemgetter(*non_base_index)(temp_list)
        non_base_count = int(len(non_base_index))

        if non_base_count == 1:
            non_base_taxa_string = non_base_taxa
        if non_base_count > 1:
            non_base_taxa_string = ";".join(sorted(non_base_taxa))
        
        # If fewer than three bases remain after removing non-bases, return 'non_base'result
        if len(temp_list) - non_base_count < 3:
            return(pd.DataFrame([[pos,seq_bases,'non_base',non_base_taxa_string,non_base_count,np.nan,np.nan,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
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
        non_base_taxa_string = np.nan
        non_base_count = 0

    # If column contains any singletons, remove them and reassess pattern unless fewer than three bases remain
    if(1 in allele_counts):
        singleton_list = [base for base, count in site_set.items() if count == 1]
        singelton_index = findOccurrences(seq_string, singleton_list)
        singleton_taxa = itemgetter(*singelton_index)(temp_list)
        singleton_count = int(len(singelton_index))

        if singleton_count == 1:
            singleton_taxa_string = singleton_taxa
        if singleton_count > 1:
            singleton_taxa_string = ";".join(sorted(singleton_taxa))

        # If fewer than three bases remain after removing singletons, return 'non_base'result
        if len(temp_list) - singleton_count < 3:
            return(pd.DataFrame([[pos,seq_bases,'non_base',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        
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
        singleton_taxa_string = np.nan
        singleton_count = 0

    # If only a single allele is present, return 'singleton' or 'invariant' result
    if(allele_count == 1):
        if singleton_count > 0:
            return(pd.DataFrame([[pos,seq_bases,'singleton',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
        else:
            return(pd.DataFrame([[pos,seq_bases,'invariant',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,np.nan,np.nan,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))
    
    # Process biallelic sites
    if(allele_count == 2):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'biallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,np.nan,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    # Process triallelic sites
    if(allele_count == 3):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        split_3_base = sorted(allele_list)[2]
        base_3_index = findOccurrences(seq_string, split_3_base)
        base_3_taxa = ";".join(itemgetter(*base_3_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'triallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,base_3_taxa,np.nan,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    # Process quadallelic sites
    if(allele_count == 4):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        split_3_base = sorted(allele_list)[2]
        base_3_index = findOccurrences(seq_string, split_3_base)
        base_3_taxa = ";".join(itemgetter(*base_3_index)(temp_list))
        
        split_4_base = sorted(allele_list)[3]
        base_4_index = findOccurrences(seq_string, split_4_base)
        base_4_taxa = ";".join(itemgetter(*base_4_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'quadallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,np.nan]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

    # Process pentallelic sites
    if(allele_count == 5):
        # Get alleles and their associated taxa
        split_1_base = sorted(allele_list)[0]
        base_1_index = findOccurrences(seq_string, split_1_base)
        base_1_taxa = ";".join(itemgetter(*base_1_index)(temp_list))
        
        split_2_base = sorted(allele_list)[1]
        base_2_index = findOccurrences(seq_string, split_2_base)
        base_2_taxa = ";".join(itemgetter(*base_2_index)(temp_list))
        
        split_3_base = sorted(allele_list)[2]
        base_3_index = findOccurrences(seq_string, split_3_base)
        base_3_taxa = ";".join(itemgetter(*base_3_index)(temp_list))
        
        split_4_base = sorted(allele_list)[3]
        base_4_index = findOccurrences(seq_string, split_4_base)
        base_4_taxa = ";".join(itemgetter(*base_4_index)(temp_list))
        
        split_5_base = sorted(allele_list)[4]
        base_5_index = findOccurrences(seq_string, split_5_base)
        base_5_taxa = ";".join(itemgetter(*base_5_index)(temp_list))
        
        return(pd.DataFrame([[pos,seq_bases,'pentallelic',non_base_taxa_string,non_base_count,singleton_taxa_string,singleton_count,gap_taxa_string,base_1_taxa,base_2_taxa,base_3_taxa,base_4_taxa,base_5_taxa]], columns=['Alignment_Position','All_Site_Bases','Site_Pattern','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Singleton_Count','Gap_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5']))

def splitsProcessor(path_to_align,spp_info,use_gaps,align_name):

    # Set path to alignment
    alignment_path = str(path_to_align)

    # Get species list from semicolon-separated string
    global spp_list
    spp_list = sorted(str(spp_info).split(";"))

    # Set alignment name
    alignment_name = str(align_name)

    # Set gap data
    info_gap = str(use_gaps)

    # Set valid bases based on info_gap
    global bases
    
    if use_gaps != "0" and use_gaps != "1":
        sys.exit("ERROR: info_gap must be '0' or '1'")

    if use_gaps == "0":
        bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c']

    else:
        bases = ['A', 'C', 'T', 'G', 'a', 'g', 't', 'c','-']

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

    # If fewer than three species exist in alignment or species list, raise exception
    if (len(spp_list) < 3) or (len(raw_spp) < 3):
        sys.exit("ERROR: Cannot process fewer than 3 species...")
        
    # If requested species are not in alignment, raise exception
    if not all(elem in raw_spp for elem in spp_list):
        sys.exit("ERROR: Requested species not found in "+os.path.basename(alignment_path)+"...")

    # Re-order alignment to sorted order of species list
    else:
        # Create dummy alignment  
        global pruned_alignment
        pruned_alignment = raw_alignment[0:0]
        
        # Populate alignment by adding taxa sorted by taxon ID
        for i in spp_list:
            #pruned_alignment.add_sequence(str(raw_alignment[raw_spp.index(i)].id), str(raw_alignment[raw_spp.index(i)].seq))
            pruned_alignment.append(SeqRecord(Seq(str(raw_alignment[raw_spp.index(i)].seq), generic_dna), id=str(raw_alignment[raw_spp.index(i)].id)))

        
        # If resulting alignment is empty, raise exception
        if int(pruned_alignment.get_alignment_length()) == 0:
            sys.exit("ERROR: Alignment processed, but appears to have no bases...")

    # Detect CPUs and create a pool of workers
    if mp.cpu_count() == 1:
        pool_cpu = 1
    else:
        pool_cpu = mp.cpu_count() - 1

    # Get a 0-based coordinate map for alignment
    alignment_positions = range(0, pruned_alignment.get_alignment_length())

    if os.name == 'nt':
        split_results = windowsPool(pool_cpu,alignment_positions,pruned_alignment,spp_list,bases)
    else:
        with mp.Pool(pool_cpu) as pool:
            split_results = pool.map(getSiteSplits, alignment_positions)
    
    split_results = pd.concat(split_results).sort_values(by=['Alignment_Position']).reset_index()
    
    # Add 1 to site position to convert from 0 to 1-based
    split_results.Alignment_Position = split_results.Alignment_Position + 1
    split_results['Alignment_Name'] = alignment_name

    temp_df = tempfile.NamedTemporaryFile()
    temp_df_name = temp_df.name
    temp_df.close()
    with open(temp_df_name, "wb") as f:
        pickle.dump(split_results, f)
    f.close()

    print(temp_df_name)

def main(path_to_align,spp_info,use_gaps,align_name):
    splitsProcessor(path_to_align,spp_info,use_gaps,align_name)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
