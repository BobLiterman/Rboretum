.. _getSpeciesComposition:

##########################
**getSpeciesComposition**
##########################

*getSpeciesComposition* takes one or more alignment files and returns information about the sequence content by species, including base counts, 'N' counts, and gap counts

=======================
Function and Arguments
=======================

**Usage**:
::

  getSpeciesComposition <- function(alignment_path,species_info,alignment_name,prefix,suffix)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**alignment_path**				        An absolute or relative path to an alignment file or a directory containing multiple alignments
**species_info**                  A list of species to consider, provided as a character vector or labels or a semicolon-delimited string [Default: all species]
**alignment_name**                A character vector containing desired alignment IDs [Default: Derive name from filename]
**prefix**                        If **alignment_path** points to a directory, only query files starting with **prefix** (e.g. 'Alignment') [Default: Use all files in directory]
**suffix**                        If **alignment_path** points to a directory, only query files ending wtih **suffix** (e.g. '.nex') [Default: Use all files in directory]
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getSpeciesComposition.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()
  
  # Set path to alignment data
  myAlignmentFile <- rb_align1_path
  myAlignmentDir <- rb_alignment_dir

  # Get species composition information
  getSpeciesComposition(alignment_path = myAlignmentFile)

  Taxon Total_Bases Total_N Total_Gaps Percent_GC Percent_N Percent_Gap Alignment_Name
  1  Species_A        2000       0          0     0.4780         0           0  Gene_1.phylip
  2  Species_B        2000       0          0     0.4935         0           0  Gene_1.phylip
  3  Species_C        2000       0          0     0.4940         0           0  Gene_1.phylip
  4  Species_D        2000       0          0     0.4880         0           0  Gene_1.phylip
  5  Species_E        2000       0          0     0.4925         0           0  Gene_1.phylip
  6  Species_F        2000       0          0     0.4820         0           0  Gene_1.phylip
  7  Species_G        2000       0          0     0.4905         0           0  Gene_1.phylip
  8  Species_H        2000       0          0     0.4950         0           0  Gene_1.phylip
  9  Species_I        2000       0          0     0.5035         0           0  Gene_1.phylip
  10 Species_J        2000       0          0     0.5100         0           0  Gene_1.phylip
  11 Species_K        2000       0          0     0.5025         0           0  Gene_1.phylip
  12 Species_L        2000       0          0     0.4950         0           0  Gene_1.phylip
  13 Species_M        2000       0          0     0.5005         0           0  Gene_1.phylip
  14 Species_N        2000       0          0     0.4970         0           0  Gene_1.phylip
  15 Species_O        2000       0          0     0.4885         0           0  Gene_1.phylip
  
  # Get species composition information from all .phylip files in a directory, providing new names
  getSpeciesComposition(alignment_path = myAlignmentDir,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))
  
         Taxon Total_Bases Total_N Total_Gaps Percent_GC Percent_N Percent_Gap Alignment_Name
  1  Species_A        2000       0          0  0.4780000         0           0         Gene_A
  2  Species_B        2000       0          0  0.4935000         0           0         Gene_A
  3  Species_C        2000       0          0  0.4940000         0           0         Gene_A
  4  Species_D        2000       0          0  0.4880000         0           0         Gene_A
  5  Species_E        2000       0          0  0.4925000         0           0         Gene_A
  6  Species_F        2000       0          0  0.4820000         0           0         Gene_A
  7  Species_G        2000       0          0  0.4905000         0           0         Gene_A
  8  Species_H        2000       0          0  0.4950000         0           0         Gene_A
  9  Species_I        2000       0          0  0.5035000         0           0         Gene_A
  10 Species_J        2000       0          0  0.5100000         0           0         Gene_A
  11 Species_K        2000       0          0  0.5025000         0           0         Gene_A
  12 Species_L        2000       0          0  0.4950000         0           0         Gene_A
  13 Species_M        2000       0          0  0.5005000         0           0         Gene_A
  14 Species_N        2000       0          0  0.4970000         0           0         Gene_A
  15 Species_O        2000       0          0  0.4885000         0           0         Gene_A
  16 Species_A        2500       0          0  0.4956000         0           0         Gene_B
  17 Species_B        2500       0          0  0.4840000         0           0         Gene_B
  18 Species_C        2500       0          0  0.4896000         0           0         Gene_B
  .
  .
  .
  
  
  # Get species composition from dummy alignment
  getSpeciesComposition(alignment_path = rb_dummy_align_path)
  
  # Dummy alignment
  
  >Species_A_NoGap_NoN_NoGC
  ATATATATATATATATATATATATATATATATATAT
  >Species_B_10Gap_NoN_NoGC
  ATATATATATATATATATATATATAT----------
  >Species_C_10Gap_10N_NoGC
  ATATATATATATATATNNNNNNNNNN----------
  >Species_D_NoGap_NoN_50GC
  ATATATATATATATATATGCGCGCGCGCGCGCGCGC
  >Species_E_10Gap_NoN_50GC
  ATATATATATATA----------CGCGCGCGCGCGC
  >Species_F_10Gap_10N_50GC
  ATATATATNNNNN----------NNNNNGCGCGCGC

  Taxon Total_Bases Total_N Total_Gaps Percent_GC Percent_N Percent_Gap Alignment_Name
  1 Species_A_NoGap_NoN_NoGC          36       0          0        0.0 0.0000000   0.0000000    Gap_GC_N.fa
  2 Species_B_10Gap_NoN_NoGC          26       0         10        0.0 0.0000000   0.2777778    Gap_GC_N.fa
  3 Species_C_10Gap_10N_NoGC          16      10         10        0.0 0.2777778   0.2777778    Gap_GC_N.fa
  4 Species_D_NoGap_NoN_50GC          36       0          0        0.5 0.0000000   0.0000000    Gap_GC_N.fa
  5 Species_E_10Gap_NoN_50GC          26       0         10        0.5 0.0000000   0.2777778    Gap_GC_N.fa
  6 Species_F_10Gap_10N_50GC          16      10         10        0.5 0.2777778   0.2777778    Gap_GC_N.fa
    
