.. _getAlignmentComposition:

############################
**getAlignmentComposition**
############################

*getAlignmentComposition* takes one or more alignment files and returns information about the sequence content, including base counts, 'N' counts, and gap counts

=======================
Function and Arguments
=======================

**Usage**:
::

  getAlignmentComposition <- function(alignment_path,species_info,alignment_name,prefix,suffix)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**alignment_path**				        An absolute or relative path to an alignment file or a directory containing multiple alignments
**species_info**                  A list of species to consider, provided as a phylo, multiPhlyo, character vector or labels or a semicolon-delimited string [Default: all species]
**alignment_name**                A character vector containing desired alignment IDs [Default: Derive name from filename]
**prefix**                        If **alignment_path** points to a directory, only query files starting with **prefix** (e.g. 'Alignment') [Default: Use all files in directory]
**suffix**                        If **alignment_path** points to a directory, only query files ending wtih **suffix** (e.g. '.nex') [Default: Use all files in directory]
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getAlignmentComposition.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()
  
  # Set path to alignment data
  myAlignmentFile <- rb_align1_path
  myAlignmentDir <- rb_alignment_dir

  # Get alignment composition information for a single alignment
  getAlignmentComposition(alignment_path = myAlignmentFile)
  
    Alignment_Name Alignment_Length Percent_GC Percent_N Percent_Gap
  1  Gene_1.phylip             2000  0.4940333         0           0

  # Get alignment composition information from all .phylip files in a directory, providing new names
  getAlignmentComposition(alignment_path = myAlignmentDir,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

    Alignment_Name Alignment_Length Percent_GC Percent_N Percent_Gap
  1         Gene_A             2000  0.4940333         0           0
  2         Gene_B             2500  0.4978400         0           0
  3         Gene_C             1500  0.4921778         0           0
  4         Gene_D             3000  0.5011556         0           0
  5         Gene_E             2500  0.5003733         0           0
  
  # Get species composition from dummy alignment
  getAlignmentComposition(alignment_path = rb_dummy_align_path)
  
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

    Alignment_Name Alignment_Length Percent_GC  Percent_N Percent_Gap
  1    Gap_GC_N.fa               36       0.25 0.09259259   0.1851852
