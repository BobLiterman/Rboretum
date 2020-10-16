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
  
  # Dummy Alignment

  >Species_A
  NAAAAAAAAAAAA
  >Species_B
  NAAAAAAAAAAAA
  >Species_C
  NAAAAAAAATAAA
  >Species_D
  NNAAAAATTTAAA
  >Species_E
  NNAAAAATTCAAA
  >Species_F
  NNAAAATTTCTTW
  >Species_G
  NNAAA-TCCGTTK
  >Species_H
  NNAAAGTCCGTTM
  >Species_I
  NNAAACTCG-TTS
  >Species_J
  NNAT-TTCG--GR
  
  getAlignmentComposition(alignment_path = rb_dummy_align_path)

        Alignment_Name Alignment_Length Percent_GC Percent_Degenerate Percent_N Percent_Gap
  1 Dummy_Alignment.fa               13  0.1456311         0.03846154 0.1307692  0.03846154
  
