.. _getAlignmentBreakdown:

##########################
**getAlignmentBreakdown**
##########################

*getAlignmentBreakdown* returns the species IDs from an alignment file, specified by *alignment_path*

.. warning::
  
  **Note:** Rboretum currently does not support handling of degenerate bases, and classifies them broadly as missing. Support for degenerate nucleotides is in progress and set for a future release. 

=======================
Function and Arguments
=======================

**Usage**:
::

  getAlignmentBreakdown <- function(alignment_path,species_info,use_gaps,alignment_name,prefix,suffix)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**alignment_path**				        An absolute or relative path to an alignment file(s) or a directory containing multiple alignments
**species_info**                  A list of species to consider, provided as a phylo, multiPhlyo, character vector or labels or a semicolon-delimited string [Default: Shared species among all alignments]
**use_gaps**                      If FALSE, consider gaps (-) in alignments as missing data. [Default: TRUE, treat gaps as indel characters]
**alignment_name**                A character vector containing desired alignment IDs [Default: Derive name from filename]
**prefix**                        If **alignment_path** points to a directory, only query files starting with **prefix** (e.g. 'Alignment') [Default: Use all files in directory]
**suffix**                        If **alignment_path** points to a directory, only query files ending wtih **suffix** (e.g. '.nex') [Default: Use all files in directory]
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getAlignmentBreakdown.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()

  # Set path to alignment data
  myAlignmentFile <- rb_align1_path
  mySpecies <- getAlignmentSpecies(myAlignmentFile)
  myAlignmentDir <- rb_alignment_dir
  
  # Get alignment Breakdown information for a single alignment
  getAlignmentBreakdown(alignment_path = myAlignmentFile)
  
    Alignment_Name Alignment_Length Alignment_Percent_GC Percent_Degenerate Percent_N Percent_Gap Species_GC_Mean Species_GC_StdDev Percent_Nonbase Percent_Invariant Percent_Singleton Percent_Parsimony_Informative  Percent_Biallelic Percent_Triallelic Percent_Quadallelic Percent_Pentallelic
  1  Gene_1.phylip             2000            0.4940333                  0         0           0       0.4940333       0.008249387               0             0.009            0.0355                        0.5945			   0.342             0.4685               0.145            		  0
  
  # Get alignment Breakdown information from all .phylip files in a directory, providing new names, consider gaps as missing data
  getAlignmentBreakdown(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))
  
    Alignment_Name Alignment_Length Alignment_Percent_GC Percent_Degenerate Percent_N Percent_Gap Species_GC_Mean Species_GC_StdDev Percent_Nonbase Percent_Invariant Percent_Singleton Percent_Parsimony_Informative  Percent_Biallelic Percent_Triallelic Percent_Quadallelic Percent_Pentallelic
  1         Gene_A             2000            0.4940333                  0         0           0       0.4940333       0.008249387               0            0.0090            0.0355                     0.5945000		   0.3420000          0.4685000              0.1450                   0
  2         Gene_B             2500            0.4978400                  0         0           0       0.4978400       0.012087112               0            0.0912            0.1028                     0.5788000		   0.5004000          0.2768000              0.0288                   0
  3         Gene_C             1500            0.4921778                  0         0           0       0.4921778       0.010674005               0            0.0120            0.0320                     0.5846667		   0.3653333          0.4586667              0.1320                   0
  4         Gene_D             3000            0.5011556                  0         0           0       0.5011556       0.008174551               0            0.0100            0.0380                     0.6080000		   0.3326667          0.4853333              0.1340                   0
  5         Gene_E             2500            0.5003733                  0         0           0       0.5003733       0.010630916               0            0.0000            0.0048                     0.6144000		   0.1132000          0.5200000              0.3620                   0

  # Get alignment Breakdown from dummy alignment, with and without gap support
  getAlignmentBreakdown(alignment_path = rb_dummy_align_path)
  
        Alignment_Name Alignment_Length Alignment_Percent_GC Percent_Degenerate Percent_N Percent_Gap Species_GC_Mean Species_GC_StdDev Percent_Nonbase Percent_Invariant Percent_Singleton Percent_Parsimony_Informative  Percent_Biallelic Percent_Triallelic Percent_Quadallelic Percent_Pentallelic
  1 Dummy_Alignment.fa               13            0.1470588         0.03846154 0.1384615  0.03846154       0.1686147         0.1825288       0.1538462         0.1538462         0.2307692                     0.3076923          0.2307692         0.07692308          0.07692308          0.07692308
  
  getAlignmentBreakdown(alignment_path = rb_dummy_align_path,use_gaps = FALSE)
  
        Alignment_Name Alignment_Length Alignment_Percent_GC Percent_Degenerate Percent_N Percent_Gap Species_GC_Mean Species_GC_StdDev Percent_Nonbase Percent_Invariant Percent_Singleton Percent_Parsimony_Informative  Percent_Biallelic Percent_Triallelic Percent_Quadallelic Percent_Pentallelic
  1 Dummy_Alignment.fa               13            0.1470588         0.03846154 0.1384615  0.03846154       0.1686147         0.1825288       0.1538462         0.2307692         0.1538462                     0.3846154          0.2307692         0.07692308           0.1538462                   0

**Dummy Alignment**

.. image:: ../images/Dummy_Align.png
  :width: 600
