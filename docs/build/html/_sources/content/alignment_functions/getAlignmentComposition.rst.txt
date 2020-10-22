.. _getAlignmentComposition:

############################
**getAlignmentComposition**
############################

**getAlignmentComposition** takes one or more alignment files and returns information about the sequence content, including base counts, 'N' counts, and gap counts

.. warning::
  
  **Note:** Rboretum currently does not support handling of degenerate bases, and classifies them broadly as N/missing. Support for degenerate nucleotides is in progress and set for a future release. 

=======================
Function and Arguments
=======================

**Usage**:
::

  getAlignmentComposition <- function(alignment_path,species_info,alignment_name,prefix,suffix)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**alignment_path**				        An absolute or relative path to an alignment file(s) or a directory containing multiple alignments
**species_info**                  A list of species to consider, provided as a phylo, multiPhlyo, character vector or labels or a semicolon-delimited string [Default: Shared species among all alignments]
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
  
    Alignment_Name Alignment_Length Percent_GC Percent_Degenerate  Percent_N Percent_Gap
  1     Gene_1.phy             1551  0.3982911                  0 0.02677842  0.03249516

  # Get alignment composition information from all .phylip files in a directory, providing new names
  getAlignmentComposition(alignment_path = myAlignmentDir,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

    Alignment_Name Alignment_Length Percent_GC Percent_Degenerate  Percent_N Percent_Gap
  1         Gene_A             1551  0.3982911                  0 0.02677842  0.03249516
  2         Gene_B             2804  0.5219202                  0 0.02691393  0.10537328
  3         Gene_C             1031  0.5919662                  0 0.02644682  0.27009376
  4         Gene_D             2219  0.4160601                  0 0.02679886  0.09397627
  5         Gene_E             1500  0.4895705                  0 0.02626667  0.00000000
  
  # Get alignment composition information from all .phy files in a directory, providing new names, considering only Species A - E
  getAlignmentComposition(alignment_path = myAlignmentDir,species_info = 'Species_A;Species_B;Species_C;Species_D;Species_E',suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

    Alignment_Name Alignment_Length Percent_GC Percent_Degenerate  Percent_N Percent_Gap
  1         Gene_A             1551  0.4011516                  0 0.02398453  0.03546099
  2         Gene_B             2804  0.5211614                  0 0.02410842  0.10627675
  3         Gene_C             1031  0.5859697                  0 0.02347236  0.27138700
  4         Gene_D             2219  0.4158223                  0 0.02424516  0.09508788
  5         Gene_E             1500  0.4852580                  0 0.02320000  0.00000000

  # Get species composition from dummy alignment
  getAlignmentComposition(alignment_path = rb_dummy_align_path)

        Alignment_Name Alignment_Length Percent_GC Percent_Degenerate Percent_N Percent_Gap
  1 Dummy_Alignment.fa               13  0.1470588         0.03846154 0.1384615  0.03846154

**Dummy Alignment**

.. image:: ../images/Dummy_Align.png
  :width: 600
