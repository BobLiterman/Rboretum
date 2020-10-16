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
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getSpeciesComposition.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()
  
  # Set path to alignment data
  myAlignmentFile <- rb_align1_path
  myAlignmentDir <- rb_alignment_dir

  # Get species composition information
  getSpeciesComposition(alignment_path = myAlignmentFile)

         Taxon Total_Bases Total_Degenerate Total_N Total_Gaps Percent_GC Percent_Degenerate Percent_N Percent_Gap Alignment_Name
  1  Species_A        2000                0       0          0     0.4950                  0         0           0  Gene_1.phylip
  2  Species_B        2000                0       0          0     0.4940                  0         0           0  Gene_1.phylip
  3  Species_C        2000                0       0          0     0.5025                  0         0           0  Gene_1.phylip
  4  Species_D        2000                0       0          0     0.4950                  0         0           0  Gene_1.phylip
  5  Species_E        2000                0       0          0     0.4780                  0         0           0  Gene_1.phylip
  6  Species_F        2000                0       0          0     0.4820                  0         0           0  Gene_1.phylip
  7  Species_G        2000                0       0          0     0.4925                  0         0           0  Gene_1.phylip
  8  Species_H        2000                0       0          0     0.5005                  0         0           0  Gene_1.phylip
  9  Species_I        2000                0       0          0     0.5100                  0         0           0  Gene_1.phylip
  10 Species_J        2000                0       0          0     0.5035                  0         0           0  Gene_1.phylip
  11 Species_K        2000                0       0          0     0.4970                  0         0           0  Gene_1.phylip
  12 Species_L        2000                0       0          0     0.4905                  0         0           0  Gene_1.phylip
  13 Species_M        2000                0       0          0     0.4935                  0         0           0  Gene_1.phylip
  14 Species_N        2000                0       0          0     0.4885                  0         0           0  Gene_1.phylip
  15 Species_O        2000                0       0          0     0.4880                  0         0           0  Gene_1.phylip
  
  # Get species composition information from all .phylip files in a directory, providing new names
  getSpeciesComposition(alignment_path = myAlignmentDir,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))
  
         Taxon Total_Bases Total_Degenerate Total_N Total_Gaps Percent_GC Percent_Degenerate Percent_N Percent_Gap Alignment_Name
  1  Species_A        2000                0       0          0  0.4950000                  0         0           0         Gene_A
  2  Species_B        2000                0       0          0  0.4940000                  0         0           0         Gene_A
  3  Species_C        2000                0       0          0  0.5025000                  0         0           0         Gene_A
  4  Species_D        2000                0       0          0  0.4950000                  0         0           0         Gene_A
  5  Species_E        2000                0       0          0  0.4780000                  0         0           0         Gene_A
  6  Species_F        2000                0       0          0  0.4820000                  0         0           0         Gene_A
  7  Species_G        2000                0       0          0  0.4925000                  0         0           0         Gene_A
  8  Species_H        2000                0       0          0  0.5005000                  0         0           0         Gene_A
  9  Species_I        2000                0       0          0  0.5100000                  0         0           0         Gene_A
  10 Species_J        2000                0       0          0  0.5035000                  0         0           0         Gene_A
  11 Species_K        2000                0       0          0  0.4970000                  0         0           0         Gene_A
  12 Species_L        2000                0       0          0  0.4905000                  0         0           0         Gene_A
  13 Species_M        2000                0       0          0  0.4935000                  0         0           0         Gene_A
  14 Species_N        2000                0       0          0  0.4885000                  0         0           0         Gene_A
  15 Species_O        2000                0       0          0  0.4880000                  0         0           0         Gene_A
  16 Species_A        2500                0       0          0  0.4916000                  0         0           0         Gene_B
  17 Species_B        2500                0       0          0  0.4896000                  0         0           0         Gene_B
  18 Species_C        2500                0       0          0  0.4904000                  0         0           0         Gene_B
  19 Species_D        2500                0       0          0  0.4920000                  0         0           0         Gene_B
  .
  .
  .
  
  # Get species composition from dummy alignment
  getSpeciesComposition(alignment_path = rb_dummy_align_path)
  
         Taxon Total_Bases Total_Degenerate Total_N Total_Gaps Percent_GC Percent_Degenerate  Percent_N Percent_Gap     Alignment_Name
  1  Species_A          12                0       1          0 0.00000000         0.00000000 0.07692308  0.00000000 Dummy_Alignment.fa
  2  Species_B          12                0       1          0 0.00000000         0.00000000 0.07692308  0.00000000 Dummy_Alignment.fa
  3  Species_C          11                0       2          0 0.00000000         0.00000000 0.15384615  0.00000000 Dummy_Alignment.fa
  4  Species_D          11                0       2          0 0.00000000         0.00000000 0.15384615  0.00000000 Dummy_Alignment.fa
  5  Species_E          11                0       2          0 0.09090909         0.00000000 0.15384615  0.00000000 Dummy_Alignment.fa
  6  Species_F          10                1       2          0 0.10000000         0.07692308 0.15384615  0.00000000 Dummy_Alignment.fa
  7  Species_G           9                1       2          1 0.33333333         0.07692308 0.15384615  0.07692308 Dummy_Alignment.fa
  8  Species_H          10                1       2          0 0.40000000         0.07692308 0.15384615  0.00000000 Dummy_Alignment.fa
  9  Species_I           9                1       2          1 0.33333333         0.07692308 0.15384615  0.07692308 Dummy_Alignment.fa
  10 Species_J           7                1       2          3 0.42857143         0.07692308 0.15384615  0.23076923 Dummy_Alignment.fa

**Dummy Alignment**

.. image:: ../images/Dummy_Align.png
  :width: 600
