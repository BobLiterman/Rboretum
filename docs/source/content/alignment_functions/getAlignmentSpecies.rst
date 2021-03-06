.. _getAlignmentSpecies:

########################
**getAlignmentSpecies**
########################

**getAlignmentSpecies** returns the species IDs from an alignment file, specified by *alignment_path*

=======================
Function and Arguments
=======================

**Usage**:
::

  getAlignmentSpecies <- function(alignment_path)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**alignment_path**				        An absolute or relative path to the alignment file from which you want to extract labels
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getAlignmentSpecies.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()

  # Set alignment path
  myAlignmentFile <- rb_align1_path
  
  myAlignmentFile
  [1] "<PACKAGE_DIR>/extdata/alignments/Gene_1.phylip"


  # Get sample IDs from alignment
  getAlignmentSpecies(alignment_path = myAlignmentFile)
  
  [1] "Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  
