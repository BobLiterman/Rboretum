.. _summarizeMultiphylo:

########################
**summarizeMultiphylo**
########################

*summarizeMultiphylo* prints a brief, useful summary of a multiPhylo object

=======================
Function and Arguments
=======================

**Usage**:
::

  summarizeMultiphylo <- function(trees)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**trees**				                  A rooted multiPhylo object where all trees share 3+ taxa
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/summarizeMultiphylo.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()

  # Read in trees with two topologies (1=2, 1!=3)

  myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
  myTree_2 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))
  myTrees <- c(myTree_1,myTree_2)

  summarizeMultiPhylo(myTrees)

  [1] "Trees unnamed, adding dummy names..."
  [1] "Command: treeNamer(trees)"
  [1] "Read in 2 trees, that all share 15 tip labels..."
  [1] "All trees contain identical tip labels..."
  [1] "Among all trees, including root splits there are 17 unique monophyletic clades..."
  [1] "Command: getTreeClades(trees,return_counts = TRUE)"
  [1] "All trees share 9 clades, including:"
  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_A;Species_F"                                                                                                              
  [3] "Species_B;Species_O"                                                                                                              
  [4] "Species_C;Species_H"                                                                                                              
  [5] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                                                
  [6] "Species_G;Species_I;Species_N"                                                                                                    
  [7] "Species_G;Species_N"                                                                                                              
  [8] "Species_J;Species_M"                                                                                                              
  [9] "Species_K;Species_L"                                                                                                              
  [1] "Command: getTreeClades(trees,return_shared = TRUE)"
  [1] "All trees have a unique topology..."
  
