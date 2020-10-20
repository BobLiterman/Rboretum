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
  sourceRboretum()
  
  # Read in trees
  myTrees <- readRooted(rb_unroot_dir,root_taxa = c('Species_C','Species_H'))
  summarizeMultiPhylo(myTrees)
  
  [1] "Read in 5 trees, that all share 15 tip labels..."
  [1] "All trees contain identical tip labels..."
  [1] "Among all trees, including root splits there are 21 unique monophyletic clades..."
  [1] "Command: getTreeClades(trees,return_counts = TRUE)"
  [1] "All trees share 5 clades, including:"
  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_C;Species_H"                                                                                                              
  [3] "Species_G;Species_N"                                                                                                              
  [4] "Species_J;Species_M"                                                                                                              
  [5] "Species_K;Species_L"                                                                                                              
  [1] "Command: getTreeClades(trees,return_shared = TRUE)"
  [1] "Of 5 raw trees, there were 3 unique topologies..."
  [1] "Command: getUniqueTopologies(trees,print_table = TRUE)"
     Tree_Name              Trees_with_Topology Tree_Count Tree_Percent
  1 Topology_1 Gene_1.nwk;Gene_2.nwk;Gene_3.nwk          3           60
  2 Topology_2                       Gene_4.nwk          1           20
  3 Topology_3                       Gene_5.nwk          1           20
  
  # Simulate trees with different numbers of taxa
  mySimTrees <- c(rtree(10),rtree(10),rtree(20))
  summarizeMultiPhylo(mySimTrees)

  [1] "Trees unnamed, adding dummy names..."
  [1] "Command: treeNamer(trees)"
  [1] "Read in 3 trees, that all share 10 tip labels..."
  [1] "Trees do not have the same taxa, trimming to common tip labels..."
  [1] "Command: treeTrimmer(trees,getSharedTaxa(trees)"
  [1] "Among all trees, including root splits there are 24 unique monophyletic clades..."
  [1] "Command: getTreeClades(trees,return_counts = TRUE)"
  [1] "Trees in 'trees' share no clades in common..."
  [1] "All trees have a unique topology..."  
