.. _getTreeSplits:

##################
**getTreeSplits**
##################

*getTreeSplits* breaks rooted phylo and multiPhylo objects down into their respective splits. 

  - If a single topology is supplied (via 1 or more trees), a single dataframe is returned.
  - If multiple unqique toplogies are present, a list of dataframes is returned.

=======================
Function and Arguments
=======================

**Usage**:
::

  getTreeSplits <- function(tree)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**tree**				                  A rooted phylo object, or a rooted multiPhylo object where all trees share 3+ taxa
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getTreeSplits.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()

  # Read in trees with two topologies (1=2, 1!=3)

  myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
  myTree_2 <- readRooted(rb_tree2_path,root_taxa = c('Species_C','Species_H'))
  myTree_3 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

  mySameTrees <- c(myTree_1,myTree_2)
  myDifferentTrees <- c(myTree_1,myTree_3)

  # Topology A

  getTreeSplits(myTree_1)

  # A tibble: 12 x 4
   Clade                                                                                           Mirror_Clade                                                                                                                 Split_Node Root 
   <chr>                                                                                           <chr>                                                                                                                             <int> <lgl>
  1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Spec~ Species_B;Species_C;Species_D;Species_H;Species_O                                                                                    18 FALSE
  2 Species_A;Species_E;Species_F;Species_K;Species_L                                               Species_B;Species_C;Species_D;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                                  19 FALSE
  3 Species_A;Species_F;Species_K;Species_L                                                         Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                        20 FALSE
  4 Species_K;Species_L                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Spe~         21 FALSE
  5 Species_A;Species_F                                                                             Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         22 FALSE
  6 Species_G;Species_I;Species_J;Species_M;Species_N                                               Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_K;Species_L;Species_O                                  23 FALSE
  7 Species_J;Species_M                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_K;Species_L;Species_N;Spe~         24 FALSE
  8 Species_G;Species_I;Species_N                                                                   Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_J;Species_K;Species_L;Species_M;Species_O              25 FALSE
  9 Species_G;Species_N                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         26 FALSE
  10 Species_B;Species_D;Species_O                                                                   Species_A;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N              27 FALSE
  11 Species_B;Species_O                                                                             Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         28 FALSE
  12 Species_C;Species_H                                                                             Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         16 TRUE 

  getTreeSplits(myTree_2)

  # A tibble: 12 x 4
   Clade                                                                                           Mirror_Clade                                                                                                                 Split_Node Root 
   <chr>                                                                                           <chr>                                                                                                                             <int> <lgl>
  1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Spec~ Species_B;Species_C;Species_D;Species_H;Species_O                                                                                    18 FALSE
  2 Species_A;Species_E;Species_F;Species_K;Species_L                                               Species_B;Species_C;Species_D;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                                  19 FALSE
  3 Species_A;Species_F;Species_K;Species_L                                                         Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                        20 FALSE
  4 Species_K;Species_L                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Spe~         21 FALSE
  5 Species_A;Species_F                                                                             Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         22 FALSE
  6 Species_G;Species_I;Species_J;Species_M;Species_N                                               Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_K;Species_L;Species_O                                  23 FALSE
  7 Species_J;Species_M                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_K;Species_L;Species_N;Spe~         24 FALSE
  8 Species_G;Species_I;Species_N                                                                   Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_J;Species_K;Species_L;Species_M;Species_O              25 FALSE
  9 Species_G;Species_N                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         26 FALSE
  10 Species_B;Species_D;Species_O                                                                   Species_A;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N              27 FALSE
  11 Species_B;Species_O                                                                             Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         28 FALSE
  12 Species_C;Species_H                                                                             Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         16 TRUE 

  # Topology B

  getTreeSplits(myTree_3)

  # A tibble: 12 x 4
   Clade                                                                           Mirror_Clade                                                                                                                      Split_Node Root 
   <chr>                                                                           <chr>                                                                                                                                  <int> <lgl>
  1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_M;Species_N Species_B;Species_C;Species_D;Species_H;Species_K;Species_L;Species_O                                                                     18 FALSE
  2 Species_A;Species_E;Species_F                                                   Species_B;Species_C;Species_D;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                   19 FALSE
  3 Species_A;Species_F                                                             Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         20 FALSE
  4 Species_G;Species_I;Species_J;Species_M;Species_N                               Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_K;Species_L;Species_O                                       21 FALSE
  5 Species_J;Species_M                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_K;Species_L;Species_N;Species_O         22 FALSE
  6 Species_G;Species_I;Species_N                                                   Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_J;Species_K;Species_L;Species_M;Species_O                   23 FALSE
  7 Species_G;Species_N                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_O         24 FALSE
  8 Species_B;Species_D;Species_K;Species_L;Species_O                               Species_A;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N                                       25 FALSE
  9 Species_B;Species_O                                                             Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N         26 FALSE
  10 Species_D;Species_K;Species_L                                                   Species_A;Species_B;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                   27 FALSE
  11 Species_K;Species_L                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O         28 FALSE
  12 Species_C;Species_H                                                             Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         16 TRUE 

  # Get splits when trees are identical

  getTreeSplits(mySameTrees)

  [1] "All trees supplied to getTreeSplits share a common topology...returning results from common topology..."

  # A tibble: 12 x 4
   Clade                                                                                           Mirror_Clade                                                                                                                 Split_Node Root 
   <chr>                                                                                           <chr>                                                                                                                             <int> <lgl>
  1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Spec~ Species_B;Species_C;Species_D;Species_H;Species_O                                                                                    18 FALSE
  2 Species_A;Species_E;Species_F;Species_K;Species_L                                               Species_B;Species_C;Species_D;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                                  19 FALSE
  3 Species_A;Species_F;Species_K;Species_L                                                         Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                        20 FALSE
  4 Species_K;Species_L                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Spe~         21 FALSE
  5 Species_A;Species_F                                                                             Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         22 FALSE
  6 Species_G;Species_I;Species_J;Species_M;Species_N                                               Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_K;Species_L;Species_O                                  23 FALSE
  7 Species_J;Species_M                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_K;Species_L;Species_N;Spe~         24 FALSE
  8 Species_G;Species_I;Species_N                                                                   Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_J;Species_K;Species_L;Species_M;Species_O              25 FALSE
  9 Species_G;Species_N                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         26 FALSE
  10 Species_B;Species_D;Species_O                                                                   Species_A;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N              27 FALSE
  11 Species_B;Species_O                                                                             Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         28 FALSE
  12 Species_C;Species_H                                                                             Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         16 TRUE 

  # Get splits when trees are not identical

  getTreeSplits(myDifferentTrees)

  $Tree_1

  # A tibble: 12 x 4
   Clade                                                                                           Mirror_Clade                                                                                                                 Split_Node Root 
   <chr>                                                                                           <chr>                                                                                                                             <int> <lgl>
  1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Spec~ Species_B;Species_C;Species_D;Species_H;Species_O                                                                                    18 FALSE
  2 Species_A;Species_E;Species_F;Species_K;Species_L                                               Species_B;Species_C;Species_D;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                                  19 FALSE
  3 Species_A;Species_F;Species_K;Species_L                                                         Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                        20 FALSE
  4 Species_K;Species_L                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Spe~         21 FALSE
  5 Species_A;Species_F                                                                             Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         22 FALSE
  6 Species_G;Species_I;Species_J;Species_M;Species_N                                               Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_K;Species_L;Species_O                                  23 FALSE
  7 Species_J;Species_M                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_K;Species_L;Species_N;Spe~         24 FALSE
  8 Species_G;Species_I;Species_N                                                                   Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_J;Species_K;Species_L;Species_M;Species_O              25 FALSE
  9 Species_G;Species_N                                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         26 FALSE
  10 Species_B;Species_D;Species_O                                                                   Species_A;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N              27 FALSE
  11 Species_B;Species_O                                                                             Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Spe~         28 FALSE
  12 Species_C;Species_H                                                                             Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Spe~         16 TRUE 

  $Tree_2

  # A tibble: 12 x 4
   Clade                                                                           Mirror_Clade                                                                                                                      Split_Node Root 
   <chr>                                                                           <chr>                                                                                                                                  <int> <lgl>
  1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_M;Species_N Species_B;Species_C;Species_D;Species_H;Species_K;Species_L;Species_O                                                                     18 FALSE
  2 Species_A;Species_E;Species_F                                                   Species_B;Species_C;Species_D;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                   19 FALSE
  3 Species_A;Species_F                                                             Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         20 FALSE
  4 Species_G;Species_I;Species_J;Species_M;Species_N                               Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_K;Species_L;Species_O                                       21 FALSE
  5 Species_J;Species_M                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_K;Species_L;Species_N;Species_O         22 FALSE
  6 Species_G;Species_I;Species_N                                                   Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_J;Species_K;Species_L;Species_M;Species_O                   23 FALSE
  7 Species_G;Species_N                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_O         24 FALSE
  8 Species_B;Species_D;Species_K;Species_L;Species_O                               Species_A;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N                                       25 FALSE
  9 Species_B;Species_O                                                             Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N         26 FALSE
  10 Species_D;Species_K;Species_L                                                   Species_A;Species_B;Species_C;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O                   27 FALSE
  11 Species_K;Species_L                                                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_M;Species_N;Species_O         28 FALSE
  12 Species_C;Species_H                                                             Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         16 TRUE   
  
