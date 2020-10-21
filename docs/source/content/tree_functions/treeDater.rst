.. _treeDater:

###############
**treeDater**
###############

**treeDater** takes phylo/multiPhylo objects where branch lengths are substitution rates, and calls *chronos* from the **ape** package to create an ultrametric, node-calibrated dated phylo/multiPhylo.

=======================
Function and Arguments
=======================

**Usage**:
::

  treeDater <- function(tree,calibration_df,taxa,min_max,iterations)

===========================      ===================================================================================================================================================================================================================================
 Argument                         Description
===========================      ===================================================================================================================================================================================================================================
**tree**				                  A rooted phylo or multiPhylo object (all trees sharing the same topology) where branch lengths represent substitution rates
**calibration_df**                A 4-column dataframe/tibble with calibration information: (1) Taxon 1 (2) Taxon 2 [to get MRCA] (3) Min divergence time bound (4) Max divergence time bound; Multiple calibration points are allowed [Supercedes 'taxa' argument]
**taxa**                          A character vector (or semicolon-delimited set) of taxon IDs from which to find the MRCA and calibrate [One calibration point allowed; Superceded by 'calibration_df']
**min_max**                       If using 'taxa', a two-element numeric vector [e.g. c(50,75)] that provides the minimum and maximum age estimates for the focal calibration node [Superceded by 'calibration_df'; min <= max]
**iterations**                    How many times to estimate the age of each node prior to summarizing [Default: 1000]
===========================      ===================================================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/treeDater.R
  
  library(Rboretum)
  sourceRboretum()

  # Date a phylo object where branch lengths are substitution rates
  myTree <- readRooted(rb_tree1_path,c('Species_C','Species_H'))
  is.ultrametric(myTree)
  [1] FALSE

  # Estimate node ages by calibrating the root node to between 100MY and 120MY, iterating estimates 100 times
  myDatedTree <- treeDater(tree = myTree, taxa='Species_C;Species_M',min_max = c(100,120),iterations = 100)
  is.ultrametric(myDatedTree)
  [1] TRUE

  extractNodeAges(myDatedTree)
  
  # A tibble: 12 x 2
     Clade                                                                                               Node_Age
     <chr>                                                                                                  <dbl>
   1 Species_A;Species_F                                                                                     18.4
   2 Species_A;Species_F;Species_K;Species_L                                                                 36.8
   3 Species_A;Species_E;Species_F;Species_K;Species_L                                                       55.3
   4 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N     73.7
   5 Species_B;Species_D;Species_O                                                                           61.4
   6 Species_B;Species_O                                                                                     30.7
   7 Species_G;Species_I;Species_J;Species_M;Species_N                                                       55.3
   8 Species_G;Species_I;Species_N                                                                           36.8
   9 Species_G;Species_N                                                                                     18.4
  10 Species_J;Species_M                                                                                     27.6
  11 Species_K;Species_L                                                                                     18.4
  12 Species_C;Species_H                                                                                    111. 
  
  # Estimate node ages by calibrating at two nodes
  myCalibration <- tibble(Taxon_A = c('Species_C','Species_A'),
                          Taxon_B = c('Species_M','Species_F'),
                          Min = c(100,15),
                          Max = c(120,17))
  myCalibration

  # A tibble: 2 x 4
    Taxon_A   Taxon_B     Min   Max
    <chr>     <chr>     <dbl> <dbl>
  1 Species_C Species_M   100   120
  2 Species_A Species_F    15    17
  
  myRedatedTree <- treeDater(tree = myTree, calibration_df = myCalibration,iterations = 100)

  is.ultrametric(myRedatedTree)
  [1] TRUE

  extractNodeAges(myRedatedTree)
  # A tibble: 12 x 2
     Clade                                                                                               Node_Age
     <chr>                                                                                                  <dbl>
   1 Species_A;Species_F                                                                                     16.1
   2 Species_A;Species_F;Species_K;Species_L                                                                 34.9
   3 Species_A;Species_E;Species_F;Species_K;Species_L                                                       53.7
   4 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N     72.5
   5 Species_B;Species_D;Species_O                                                                           60.9
   6 Species_B;Species_O                                                                                     30.5
   7 Species_G;Species_I;Species_J;Species_M;Species_N                                                       54.4
   8 Species_G;Species_I;Species_N                                                                           36.3
   9 Species_G;Species_N                                                                                     18.1
  10 Species_J;Species_M                                                                                     27.2
  11 Species_K;Species_L                                                                                     17.5
  12 Species_C;Species_H                                                                                    110. 
  
  # Date a multiPhylo object where all trees share a common topology
  myTrees <- readRooted(c(rb_tree1_path,rb_tree2_path,rb_tree3_path),root_taxa = c('Species_C','Species_H'))

  myDatedTrees <- treeDater(tree = myTrees, taxa='Species_C;Species_M',min_max = c(100,120),iterations = 100)
  myDatedTrees
  3 phylogenetic trees

  all(is.ultrametric(myDatedTrees))
  [1] TRUE

  extractNodeAges(myDatedTrees,return_summary = 'both')
  # A tibble: 12 x 5
     Clade                                                                                               Mean_Node_Age Median_Node_Age StdDev_Node_Age MAD_Node_Age
     <chr>                                                                                                       <dbl>           <dbl>           <dbl>        <dbl>
   1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N          73.2            73.0           1.46         1.62 
   2 Species_A;Species_E;Species_F;Species_K;Species_L                                                            54.9            54.7           1.09         1.22 
   3 Species_A;Species_F                                                                                          18.3            18.2           0.364        0.405
   4 Species_A;Species_F;Species_K;Species_L                                                                      36.6            36.5           0.728        0.811
   5 Species_B;Species_D;Species_O                                                                                61.0            60.8           1.21         1.35 
   6 Species_B;Species_O                                                                                          30.5            30.4           0.606        0.676
   7 Species_C;Species_H                                                                                         110.            109.            2.18         2.43 
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                            54.9            54.7           1.09         1.22 
   9 Species_G;Species_I;Species_N                                                                                36.6            36.5           0.728        0.811
  10 Species_G;Species_N                                                                                          18.3            18.2           0.364        0.405
  11 Species_J;Species_M                                                                                          27.5            27.4           0.546        0.608
  12 Species_K;Species_L                                                                                          18.3            18.2           0.364        0.405
