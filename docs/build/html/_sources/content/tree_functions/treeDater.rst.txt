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

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    111. 
   2 Species_A;Species_F                                                                                                                                       18.5
   3 Species_A;Species_F;Species_K;Species_L                                                                                                                   37.0
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         55.5
   5 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       74.1
   6 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         92.6
   7 Species_B;Species_D;Species_O                                                                                                                             61.7
   8 Species_B;Species_O                                                                                                                                       30.9
   9 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         55.5
  10 Species_G;Species_I;Species_N                                                                                                                             37.0
  11 Species_G;Species_N                                                                                                                                       18.5
  12 Species_J;Species_M                                                                                                                                       27.8
  13 Species_K;Species_L                                                                                                                                       18.5
  14 Species_C;Species_H                                                                                                                                       55.5
  
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

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110. 
   2 Species_A;Species_F                                                                                                                                       16.2
   3 Species_A;Species_F;Species_K;Species_L                                                                                                                   34.9
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         53.6
   5 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       72.2
   6 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         90.9
   7 Species_B;Species_D;Species_O                                                                                                                             60.7
   8 Species_B;Species_O                                                                                                                                       30.4
   9 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         54.2
  10 Species_G;Species_I;Species_N                                                                                                                             36.2
  11 Species_G;Species_N                                                                                                                                       18.2
  12 Species_J;Species_M                                                                                                                                       27.2
  13 Species_K;Species_L                                                                                                                                       17.6
  14 Species_C;Species_H                                                                                                                                       54.9
  
  # Date a multiPhylo object where all trees share a common topology
  myTrees <- readRooted(c(rb_tree1_path,rb_tree2_path,rb_tree3_path),root_taxa = c('Species_C','Species_H'))

  myDatedTrees <- treeDater(tree = myTrees, taxa='Species_C;Species_M',min_max = c(100,120),iterations = 100)
  
  myDatedTrees
  3 phylogenetic trees

  all(is.ultrametric(myDatedTrees))
  [1] TRUE

  extractNodeAges(myDatedTrees,return_summary = 'both')

  # A tibble: 14 x 5
     Clade                                                                                                                                                 Mean_Node_Age Median_Node_Age StdDev_Node_Age MAD_Node_Age
     <chr>                                                                                                                                                         <dbl>           <dbl>           <dbl>        <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         110.            109.            0.847       0.0740
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                              91.3            90.9           0.706       0.0616
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                            73.1            72.8           0.565       0.0493
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                              54.8            54.6           0.423       0.0370
   5 Species_A;Species_F                                                                                                                                            18.3            18.2           0.141       0.0123
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                        36.5            36.4           0.282       0.0247
   7 Species_B;Species_D;Species_O                                                                                                                                  60.9            60.6           0.471       0.0411
   8 Species_B;Species_O                                                                                                                                            30.4            30.3           0.235       0.0205
   9 Species_C;Species_H                                                                                                                                            54.8            54.6           0.423       0.0370
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                              54.8            54.6           0.423       0.0370
  11 Species_G;Species_I;Species_N                                                                                                                                  36.5            36.4           0.282       0.0247
  12 Species_G;Species_N                                                                                                                                            18.3            18.2           0.141       0.0123
  13 Species_J;Species_M                                                                                                                                            27.4            27.3           0.212       0.0185
  14 Species_K;Species_L                                                                                                                                            18.3            18.2           0.141       0.0123
