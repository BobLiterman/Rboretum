.. _treeDater:

###############
**treeDater**
###############

*treeDater* takes phylo/multiPhylo objects where branch lengths are substitution rates, and calls *chronos* from the **ape** package to create an ultrametric, node-calibrated dated phylo/multiPhylo.

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

  # Set test data directory
  sourceRboretum()
  
  # Date a phylo object where branch lengths are substitution rates
  myTree <- readRooted(rb_tree1_path,c('Species_C','Species_H'))

  # Estimate node ages by calibrating the root node to between 100MY and 120MY, iterating estimates 100 times
  myDatedTree <- treeDater(tree = myTree, taxa='Species_C;Species_M',min_max = c(100,120),iterations = 100)

  is.ultrametric(myDatedTree)
  TRUE

  extractNodeAges(myDatedTree)

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O   121.  
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                        64.2 
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                      56.6 
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                        30.5 
   5 Species_A;Species_F;Species_K;Species_L                                                                                                                  24.4 
   6 Species_K;Species_L                                                                                                                                      11.1 
   7 Species_A;Species_F                                                                                                                                       9.90
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                        38.3 
   9 Species_J;Species_M                                                                                                                                      13.3 
  10 Species_G;Species_I;Species_N                                                                                                                            25.7 
  11 Species_G;Species_N                                                                                                                                      17.8 
  12 Species_B;Species_D;Species_O                                                                                                                            23.2 
  13 Species_B;Species_O                                                                                                                                      18.0 
  14 Species_C;Species_H                                                                                                                                      55.7 

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
  TRUE

  extractNodeAges(myRedatedTree)

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    121. 
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         69.5
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       62.1
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         36.2
   5 Species_A;Species_F;Species_K;Species_L                                                                                                                   30.6
   6 Species_K;Species_L                                                                                                                                       11.3
   7 Species_A;Species_F                                                                                                                                       16.8
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         42.2
   9 Species_J;Species_M                                                                                                                                       13.5
  10 Species_G;Species_I;Species_N                                                                                                                             26.8
  11 Species_G;Species_N                                                                                                                                       18.4
  12 Species_B;Species_D;Species_O                                                                                                                             18.7
  13 Species_B;Species_O                                                                                                                                       12.5
  14 Species_C;Species_H                                                                                                                                       56.0
   

  # Date a multiPhylo object
  myTrees <- readRooted(c(rb_tree1_path,rb_tree2_path,rb_timeTree3_path),c('Species_C','Species_H'))

  myDatedTrees <- treeDater(tree = myTrees, taxa='Species_C;Species_M',min_max = c(100,120),iterations = 100)

  all(is.ultrametric(myDatedTrees))
  TRUE

  extractNodeAges(myDatedTrees,return_summary = 'both')

  # A tibble: 14 x 5
     Clade                                                                                                                                                 Mean_Node_Age Median_Node_Age StdDev_Node_Age MAD_Node_Age
     <chr>                                                                                                                                                         <dbl>           <dbl>           <dbl>        <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         121.           120.             1.01         0.304
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                              64.7           60.6           21.8         22.7  
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                            56.3           52.6           19.5         20.4  
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                              31.4           26.9           17.3         15.0  
   5 Species_A;Species_F                                                                                                                                            10.3            8.36           7.00         5.76 
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                        24.4           21.4           14.0         13.8  
   7 Species_B;Species_D;Species_O                                                                                                                                  25.0           19.6           20.1         17.1  
   8 Species_B;Species_O                                                                                                                                            16.3           14.7           13.9         16.7  
   9 Species_C;Species_H                                                                                                                                            56.1           56.0            0.553        0.420
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                              33.8           34.4           15.0         20.8  
  11 Species_G;Species_I;Species_N                                                                                                                                  19.0           22.6            9.20         4.71 
  12 Species_G;Species_N                                                                                                                                            11.1           13.0            6.26         4.85 
  13 Species_J;Species_M                                                                                                                                            13.9           11.5            9.74         8.65 
  14 Species_K;Species_L                                                                                                                                            10.7            9.72           7.38         8.64 
