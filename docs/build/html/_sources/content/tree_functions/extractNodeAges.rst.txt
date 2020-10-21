.. _extractNodeAges:

####################
**extractNodeAges**
####################
.. |br| raw:: html

   <br />
   
**extractNodeAges** takes ultrametric phylo or multiPhylo objects and returns information about branching times.
|br|
|br|
If a multiPhylo is provided (where all trees share a single topology), branching times will be extracted from each tree and can be:

  - Returned individually, or
  - Summarized via *mean* or *median* (or *both*) and returned in a table.

=======================
Function and Arguments
=======================

**Usage**:
::

  extractNodeAges <- function(tree,return_summary)

===========================      =========================================================================================================================================================================================================================
 Argument                         Description
===========================      =========================================================================================================================================================================================================================
**tree**				                  A rooted, ultrametric phylo object, or a rooted, ultrametric multiPhylo object where all trees share 3+ taxa and a common topology (after pruning if necessary)
**return_summary**                **mean**, **median**, or **both**; If  and a multiPhylo is provided, return dates summarized across trees using the mean, median, or a full summary ('both') [Default: Return raw results from each tree, unsummarized] 
===========================      =========================================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/extractNodeAges.R

  library(Rboretum)
  sourceRboretum()
  
  # Read in ultrametric timetrees (branch lengths ~ time)
  timeTree_1 <- readRooted(to_root = rb_timeTree1_path,root_taxa = c('Species_C','Species_H'))
  
  extractNodeAges(timeTree_1)

  # A tibble: 12 x 2
     Clade                                                                                               Node_Age
     <chr>                                                                                                  <dbl>
   1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N     73.1
   2 Species_A;Species_E;Species_F;Species_K;Species_L                                                       54.9
   3 Species_A;Species_F;Species_K;Species_L                                                                 36.6
   4 Species_A;Species_F                                                                                     18.3
   5 Species_K;Species_L                                                                                     18.3
   6 Species_G;Species_I;Species_J;Species_M;Species_N                                                       54.9
   7 Species_G;Species_I;Species_N                                                                           36.6
   8 Species_G;Species_N                                                                                     18.3
   9 Species_J;Species_M                                                                                     27.4
  10 Species_B;Species_D;Species_O                                                                           60.9
  11 Species_B;Species_O                                                                                     30.5
  12 Species_C;Species_H                                                                                    110. 
  
  # Create multiPhylo of trees that share a topology
  timeTrees <- c(timeTree_1,timeTree_2,timeTree_3) %>% treeNamer()
  
  extractNodeAges(timeTrees) %>% arrange(Clade)

  # A tibble: 36 x 3
     Clade                                                                                               Node_Age Tree_Name
     <chr>                                                                                                  <dbl> <chr>    
   1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N     73.1 Tree_1   
   2 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N     73.1 Tree_2   
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N     73.2 Tree_3   
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                       54.9 Tree_1   
   5 Species_A;Species_E;Species_F;Species_K;Species_L                                                       54.8 Tree_2   
   6 Species_A;Species_E;Species_F;Species_K;Species_L                                                       54.9 Tree_3   
   7 Species_A;Species_F                                                                                     18.3 Tree_1   
   8 Species_A;Species_F                                                                                     18.3 Tree_2   
   9 Species_A;Species_F                                                                                     18.3 Tree_3   
  10 Species_A;Species_F;Species_K;Species_L                                                                 36.6 Tree_1   
  # ... with 26 more rows
  
  extractNodeAges(timeTrees,return_summary = 'mean')

  # A tibble: 12 x 3
     Clade                                                                                               Mean_Node_Age StdDev_Node_Age
     <chr>                                                                                                       <dbl>           <dbl>
   1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N          73.1          0.0541
   2 Species_A;Species_E;Species_F;Species_K;Species_L                                                            54.9          0.0406
   3 Species_A;Species_F                                                                                          18.3          0.0135
   4 Species_A;Species_F;Species_K;Species_L                                                                      36.6          0.0270
   5 Species_B;Species_D;Species_O                                                                                60.9          0.0451
   6 Species_B;Species_O                                                                                          30.5          0.0225
   7 Species_C;Species_H                                                                                         110.           0.0811
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                            54.9          0.0406
   9 Species_G;Species_I;Species_N                                                                                36.6          0.0270
  10 Species_G;Species_N                                                                                          18.3          0.0135
  11 Species_J;Species_M                                                                                          27.4          0.0203
  12 Species_K;Species_L                                                                                          18.3          0.0135
  
  extractNodeAges(timeTrees,return_summary = 'median')

  # A tibble: 12 x 3
     Clade                                                                                               Median_Node_Age MAD_Node_Age
     <chr>                                                                                                         <dbl>        <dbl>
   1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N            73.1       0.0712
   2 Species_A;Species_E;Species_F;Species_K;Species_L                                                              54.9       0.0534
   3 Species_A;Species_F                                                                                            18.3       0.0178
   4 Species_A;Species_F;Species_K;Species_L                                                                        36.6       0.0356
   5 Species_B;Species_D;Species_O                                                                                  60.9       0.0594
   6 Species_B;Species_O                                                                                            30.5       0.0297
   7 Species_C;Species_H                                                                                           110.        0.107 
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                              54.9       0.0534
   9 Species_G;Species_I;Species_N                                                                                  36.6       0.0356
  10 Species_G;Species_N                                                                                            18.3       0.0178
  11 Species_J;Species_M                                                                                            27.4       0.0267
  12 Species_K;Species_L                                                                                            18.3       0.0178
  
  extractNodeAges(timeTrees,return_summary = 'both')

  # A tibble: 12 x 5
     Clade                                                                                               Mean_Node_Age Median_Node_Age StdDev_Node_Age MAD_Node_Age
     <chr>                                                                                                       <dbl>           <dbl>           <dbl>        <dbl>
   1 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N          73.1            73.1          0.0541       0.0712
   2 Species_A;Species_E;Species_F;Species_K;Species_L                                                            54.9            54.9          0.0406       0.0534
   3 Species_A;Species_F                                                                                          18.3            18.3          0.0135       0.0178
   4 Species_A;Species_F;Species_K;Species_L                                                                      36.6            36.6          0.0270       0.0356
   5 Species_B;Species_D;Species_O                                                                                60.9            60.9          0.0451       0.0594
   6 Species_B;Species_O                                                                                          30.5            30.5          0.0225       0.0297
   7 Species_C;Species_H                                                                                         110.            110.           0.0811       0.107 
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                            54.9            54.9          0.0406       0.0534
   9 Species_G;Species_I;Species_N                                                                                36.6            36.6          0.0270       0.0356
  10 Species_G;Species_N                                                                                          18.3            18.3          0.0135       0.0178
  11 Species_J;Species_M                                                                                          27.4            27.4          0.0203       0.0267
  12 Species_K;Species_L                                                                                          18.3            18.3          0.0135       0.0178
