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
  timeTree_2 <- readRooted(to_root = rb_timeTree2_path,root_taxa = c('Species_C','Species_H'))
  timeTree_3 <- readRooted(to_root = rb_timeTree3_path,root_taxa = c('Species_C','Species_H'))
  timeTree_4 <- readRooted(to_root = rb_timeTree4_path,root_taxa = c('Species_C','Species_H'))
  timeTree_5 <- readRooted(to_root = rb_timeTree5_path,root_taxa = c('Species_C','Species_H'))

  extractNodeAges(timeTree_1)
  
  # A tibble: 14 x 2
   Clade                                                                                                                                                 Node_Age
   <chr>                                                                                                                                                    <dbl>
  1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110. 
  2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         91.4
  3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       73.1
  4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         54.9
  5 Species_A;Species_F;Species_K;Species_L                                                                                                                   36.6
  6 Species_A;Species_F                                                                                                                                       18.3
  7 Species_K;Species_L                                                                                                                                       18.3
  8 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         54.9
  9 Species_G;Species_I;Species_N                                                                                                                             36.6
  10 Species_G;Species_N                                                                                                                                       18.3
  11 Species_J;Species_M                                                                                                                                       27.4
  12 Species_B;Species_D;Species_O                                                                                                                             60.9
  13 Species_B;Species_O                                                                                                                                       30.5
  14 Species_C;Species_H                                                                                                                                       54.9

  extractNodeAges(timeTree_2)

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110. 
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         91.3
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       73.1
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         54.8
   5 Species_A;Species_F;Species_K;Species_L                                                                                                                   36.5
   6 Species_A;Species_F                                                                                                                                       18.3
   7 Species_K;Species_L                                                                                                                                       18.3
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         54.8
   9 Species_G;Species_I;Species_N                                                                                                                             36.5
  10 Species_G;Species_N                                                                                                                                       18.3
  11 Species_J;Species_M                                                                                                                                       27.4
  12 Species_B;Species_D;Species_O                                                                                                                             60.9
  13 Species_B;Species_O                                                                                                                                       30.4
  14 Species_C;Species_H                                                                                                                                       54.8
  
  extractNodeAges(timeTree_3)
  
  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110. 
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         91.5
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       73.2
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         54.9
   5 Species_A;Species_F;Species_K;Species_L                                                                                                                   36.6
   6 Species_A;Species_F                                                                                                                                       18.3
   7 Species_K;Species_L                                                                                                                                       18.3
   8 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         54.9
   9 Species_G;Species_I;Species_N                                                                                                                             36.6
  10 Species_G;Species_N                                                                                                                                       18.3
  11 Species_J;Species_M                                                                                                                                       27.4
  12 Species_B;Species_D;Species_O                                                                                                                             61.0
  13 Species_B;Species_O                                                                                                                                       30.5
  14 Species_C;Species_H                                                                                                                                       54.9

  extractNodeAges(timeTree_4)

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110. 
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         88.0
   3 Species_A;Species_B;Species_D;Species_E;Species_F;Species_K;Species_L;Species_O                                                                           66.0
   4 Species_A;Species_E;Species_F                                                                                                                             44.0
   5 Species_A;Species_F                                                                                                                                       22.0
   6 Species_B;Species_D;Species_K;Species_L;Species_O                                                                                                         44.0
   7 Species_B;Species_O                                                                                                                                       22.0
   8 Species_D;Species_K;Species_L                                                                                                                             22.0
   9 Species_K;Species_L                                                                                                                                       11.0
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         66.0
  11 Species_G;Species_I;Species_N                                                                                                                             44.0
  12 Species_G;Species_N                                                                                                                                       22.0
  13 Species_J;Species_M                                                                                                                                       33.0
  14 Species_C;Species_H                                                                                                                                       55.0

  extractNodeAges(timeTree_5)
  
  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110. 
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         82.5
   3 Species_A;Species_F;Species_K;Species_L                                                                                                                   55.0
   4 Species_A;Species_F                                                                                                                                       27.5
   5 Species_K;Species_L                                                                                                                                       27.5
   6 Species_B;Species_D;Species_E;Species_G;Species_I;Species_J;Species_M;Species_N;Species_O                                                                 61.9
   7 Species_B;Species_D;Species_O                                                                                                                             41.3
   8 Species_B;Species_D                                                                                                                                       20.6
   9 Species_E;Species_G;Species_I;Species_J;Species_M;Species_N                                                                                               41.3
  10 Species_E;Species_I                                                                                                                                       20.6
  11 Species_G;Species_J;Species_M;Species_N                                                                                                                   27.5
  12 Species_G;Species_N                                                                                                                                       13.8
  13 Species_J;Species_M                                                                                                                                       13.8
  14 Species_C;Species_H                                                                                                                                       55.0

  # Create multiPhylo of trees that share a topology
  timeTrees <- c(timeTree_1,timeTree_2,timeTree_3) %>% treeNamer()

  extractNodeAges(timeTrees) %>% arrange(Clade)

  # A tibble: 42 x 3
     Clade                                                                                                                                                 Node_Age Tree_Name
     <chr>                                                                                                                                                    <dbl> <chr>    
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110.  Tree_1   
   2 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110.  Tree_2   
   3 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    110.  Tree_3   
   4 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         91.4 Tree_1   
   5 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         91.3 Tree_2   
   6 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                         91.5 Tree_3   
   7 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       73.1 Tree_1   
   8 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       73.1 Tree_2   
   9 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                       73.2 Tree_3   
  10 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         54.9 Tree_1   
  # ... with 32 more rows

  extractNodeAges(timeTrees,return_summary = 'mean')
  
  # A tibble: 14 x 3
     Clade                                                                                                                                                 Mean_Node_Age StdDev_Node_Age
     <chr>                                                                                                                                                         <dbl>           <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         110.           0.0811
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                              91.4          0.0676
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                            73.1          0.0541
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                              54.9          0.0406
   5 Species_A;Species_F                                                                                                                                            18.3          0.0135
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                        36.6          0.0270
   7 Species_B;Species_D;Species_O                                                                                                                                  60.9          0.0451
   8 Species_B;Species_O                                                                                                                                            30.5          0.0225
   9 Species_C;Species_H                                                                                                                                            54.9          0.0406
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                              54.9          0.0406
  11 Species_G;Species_I;Species_N                                                                                                                                  36.6          0.0270
  12 Species_G;Species_N                                                                                                                                            18.3          0.0135
  13 Species_J;Species_M                                                                                                                                            27.4          0.0203
  14 Species_K;Species_L                                                                                                                                            18.3          0.0135

  extractNodeAges(timeTrees,return_summary = 'median')

  # A tibble: 14 x 3
     Clade                                                                                                                                                 Median_Node_Age MAD_Node_Age
     <chr>                                                                                                                                                           <dbl>        <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O           110.        0.107 
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                                91.4       0.0890
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                              73.1       0.0712
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                                54.9       0.0534
   5 Species_A;Species_F                                                                                                                                              18.3       0.0178
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                          36.6       0.0356
   7 Species_B;Species_D;Species_O                                                                                                                                    60.9       0.0594
   8 Species_B;Species_O                                                                                                                                              30.5       0.0297
   9 Species_C;Species_H                                                                                                                                              54.9       0.0534
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                                54.9       0.0534
  11 Species_G;Species_I;Species_N                                                                                                                                    36.6       0.0356
  12 Species_G;Species_N                                                                                                                                              18.3       0.0178
  13 Species_J;Species_M                                                                                                                                              27.4       0.0267
  14 Species_K;Species_L                                                                                                                                              18.3       0.0178

  extractNodeAges(timeTrees,return_summary = 'both')
  
  # A tibble: 14 x 5
     Clade                                                                                                                                                 Mean_Node_Age Median_Node_Age StdDev_Node_Age MAD_Node_Age
     <chr>                                                                                                                                                         <dbl>           <dbl>           <dbl>        <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         110.            110.           0.0811       0.107 
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                              91.4            91.4          0.0676       0.0890
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                            73.1            73.1          0.0541       0.0712
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                              54.9            54.9          0.0406       0.0534
   5 Species_A;Species_F                                                                                                                                            18.3            18.3          0.0135       0.0178
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                        36.6            36.6          0.0270       0.0356
   7 Species_B;Species_D;Species_O                                                                                                                                  60.9            60.9          0.0451       0.0594
   8 Species_B;Species_O                                                                                                                                            30.5            30.5          0.0225       0.0297
   9 Species_C;Species_H                                                                                                                                            54.9            54.9          0.0406       0.0534
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                              54.9            54.9          0.0406       0.0534
  11 Species_G;Species_I;Species_N                                                                                                                                  36.6            36.6          0.0270       0.0356
  12 Species_G;Species_N                                                                                                                                            18.3            18.3          0.0135       0.0178
  13 Species_J;Species_M                                                                                                                                            27.4            27.4          0.0203       0.0267
  14 Species_K;Species_L                                                                                                                                            18.3            18.3          0.0135       0.0178
