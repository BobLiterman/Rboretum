.. _extractNodeAges:

####################
**extractNodeAges**
####################

*extractNodeAges* takes ultrametric phylo or multiPhylo objects and returns information about branching times.

=======================
Function and Arguments
=======================

**Usage**:
::

  extractNodeAges <- function(tree,return_summary)

===========================      =====================================================================================================================================================================================================================
 Argument                         Description
===========================      =====================================================================================================================================================================================================================
**tree**				                  A rooted, ultrametric phylo object, or a rooted, ultrametric multiPhylo object where all trees share 3+ taxa and a common topology (after pruning if necessary)
**return_summary**                **mean**, **median**, or **both**; If  and a multiPhylo is provided, return dates summarized across trees using the mean, median, or a full summary ('both') [Default: Return raw results from each tree as a list] 
===========================      =====================================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/extractNodeAges.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()
  
  # Read in ultrametric timetrees (branch lengths ~ time)
  rb_timeTree1_path
  [1] "C:/Users/Robert.Literman/Documents/R/win-library/4.0/Rboretum/extdata/timeTrees/Chronogram_1.nwk"
  rb_timeTree2_path
  [1] "C:/Users/Robert.Literman/Documents/R/win-library/4.0/Rboretum/extdata/timeTrees/Chronogram_2.nwk"
  rb_timeTree3_path
  [1] "C:/Users/Robert.Literman/Documents/R/win-library/4.0/Rboretum/extdata/timeTrees/Chronogram_3.nwk"

  timeTree_1 <- readRooted(to_root = rb_timeTree1_path,root_taxa = c('Species_C','Species_H'))
  timeTree_2 <- readRooted(to_root = rb_timeTree2_path,root_taxa = c('Species_C','Species_H'))
  timeTree_3 <- readRooted(to_root = rb_timeTree3_path,root_taxa = c('Species_C','Species_H'))

  # Create multiPhylo of trees
  timeTrees <- c(timeTree_1,timeTree_2,timeTree_3)

  extractNodeAges(timeTree_1)

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    60.0 
   2 Species_C;Species_H                                                                                                                                      27.8 
   3 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                        24.4 
   4 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                      20.1 
   5 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         7.97
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                   5.95
   7 Species_K;Species_L                                                                                                                                       2.09
   8 Species_A;Species_F                                                                                                                                       1.71
   9 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                        11.0 
  10 Species_J;Species_M                                                                                                                                       2.91
  11 Species_G;Species_I;Species_N                                                                                                                             5.68
  12 Species_G;Species_N                                                                                                                                       3.42
  13 Species_B;Species_D;Species_O                                                                                                                             4.23
  14 Species_B;Species_O                                                                                                                                       2.24

  extractNodeAges(timeTree_2)

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    60.0 
   2 Species_C;Species_H                                                                                                                                      27.8 
   3 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                        26.8 
   4 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                      19.7 
   5 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         9.38
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                   6.83
   7 Species_K;Species_L                                                                                                                                       2.61
   8 Species_A;Species_F                                                                                                                                       2.24
   9 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         8.38
  10 Species_J;Species_M                                                                                                                                       2.54
  11 Species_G;Species_I;Species_N                                                                                                                             4.28
  12 Species_G;Species_N                                                                                                                                       2.01
  13 Species_B;Species_D;Species_O                                                                                                                             4.42
  14 Species_B;Species_O                                                                                                                                       2.53

  extractNodeAges(timeTree_3)

  # A tibble: 14 x 2
     Clade                                                                                                                                                 Node_Age
     <chr>                                                                                                                                                    <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    60.0 
   2 Species_C;Species_H                                                                                                                                      27.8 
   3 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                        22.6 
   4 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                      19.4 
   5 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         8.41
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                   6.05
   7 Species_K;Species_L                                                                                                                                       1.95
   8 Species_A;Species_F                                                                                                                                       2.24
   9 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         9.24
  10 Species_J;Species_M                                                                                                                                       2.82
  11 Species_G;Species_I;Species_N                                                                                                                             4.27
  12 Species_G;Species_N                                                                                                                                       2.07
  13 Species_B;Species_D;Species_O                                                                                                                             4.05
  14 Species_B;Species_O                                                                                                                                       1.69

  print(extractNodeAges(timeTrees),n = 42)

  # A tibble: 42 x 3
     Clade                                                                                                                                                 Node_Age Tree_Name
     <chr>                                                                                                                                                    <dbl> <chr>    
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    60.0  Tree_1   
   2 Species_C;Species_H                                                                                                                                      27.8  Tree_1   
   3 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                        24.4  Tree_1   
   4 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                      20.1  Tree_1   
   5 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         7.97 Tree_1   
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                   5.95 Tree_1   
   7 Species_K;Species_L                                                                                                                                       2.09 Tree_1   
   8 Species_A;Species_F                                                                                                                                       1.71 Tree_1   
   9 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                        11.0  Tree_1   
  10 Species_J;Species_M                                                                                                                                       2.91 Tree_1   
  11 Species_G;Species_I;Species_N                                                                                                                             5.68 Tree_1   
  12 Species_G;Species_N                                                                                                                                       3.42 Tree_1   
  13 Species_B;Species_D;Species_O                                                                                                                             4.23 Tree_1   
  14 Species_B;Species_O                                                                                                                                       2.24 Tree_1   
  15 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    60.0  Tree_2   
  16 Species_C;Species_H                                                                                                                                      27.8  Tree_2   
  17 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                        26.8  Tree_2   
  18 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                      19.7  Tree_2   
  19 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         9.38 Tree_2   
  20 Species_A;Species_F;Species_K;Species_L                                                                                                                   6.83 Tree_2   
  21 Species_K;Species_L                                                                                                                                       2.61 Tree_2   
  22 Species_A;Species_F                                                                                                                                       2.24 Tree_2   
  23 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         8.38 Tree_2   
  24 Species_J;Species_M                                                                                                                                       2.54 Tree_2   
  25 Species_G;Species_I;Species_N                                                                                                                             4.28 Tree_2   
  26 Species_G;Species_N                                                                                                                                       2.01 Tree_2   
  27 Species_B;Species_D;Species_O                                                                                                                             4.42 Tree_2   
  28 Species_B;Species_O                                                                                                                                       2.53 Tree_2   
  29 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O    60.0  Tree_3   
  30 Species_C;Species_H                                                                                                                                      27.8  Tree_3   
  31 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                        22.6  Tree_3   
  32 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                      19.4  Tree_3   
  33 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                         8.41 Tree_3   
  34 Species_A;Species_F;Species_K;Species_L                                                                                                                   6.05 Tree_3   
  35 Species_K;Species_L                                                                                                                                       1.95 Tree_3   
  36 Species_A;Species_F                                                                                                                                       2.24 Tree_3   
  37 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                         9.24 Tree_3   
  38 Species_J;Species_M                                                                                                                                       2.82 Tree_3   
  39 Species_G;Species_I;Species_N                                                                                                                             4.27 Tree_3   
  40 Species_G;Species_N                                                                                                                                       2.07 Tree_3   
  41 Species_B;Species_D;Species_O                                                                                                                             4.05 Tree_3   
  42 Species_B;Species_O                                                                                                                                       1.69 Tree_3   
    
  extractNodeAges(timeTrees,return_summary = 'mean')

  # A tibble: 14 x 3
     Clade                                                                                                                                                 Mean_Node_Age StdDev_Node_Age
     <chr>                                                                                                                                                         <dbl>           <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         60.0          0.00841
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                             24.6          2.10   
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                           19.8          0.349  
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                              8.59         0.722  
   5 Species_A;Species_F                                                                                                                                            2.07         0.305  
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                        6.28         0.482  
   7 Species_B;Species_D;Species_O                                                                                                                                  4.23         0.189  
   8 Species_B;Species_O                                                                                                                                            2.15         0.427  
   9 Species_C;Species_H                                                                                                                                           27.8          0.0149 
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                              9.54         1.35   
  11 Species_G;Species_I;Species_N                                                                                                                                  4.74         0.813  
  12 Species_G;Species_N                                                                                                                                            2.50         0.797  
  13 Species_J;Species_M                                                                                                                                            2.76         0.191  
  14 Species_K;Species_L                                                                                                                                            2.22         0.350  

  extractNodeAges(timeTrees,return_summary = 'median')

  # A tibble: 14 x 3
     Clade                                                                                                                                                 Median_Node_Age  MAD_Node_Age
     <chr>                                                                                                                                                           <dbl>         <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O           60.0  0.00000000593
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                               24.4  2.69         
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                             19.7  0.355        
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                                8.41 0.656        
   5 Species_A;Species_F                                                                                                                                              2.24 0.00379      
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                          6.05 0.151        
   7 Species_B;Species_D;Species_O                                                                                                                                    4.23 0.270        
   8 Species_B;Species_O                                                                                                                                              2.24 0.424        
   9 Species_C;Species_H                                                                                                                                             27.8  0.0173       
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                                9.24 1.27         
  11 Species_G;Species_I;Species_N                                                                                                                                    4.28 0.0248       
  12 Species_G;Species_N                                                                                                                                              2.07 0.100        
  13 Species_J;Species_M                                                                                                                                              2.82 0.132        
  14 Species_K;Species_L                                                                                                                                              2.09 0.216        

  extractNodeAges(timeTrees,return_summary = 'both')

  # A tibble: 14 x 5
     Clade                                                                                                                                                 Mean_Node_Age Median_Node_Age StdDev_Node_Age  MAD_Node_Age
     <chr>                                                                                                                                                         <dbl>           <dbl>           <dbl>         <dbl>
   1 Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O         60.0            60.0          0.00841 0.00000000593
   2 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O                             24.6            24.4          2.10    2.69         
   3 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                                           19.8            19.7          0.349   0.355        
   4 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                                              8.59            8.41         0.722   0.656        
   5 Species_A;Species_F                                                                                                                                            2.07            2.24         0.305   0.00379      
   6 Species_A;Species_F;Species_K;Species_L                                                                                                                        6.28            6.05         0.482   0.151        
   7 Species_B;Species_D;Species_O                                                                                                                                  4.23            4.23         0.189   0.270        
   8 Species_B;Species_O                                                                                                                                            2.15            2.24         0.427   0.424        
   9 Species_C;Species_H                                                                                                                                           27.8            27.8          0.0149  0.0173       
  10 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                                              9.54            9.24         1.35    1.27         
  11 Species_G;Species_I;Species_N                                                                                                                                  4.74            4.28         0.813   0.0248       
  12 Species_G;Species_N                                                                                                                                            2.50            2.07         0.797   0.100        
  13 Species_J;Species_M                                                                                                                                            2.76            2.82         0.191   0.132        
  14 Species_K;Species_L                                                                                                                                            2.22            2.09         0.350   0.216  
