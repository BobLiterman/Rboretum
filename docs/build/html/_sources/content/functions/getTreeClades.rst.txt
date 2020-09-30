.. _getTreeClades:

##################
**getTreeClades**
##################

*getTreeClades* breaks rooted phylo and multiPhylo objects down into their respective clade sets. 

  - If a phylo or multiPhylo is passed, a character vector is returned with each unique clade (from all trees if a multiPhylo is provided), represented as a semicolon-delimited string
  - multiPhylo objects are trimmed down to common taxa prior to clade breakdown
  - Results can be filtered to:
  
    - Omit root clades (include_root = FALSE)
    - Return only shared clades (return_shared = TRUE)
    - Return a table of results with information about the splits in all trees (return_counts = TRUE)
  

=======================
Function and Arguments
=======================

**Usage**:
::

  getTreeClades <- function(tree,include_root,print_counts,return_counts,return_shared)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**tree**				                  A rooted phylo object, or a rooted multiPhylo object where all trees share 3+ taxa
**include_root**                  If FALSE, exclude root clades in the returned results [Default: TRUE, return root clades]
**print_counts**                  If TRUE and a multiPhylo is provided return clades, but also print clade information table to console [Default: FALSE, no printing]
**return_counts**                 If TRUE and a multiPhylo is provided, return clade information table rather than clades [Default: FALSE, return clades]
**return_shared**                 If TRUE and a multiPhlyo is provided, return only clades found in all trees
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getTreeClades.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()

  # Read in trees with two topologies (1=2, 1!=3)

  myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
  myTree_2 <- readRooted(rb_tree2_path,root_taxa = c('Species_C','Species_H'))
  myTree_3 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

  # Toplogy A

  getTreeClades(myTree_1)

  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N"                              
  [3] "Species_A;Species_E;Species_F;Species_K;Species_L"                                                                                
  [4] "Species_A;Species_F"                                                                                                              
  [5] "Species_A;Species_F;Species_K;Species_L"                                                                                          
  [6] "Species_B;Species_D;Species_O"                                                                                                    
  [7] "Species_B;Species_O"                                                                                                              
  [8] "Species_C;Species_H"                                                                                                              
  [9] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                                                
  [10] "Species_G;Species_I;Species_N"                                                                                                    
  [11] "Species_G;Species_N"                                                                                                              
  [12] "Species_J;Species_M"                                                                                                              
  [13] "Species_K;Species_L"                                                                                                              

  getTreeClades(myTree_2)

  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N"                              
  [3] "Species_A;Species_E;Species_F;Species_K;Species_L"                                                                                
  [4] "Species_A;Species_F"                                                                                                              
  [5] "Species_A;Species_F;Species_K;Species_L"                                                                                          
  [6] "Species_B;Species_D;Species_O"                                                                                                    
  [7] "Species_B;Species_O"                                                                                                              
  [8] "Species_C;Species_H"                                                                                                              
  [9] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                                                
  [10] "Species_G;Species_I;Species_N"                                                                                                    
  [11] "Species_G;Species_N"                                                                                                              
  [12] "Species_J;Species_M"                                                                                                              
  [13] "Species_K;Species_L"                                                                                                              

  # Topology B

  getTreeClades(myTree_3)

  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_A;Species_E;Species_F"                                                                                                    
  [3] "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_M;Species_N"                                                  
  [4] "Species_A;Species_F"                                                                                                              
  [5] "Species_B;Species_D;Species_K;Species_L;Species_O"                                                                                
  [6] "Species_B;Species_O"                                                                                                              
  [7] "Species_C;Species_H"                                                                                                              
  [8] "Species_D;Species_K;Species_L"                                                                                                    
  [9] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                                                
  [10] "Species_G;Species_I;Species_N"                                                                                                    
  [11] "Species_G;Species_N"                                                                                                              
  [12] "Species_J;Species_M"                                                                                                              
  [13] "Species_K;Species_L"                                                                                                              

  # Get splits from a multiPhylo

  myMultiPhylo <- readRooted(c(rb_tree1_path,rb_tree2_path,rb_tree4_path),root_taxa = c('Species_C','Species_H'))

  # Return all identified splits

  getTreeClades(myMultiPhylo)

  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_A;Species_E;Species_F"                                                                                                    
  [3] "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N"                              
  [4] "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_M;Species_N"                                                  
  [5] "Species_A;Species_E;Species_F;Species_K;Species_L"                                                                                
  [6] "Species_A;Species_F"                                                                                                              
  [7] "Species_A;Species_F;Species_K;Species_L"                                                                                          
  [8] "Species_B;Species_D;Species_K;Species_L;Species_O"                                                                                
  [9] "Species_B;Species_D;Species_O"                                                                                                    
  [10] "Species_B;Species_O"                                                                                                              
  [11] "Species_C;Species_H"                                                                                                              
  [12] "Species_D;Species_K;Species_L"                                                                                                    
  [13] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                                                
  [14] "Species_G;Species_I;Species_N"                                                                                                    
  [15] "Species_G;Species_N"                                                                                                              
  [16] "Species_J;Species_M"                                                                                                              
  [17] "Species_K;Species_L"                                                                                                              

  length(getTreeClades(myMultiPhylo))

  [1] 17

  # Return all identified splits, but exclude root split

  getTreeClades(myMultiPhylo,include_root = FALSE)

  [1] "Species_A;Species_E;Species_F"                                                                       "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N"
  [3] "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_M;Species_N"                     "Species_A;Species_E;Species_F;Species_K;Species_L"                                                  
  [5] "Species_A;Species_F"                                                                                 "Species_A;Species_F;Species_K;Species_L"                                                            
  [7] "Species_B;Species_D;Species_K;Species_L;Species_O"                                                   "Species_B;Species_D;Species_O"                                                                      
  [9] "Species_B;Species_O"                                                                                 "Species_D;Species_K;Species_L"                                                                      
  [11] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                   "Species_G;Species_I;Species_N"                                                                      
  [13] "Species_G;Species_N"                                                                                 "Species_J;Species_M"                                                                                
  [15] "Species_K;Species_L"                                                                                

  length(getTreeClades(myMultiPhylo,include_root = FALSE))

  [1] 15

  # Return only splits present in all trees

  getTreeClades(myMultiPhylo,return_shared = TRUE)

  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_A;Species_F"                                                                                                              
  [3] "Species_B;Species_O"                                                                                                              
  [4] "Species_C;Species_H"                                                                                                              
  [5] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                                                
  [6] "Species_G;Species_I;Species_N"                                                                                                    
  [7] "Species_G;Species_N"                                                                                                              
  [8] "Species_J;Species_M"                                                                                                              
  [9] "Species_K;Species_L"                                                                                                              

  length(getTreeClades(myMultiPhylo,return_shared = TRUE))

  [1] 9

  # Return a table of results

  getTreeClades(myMultiPhylo,return_counts = TRUE)

  # A tibble: 17 x 3
   Clade                                                                                                                             Count Trees                           
   <chr>                                                                                                                             <int> <chr>                           
  1 Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O     3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  2 Species_A;Species_F                                                                                                                   3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  3 Species_B;Species_O                                                                                                                   3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  4 Species_C;Species_H                                                                                                                   3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  5 Species_G;Species_I;Species_J;Species_M;Species_N                                                                                     3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  6 Species_G;Species_I;Species_N                                                                                                         3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  7 Species_G;Species_N                                                                                                                   3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  8 Species_J;Species_M                                                                                                                   3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  9 Species_K;Species_L                                                                                                                   3 Gene_1.nwk;Gene_2.nwk;Gene_4.nwk
  10 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N                                   2 Gene_1.nwk;Gene_2.nwk           
  11 Species_A;Species_E;Species_F;Species_K;Species_L                                                                                     2 Gene_1.nwk;Gene_2.nwk           
  12 Species_A;Species_F;Species_K;Species_L                                                                                               2 Gene_1.nwk;Gene_2.nwk           
  13 Species_B;Species_D;Species_O                                                                                                         2 Gene_1.nwk;Gene_2.nwk           
  14 Species_A;Species_E;Species_F                                                                                                         1 Gene_4.nwk                      
  15 Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_M;Species_N                                                       1 Gene_4.nwk                      
  16 Species_B;Species_D;Species_K;Species_L;Species_O                                                                                     1 Gene_4.nwk                      
  17 Species_D;Species_K;Species_L                                                                                                         1 Gene_4.nwk                      
  
