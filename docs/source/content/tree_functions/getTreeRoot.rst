.. _getTreeRoot:

##################
**getTreeRoot**
##################

**getTreeRoot** takes a rooted phylo object and returns the two root splits

=======================
Function and Arguments
=======================

**Usage**:
::

  getTreeRoot <- function(tree)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**tree**				                  A rooted phylo object
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getTreeRoot.R

  library(Rboretum)
  sourceRboretum()
  
  # Read in rooted tree
  myTree <- readRooted(rb_tree1_path,'Species_C;Species_H')
  getTreeRoot(myTree)
  [1] "Species_C;Species_H"                                                                                                              
  [2] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"

  myTree <- readRooted(rb_tree1_path,'Species_B;Species_O')
  getTreeRoot(myTree)
  [1] "Species_B;Species_O"                                                                                                              
  [2] "Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N"

  myTree <- readRooted(rb_tree1_path,'Species_F')
  getTreeRoot(myTree)
  [1] "Species_F"                                                                                                                                  
  [2] "Species_A;Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
