.. _checkTips:

##############
**checkTips**
##############

*checkTips* queries phylo or multiPhylo objects regarding the presence, monophyly, and/or root nature of a user-supplied set of tips

=======================
Function and Arguments
=======================

**Usage**:
::

  checkTips <- function(tree,taxa,check_mono,check_root)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**tree**				                  A phylo object, or a multiPhylo where all trees share 3+ taxa
**taxa**                          Tip labels to query; Can be provided as a vector of IDs, or a semicolon-delimited character
**check_mono**                    If TRUE, check if **taxa** are monophyletic in **tree**. Applies to all trees in a multiPhylo. [Default: FALSE]
**check_root**                    If TRUE, check if **taxa** form one of the two tree root clades. Applies to all trees in a multiPhylo. [Default: FALSE]
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/checkTips.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()

  myTree_1 <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
  myTree_2 <- readRooted(rb_tree4_path,root_taxa = c('Species_C','Species_H'))

  # Set some clades of interest

  coi_1 <- c('Species_J','Species_M')
  coi_2 <- c('Species_C','Species_H')
  coi_3 <- c('Species_D','Species_K','Species_L')

  # Visualize trees

  coi_list <- list("One"=coi_1,"Two"=coi_2,"Three"=coi_3)

  treePlotter(myTree_1,branch_weight = 1.5,xmax=10,to_color = coi_list,taxa_fontface = "bold")
  treePlotter(myTree_2,branch_weight = 1.5, xmax=10,to_color = coi_list,taxa_fontface = "bold")
  

.. image:: ../images/checkTips_1.png
  :width: 400
  :alt: treePlotter(myTree_1,xmax=10)


.. image:: ../images/checkTips_2.png
  :width: 400
  :alt: treePlotter(myTree_2,xmax=10)

.. code-block:: r

  # Check if taxa from Clade of Interest 1 are in the trees

  checkTips(tree=myTree_1,taxa = coi_1)
  [1] TRUE

  checkTips(tree=myTree_2,taxa = coi_1)
  [1] TRUE

  # Check if taxa from Clade of Interest 1 form a monophyletic group

  checkTips(tree=myTree_1,taxa = coi_1,check_mono = TRUE)
  [1] TRUE

  checkTips(tree=myTree_2,taxa = coi_1,check_mono = TRUE)
  [1] TRUE

  # Check if taxa from Clade of Interest 1 form a monophyletic group at the root of the tree

  checkTips(tree=myTree_1,taxa = coi_1,check_mono = TRUE, check_root=TRUE)
  [1] FALSE

  checkTips(tree=myTree_2,taxa = coi_1,check_mono = TRUE, check_root=TRUE)
  [1] FALSE

  # Check if taxa from Clade of Interest 2 form a monophyletic group at the root of the tree

  checkTips(tree=myTree_1,taxa = coi_2,check_mono = TRUE, check_root=TRUE)
  [1] TRUE

  checkTips(tree=myTree_2,taxa = coi_2,check_mono = TRUE, check_root=TRUE)
  [1] TRUE

  # Check if taxa from Clade of Interest 3 form a monophyletic group

  checkTips(tree=myTree_1,taxa = coi_3,check_mono = TRUE)
  [1] FALSE

  checkTips(tree=myTree_2,taxa = coi_3,check_mono = TRUE)
  [1] TRUE
