.. _getSharedTaxa:

##################
**getSharedTaxa**
##################

**getSharedTaxa** returns all tips shared among trees in a multiPhylo

=======================
Function and Arguments
=======================

**Usage**:
::

  getSharedTaxa <- function(trees)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**trees**				                  A multiPhylo object where all trees share 3+ taxa
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getSharedTaxa.R

  library(Rboretum)
  sourceRboretum()

  # Create a multiPhylo where all trees share all taxa

  tree_1 <- ape::rtree(25)
  tree_2 <- ape::rtree(25)
  tree_3 <- ape::rtree(25)
  trees <- c(tree_1,tree_2,tree_3)
  getSharedTaxa(trees)

  [1] "t1"  "t2"  "t3"  "t4"  "t5"  "t6"  "t7"  "t8"  "t9"  "t10" "t11" "t12" "t13" "t14" "t15" "t16" "t17" "t18" "t19" "t20" "t21"
  [22] "t22" "t23" "t24" "t25"

  length(getSharedTaxa(trees))

  [1] 25

  # Create a multiPhylo where all trees share 10 taxa

  tree_1 <- ape::rtree(30)
  tree_2 <- ape::rtree(20)
  tree_3 <- ape::rtree(10)

  getSharedTaxa(c(tree_1,tree_2,tree_3))

  [1] "t1"  "t2"  "t3"  "t4"  "t5"  "t6"  "t7"  "t8"  "t9"  "t10"

  length(getSharedTaxa(c(tree_1,tree_2,tree_3)))

  [1] 10
  
