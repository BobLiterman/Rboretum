.. _treeNamer:

###############
**treeNamer**
###############

*treeNamer* takes a multiPhylo object, and returns the same trees with 'dummy names' (e.g. Tree_1,Tree_2,etc.)

=======================
Function and Arguments
=======================

**Usage**:

::
  
  treeNamer <- function(trees)
  

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**trees**				                  A multiPhylo object 
===========================      ===============================================================================================================================================================================================================
  
==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/treeNamer.R
  
  library(Rboretum)

  # Set test data directory
  sourceRboretum()
  
  # Read in multiPhylo
  noname_multiPhylo <- readRooted(rb_all_unrooted,root_taxa = c('Species_C','Species_H'))
  names(noname_multiPhylo) # Names are based on filenames by default
  [1] "Gene_1.nwk" "Gene_2.nwk" "Gene_3.nwk" "Gene_4.nwk" "Gene_5.nwk"
  
  # Add dummy names
  named_multiPhlyo <- treeNamer(noname_multiPhylo)
  names(named_multiPhlyo)
  [1] "Tree_1" "Tree_2" "Tree_3" "Tree_4" "Tree_5"

