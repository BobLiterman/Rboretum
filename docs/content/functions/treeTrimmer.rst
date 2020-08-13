.. _treeTrimmer:

###############
**treeTrimmer**
###############

*treeTrimmer* is an ape wrapper, equivalent to:
::

  ape::drop.tip(tree,taxa_to_remove)

- *treeTrimmer*  can accept a phylo or multiPhylo argument for pruning
- User can set whether to keep or remove specified taxa
- With a multiPhylo **tree** and no other arguments, the default behavior is to trim all trees down to those taxa shared among all trees

=======================
Function and Arguments
=======================

**Usage**:
::

  treeTrimmer <- function(tree,taxa,remove)
  

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**tree**				                  Tree(s) to trim. Options include: (1) A phylo or (2) a multiPhylo object where all trees share 3+ tip labels 
**taxa**					                Character vector (or semicolon-separated list) of desired tip labels to keep (or discard if remove=TRUE)
**remove**                        OPTIONAL: If TRUE, tip labels specified by **taxa** are removed from trees rather than retained [Default: FALSE, prune **tree** down to **taxa**]
===========================      ===============================================================================================================================================================================================================

================
Function Return
================

- If a single tree is provided via **tree**, *treeTrimmer* returns a phylo object pruned down to the desired taxa specified by **taxa** and **remove**
- If multiple trees are provided via **tree**, *treeTrimmer* returns a multiPhylo where all trees are pruned down to the desired taxa specified by **taxa** and **remove**
- Function will throw an error if:

  - Fewer than 3 species are specified for retention in any tree
  - Any trees cannot be trimmed down to the final species lists (e.g. taxa are missing)
  
  
==============
Example Usage
==============

::
  
  library(Rboretum)

  # Set test data directory
  
  sourceRboretum()

  # Read in a single tree and multiple trees, rooted at the clade of Species C + Species H
  
  myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))
  naturalsort(myTree$tip.label)

  myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk")

  # Trim myTree down to Species A - I 
  
  taxa_to_keep <- c("Species_A","Species_B","Species_C","Species_D","Species_E","Species_F","Species_G","Species_H","Species_I")
  taxa_to_remove <- c("Species_J","Species_K","Species_L","Species_M","Species_N","Species_O")

  # Trim myTree by supplying a list of taxa to keep
  
  myTrimmedTree_keep <- treeTrimmer(tree = myTree,taxa = taxa_to_keep)
  naturalsort(myTrimmedTree_keep$tip.label)
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I" "Species_J" "Species_K" "Species_L" "Species_M" "Species_N" "Species_O"

  # Trim myTree by supplying a list of taxa to remove
  
  myTrimmedTree_remove <- treeTrimmer(tree = myTree,taxa = taxa_to_remove,remove = TRUE)
  naturalsort(myTrimmedTree_remove$tip.label)
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  # Test default behavior with mutltiPhylo (trim to common taxa)
  
  mixed_trees <- c(myTrees,myTrimmedTree_keep)
  myTrimmedTrees_mixed <- treeTrimmer(tree=mixed_trees)
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  # Check multiPhylo trimming
  
  purrr::map(.x = myTrimmedTrees_mixed, .f = function(x){naturalsort(x$tip.label)})
  $Tree_1
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Tree_2
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Tree_3
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Tree_4
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Tree_5
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Tree_6
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"
  
  # Trim multiPhylo given a list of taxa
  
  myTrimmedTrees <- treeTrimmer(tree=myTrees,taxa=taxa_to_keep)

  # Check multiPhylo trimming
  
  purrr::map(.x = myTrimmedTrees, .f = function(x){naturalsort(x$tip.label)})
  $Gene_1.nwk
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Gene_2.nwk
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Gene_3.nwk
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Gene_4.nwk
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"

  $Gene_5.nwk
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"
  
