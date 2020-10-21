.. _treeTrimmer:

###############
**treeTrimmer**
###############

**treeTrimmer** is an ape wrapper, equivalent to:
::

  ape::drop.tip(tree,taxa_to_remove)

- **treeTrimmer**  can accept a phylo or multiPhylo *tree* argument for pruning
- User can set whether to keep or remove specified *taxa* with the *remove* argument
- If *tree* is a multiPhylo and no other arguments are given, the default behavior is to trim all trees down to those taxa shared among all trees.

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
**taxa**					                Character vector (or semicolon-separated list) of desired tip labels to keep [default] or discard [if *remove* = TRUE]
**remove**                        OPTIONAL: If TRUE, tip labels specified by *taxa* are removed from trees rather than retained [Default: FALSE, prune *tree* down to *taxa*]
===========================      ===============================================================================================================================================================================================================

================
Function Return
================

- If a single tree is provided via *tree*, **treeTrimmer** returns a phylo object pruned down to the desired taxa specified by *taxa* and *remove*
- If multiple trees are provided via *tree*, **treeTrimmer** returns a multiPhylo where all trees are pruned down to the desired taxa specified by *taxa* and *remove*
- If multiple trees are provided via *tree* and no *taxa* are specified, **treeTrimmer** returns a multiPhylo where all trees are pruned down to the common set of taxa among trees
- Function will throw an error if:

  - Fewer than 3 species are specified for retention in any tree
  - Any trees cannot be trimmed down to the final species lists (e.g. taxa are missing)
  
  
==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/treeTrimmer.R

  library(Rboretum)
  sourceRboretum()
  
  # Read in a single tree rooted at the clade of Species C + Species H
  is.rooted(ape::read.tree(rb_tree1_path))
  [1] FALSE

  myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))
  is.rooted(myTree)
  [1] TRUE
  
  # Trim myTree, keeping only Species A - I 
  taxa_to_keep <- c("Species_A","Species_B","Species_C","Species_D","Species_E","Species_F","Species_G","Species_H","Species_I")
  taxa_to_remove <- c("Species_J","Species_K","Species_L","Species_M","Species_N","Species_O")
  
  # Trim myTree by supplying a list of taxa to keep
  myTrimmedTree_keep <- treeTrimmer(tree = myTree,taxa = taxa_to_keep)
  
  # Trim myTree by supplying a list of taxa to remove
  myTrimmedTree_remove <- treeTrimmer(tree = myTree,taxa = taxa_to_remove,remove = TRUE)
  
  # Check tip labels
  naturalsort(myTrimmedTree_keep$tip.label)
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"
  
  naturalsort(myTrimmedTree_remove$tip.label)
  [1] "Species_A" "Species_B" "Species_C" "Species_D" "Species_E" "Species_F" "Species_G" "Species_H" "Species_I"
  
  # Read in all trees from a directory that start with "Gene" and end with ".nwk", and root at the clade of Species C + Species H
  myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),prefix="Gene",suffix=".nwk",dummy_names = TRUE)
  purrr::map(.x=myTrees,.f=function(x){is.rooted(x)}) %>% unlist()
  
  Tree_1 Tree_2 Tree_3 Tree_4 Tree_5 
  TRUE   TRUE   TRUE   TRUE   TRUE 
  
  # Create a multiPhlyo where trees do not share all taxa (one tree only has Species_A - Species_I)
  mixed_trees <- c(myTrees,myTrimmedTree_keep)
  
  # Default behavior is to trim down to common taxa
  myTrimmedTrees_mixed <- treeTrimmer(tree=mixed_trees)
  
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

  # Trim all trees in the multiPhylo down to a given a list of taxa
  myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),prefix="Gene",suffix=".nwk",dummy_names = TRUE)
  myTrimmedTrees <- treeTrimmer(tree=myTrees,taxa='Species_A;Species_B;Species_C;Species_H')
  
  # Check multiPhylo trimming
  purrr::map(.x = myTrimmedTrees, .f = function(x){naturalsort(x$tip.label)})
  $Tree_1
  [1] "Species_A" "Species_B" "Species_C" "Species_H"
  
  $Tree_2
  [1] "Species_A" "Species_B" "Species_C" "Species_H"
  
  $Tree_3
  [1] "Species_A" "Species_B" "Species_C" "Species_H"
  
  $Tree_4
  [1] "Species_A" "Species_B" "Species_C" "Species_H"
  
  $Tree_5
  [1] "Species_A" "Species_B" "Species_C" "Species_H"
