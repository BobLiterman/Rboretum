.. _readRooted:

###############
**readRooted**
###############

*readRooted* is an ape wrapper, equivalent to:
::

  ape::read.tree(TREE) %>% ape::unroot.phylo(.) %>% ape::root.phylo(.,ROOT_TAXA,resolve_root=TRUE)


- *readRooted* can read in Newick and Nexus trees, as well as trees with booststrap values stored as RAxML-style branch labels
- *readRooted* can be run in with a single file or in batch mode (simultaneously reading in and rooting multiple trees with the same outgroup taxa).
- When reading in multiple trees, *readRooted* will name trees based on the file basename, but this can be overridden by supplying names via the **"tree_names"** arguments, or by setting **"dummy_names"** to TRUE.

=======================
Function and Arguments
=======================

**Usage**:
::

  readRooted(to_root,root_taxa,tree_names,dummy_names,prefix,suffix)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**to_root**				                Where to find tree files. Options include: (1) A character vector of one or more tree file paths; or (2) A path to a single directory containing all tree files 
**root_taxa**					            Outgroup species IDs. Must be in tree(s) and monophyletic. Can be provided as: (1) A character vector of one or more tip labels; or (2) A semicolon-separated list of tip labels
**tree_names**                    OPTIONAL: If multiple tree paths are provided, a character vector of names to assign to trees. Length must equal the number of trees. [**Default:** Use basename of tree file as the tree name]
**dummy_names**                   OPTIONAL: If TRUE, and multiple tree paths are provdied, trees will be named with placeholder names (e.g. Tree_1, Tree_2, etc.) [**Default:** FALSE, use file basenames]
**prefix**	                      OPTIONAL: If 'to_root' is a directory, provide a character vector of file prefixes (e.g. all trees start with "RAxML")
**suffix**	                      OPTIONAL: If 'to_root' is a directory, provide a character vector of file suffixes (e.g. all trees end with ".nwk")
===========================      ===============================================================================================================================================================================================================

================
Function Return
================

- If a single tree is provided via **"to_root"**, *readRooted* returns a phylo object rooted at **"root_taxa"**
- If multiple trees are provided via **to_root**, *readRooted* returns a named multiPhylo object, with all trees rooted at **root_taxa**
- Function will throw an error if:
  - Any tree cannot be found based on the path given
  - Any tree cannot be rooted on the specified taxa
  

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/readRooted.R
  
  library(Rboretum)
  
  # Set test data directory
  sourceRboretum()
  
  # Read in a single tree and root at the clade of Species C + Species H
  rb_tree1_path
  [1] "<PACKAGE_DIR>/Rboretum/extdata/unrootedTrees/Gene_1.nwk"
  
  raw_tree <- ape::read.tree(rb_tree1_path)
  ape::is.rooted(raw_tree)
  [1] FALSE
  
  myTree <- readRooted(to_root = rb_tree1_path, root_taxa = c('Species_C','Species_H'))
  ape::is.rooted(myTree
  [1] TRUE
  
  # Read in a multiple unrooted trees and root at the clade of Species C + Species H
  rb_all_unrooted
  [1] "<PACKAGE_DIR>/Rboretum/extdata/unrootedTrees/Gene_1.nwk"
  [2] "<PACKAGE_DIR>/Rboretum/extdata/unrootedTrees/Gene_2.nwk"
  [3] "<PACKAGE_DIR>/Rboretum/extdata/unrootedTrees/Gene_3.nwk"
  [4] "<PACKAGE_DIR>/Rboretum/extdata/unrootedTrees/Gene_4.nwk"
  [5] "<PACKAGE_DIR>/Rboretum/extdata/unrootedTrees/Gene_5.nwk"
  
  purrr::map(.x=rb_all_unrooted,.f=function(x){ape::read.tree(x) %>% ape::is.rooted(.)}) %>% unlist()
  [1] FALSE FALSE FALSE FALSE FALSE
  
  myTrees <- readRooted(to_root = rb_all_unrooted, root_taxa = c('Species_C','Species_H'))
  purrr::map(.x=myTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()
  
  Gene_1.nwk Gene_2.nwk Gene_3.nwk Gene_4.nwk Gene_5.nwk 
        TRUE       TRUE       TRUE       TRUE       TRUE 
  
  # From a directory containing multiple unrooted trees, read in all '.nwk' files and root at the clade of Species C + Species H
  rb_unroot_dir
  [1] "<PACKAGE_DIR>/Rboretum/extdata/unrootedTrees"
  
  myTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk")
  purrr::map(.x=myTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()
  
  Gene_1.nwk Gene_2.nwk Gene_3.nwk Gene_4.nwk Gene_5.nwk 
        TRUE       TRUE       TRUE       TRUE       TRUE 
  
  # Same as above, but add placeholder tree_names ('Tree_1' - 'Tree_5') as opposed to tree file basenames
  myDummyTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk",dummy_names=TRUE)
  purrr::map(.x=myDummyTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()
  
  Tree_1 Tree_2 Tree_3 Tree_4 Tree_5 
    TRUE   TRUE   TRUE   TRUE   TRUE 
  
  # Same as above, but add user-defined tree tree_names as opposed to tree file basenames
  myTreeNames <- c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E')
  
  myNamedTrees <- readRooted(to_root = rb_unroot_dir, root_taxa = c('Species_C','Species_H'),suffix=".nwk",tree_names=myTreeNames)
  purrr::map(.x=myNamedTrees,.f=function(x){ape::is.rooted(x)}) %>% unlist()
  
  Gene_A Gene_B Gene_C Gene_D Gene_E 
  TRUE   TRUE   TRUE   TRUE   TRUE 
  
