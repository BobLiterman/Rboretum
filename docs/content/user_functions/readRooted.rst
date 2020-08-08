.. _readRooted:

###############
**readRooted**
###############

*readRooted* is an ape wrapper, equivalent to:
::

  ape::read.tree(TREE) %>% ape::unroot.phylo(.) %>% ape::root.phylo(.,ROOT_TAXA,resolve_root=TRUE)


*readRooted* can be run in with a single file path or in batch mode (simultaneously reading in and rooting multiple trees with the same outgroup taxa).

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
**dummy_names**                   OPTIONAL: If TRUE, and multiple tree paths are provdied, trees will be named with placeholder names (e.g. Tree_1, Tree_2, etc.)
**prefix**	                      OPTIONAL: If 'to_root' is a directory, provide a character vector of file prefixes (e.g. all trees start with "RAxML")
**suffix**	                      OPTIONAL: If 'to_root' is a directory, provide a character vector of file suffixes (e.g. all trees end with ".nwk")
===========================      ===============================================================================================================================================================================================================

================
Function Return
================

- If a single tree is provided via **to_root**, *readRooted* returns a phylo object rooted at **root_taxa**
- If multiple trees are provided via **to_root**, *readRooted* returns a named multiPhylo object, with all trees rooted at **root_taxa**
