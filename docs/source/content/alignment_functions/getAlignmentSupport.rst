.. _getAlignmentSupport:

########################
**getAlignmentSupport**
########################
.. |br| raw:: html

   <br />

**getAlignmentSupport** queries the the parsimony-based signal from alignments profiled via **getAlignmentSignal**. Users can query support for: 

  - A specific clade or set of clades
  - All clades in a phylo object
  - All clades from a multiPhylo of trees

Users can filter the output from **getAlignmentSignal** with options including:

  - Flexible handling of missing data
  - Handling of indels
  - Separating results by variation pattern (ie. all sites, all parsimony-informative, only biallelic, only triallelic, etc.)
  
If the output from **getAlignmentSignal** contains data from mutilple alignments, results can be separated or pooled via *separate_signal*

.. warning::
  
  **Note:** Rboretum currently does not support handling of degenerate bases, and classifies them broadly as N/missing. Support for degenerate nucleotides is in progress and set for a future release. 

=======================
Function and Arguments
=======================

**Usage**:
::

  getAlignmentSupport <- function(signal,tree,include_root,clade,dataset_name,max_missing,separate_signal,only_parsinf,include_gap,only_gap,only_biallelic,only_triallelic,only_quadallelic,only_pentallelic,return_integer,return_table,existing_support){

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**signal**				                 The output from **getAlignmentSignal**, run with species specified by *tree* or *clade*
**tree**				                   A rooted Phylo object, or a rooted multiPhylo where all trees share 3+ taxa
**include_root**				           Include support about the root split [Default: TRUE]
**clade**				                   Query support for specific clades by providing them as semicolon-delimted characters (ie. c('Species_A;Species_B','Species_X;Species_Y'))
**dataset_name**				           A label for the support column(s). [Default: Alignment name + m<MISSING>]
**max_missing**				             Only record support from columns where no more than *max_missing* taxa have missing data (N or degenerate bases) [Default: All species - 3]
**separate_signal**				         If TRUE and the output from **getAlignmentSignal** has 2+ alignments, create a separate support column for each alignment [Default: TRUE]
**only_parsinf**				           Only query support from variable sites with no singletons [Default: FALSE]
**include_gap**				             Query support from variable sites where gaps are present [Default: TRUE]
**only_gap**				               Only query support from variable sites where gaps are present [Default: FALSE]
**only_biallelic**				         Only query support from biallelic sites [Default: FALSE]
**only_triallelic**				         Only query support from triallelic sites [Default: FALSE]
**only_quadallelic**				       Only query support from quadallelic sites [Default: FALSE]
**only_pentallelic**				       Only query support from pentallelic sites [Default: FALSE]
**return_integer**				         Instead of a dataframe, return a named integer with values summed across alignments
**return_table**				           Instead of a dataframe, return the output from table() (or a list), where names are semicolon-delimted clades
**existing_support**				       Output from **getAlignmentSupport** run with the same *tree* or *clade* argument
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getAlignmentSupport.R

  library(Rboretum)
  sourceRboretum()
  
  # Read in a single phylogeny
  myTree <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
  
  # Extract the signal from an alignment file
  mySignal <- getAlignmentSignal(alignment_path = rb_align1_path,species_info = myTree)
  
  # Map the signal from a single alignment onto a single tree topology
  myTreeSupport <- getAlignmentSupport(tree = myTree,signal = mySignal,include_root = TRUE)
  
  # Note the column name specifies how many taxa were allowed to have missing data (default: Total - 3)
  myTreeSupport
  
  # Map the signal from a single alignment onto a single tree topology, and only consider sites where one taxon has missing data. Add this to previous data.
  myTreeSupport_m1 <- getAlignmentSupport(tree = myTree,signal = mySignal,include_root = TRUE,max_missing = 1,existing_support = myTreeSupport)
  myTreeSupport_m1
  
  # Map the signal from a multiple alignments onto a single tree topology
  mySignals <- getAlignmentSignal(alignment_path = rb_all_align,species_info = myTree)
  
  # Get tree support, separated by alignment. Treat gaps as missing data. Only include parsimony-informative sites
  myTreeSupport_Sep <- getAlignmentSupport(tree = myTree,signal = mySignals,include_root = TRUE,include_gap = FALSE,only_parsinf = TRUE)
  
  # Get same tree support, summed across alignments
  myTreeSupport_Sum <- getAlignmentSupport(tree = myTree,signal = mySignals,include_root = TRUE,separate_signal = FALSE,include_gap = FALSE,only_parsinf = TRUE)
  
  myTreeSupport_Sep
  myTreeSupport_Sum
  
  # Get support for a specific
  myTreeSupport_Clade <- getAlignmentSupport(clade='Species_C;Species_H',signal = mySignals)
  myTreeSupport_Clade
  
  # Get support for a specific clade not found in any trees
  myTreeSupport_Clade <- getAlignmentSupport(clade='Species_C;Species_F',signal = mySignals)
  myTreeSupport_Clade
