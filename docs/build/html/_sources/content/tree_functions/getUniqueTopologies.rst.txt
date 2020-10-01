.. _getUniqueTopologies:

########################
**getUniqueTopologies**
########################

*getUniqueTopologies* takes a multiPhlyo, prunes down to common taxa, and returns unique topologies or summary information

=======================
Function and Arguments
=======================

**Usage**:
::

  getUniqueTopologies <- function(trees,tree_names,print_table,return_table)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**trees**				                  A multiPhylo object where all trees share 3+ taxa
**tree_names**                    If TRUE, use names of input trees when naming returned phylogenies (e.g. Tree Name = Gene1;Gene2;Gene3) [Default: FALSE; Return trees with names Topology_1,Topology_2,Topology_3, etc. ]
**print_table**                   If TRUE, print summary table and return dataframe [Default: FALSE, no printing]
**return_table**                  If TRUE, return summary table instead of phylo/multiPhylo
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getUniqueTopologies.R

  library(Rboretum)

  # Set test data directory
  sourceRboretum()

  # Read in multiPhylo 
  myMultiPhylo <- readRooted(to_root = rb_all_unrooted,root_taxa = c('Species_C','Species_H'))

  myMultiPhylo
  5 phylogenetic trees

  treePlotter(myMultiPhylo,basic_plot = TRUE,xmax=10,taxa_font_size = 3,node_label_font_size = 3)

.. image:: ../images/getUniqueTopologies_1.png
  :width: 800
  :alt: treePlotter(myMultiPhylo,basic_plot = TRUE,xmin=-0.5,xmax=10,taxa_font_size = 3,node_label_font_size = 3)

.. code-block:: r
  
  # Reduce to unique topologies, use tree names for new trees

  myUniqueTrees_1 <- getUniqueTopologies(trees = myMultiPhylo,tree_names = TRUE,print_table = TRUE)

                           Tree_Name              Trees_with_Topology Tree_Count Tree_Percent
  1 Gene_1.nwk;Gene_2.nwk;Gene_3.nwk Gene_1.nwk;Gene_2.nwk;Gene_3.nwk          3           60
  2                       Gene_4.nwk                       Gene_4.nwk          1           20
  3                       Gene_5.nwk                       Gene_5.nwk          1           20

  myUniqueTrees_1
  3 phylogenetic trees

  treePlotter(myUniqueTrees_1,basic_plot = TRUE,xmax=10,taxa_font_size = 3,node_label_font_size = 3)

.. image:: ../images/getUniqueTopologies_2.png
  :width: 800
  :alt: treePlotter(myUniqueTrees_1,basic_plot = TRUE,xmin=-0.5,xmax=10,taxa_font_size = 3,node_label_font_size = 3)

.. code-block:: r

  # Reduce to unique topolodies, use dummy names for new trees

  myUniqueTrees_2 <- getUniqueTopologies(trees = myMultiPhylo,tree_names = FALSE,print_table = TRUE)

     Tree_Name              Trees_with_Topology Tree_Count Tree_Percent
  1 Topology_1 Gene_1.nwk;Gene_2.nwk;Gene_3.nwk          3           60
  2 Topology_2                       Gene_4.nwk          1           20
  3 Topology_3                       Gene_5.nwk          1           20

  myUniqueTrees_2
  3 phylogenetic trees

  treePlotter(myUniqueTrees_2,basic_plot = TRUE,xmax=10,taxa_font_size = 3,node_label_font_size = 3)

.. image:: ../images/getUniqueTopologies_3.png
  :width: 800
  :alt: treePlotter(myUniqueTrees_2,basic_plot = TRUE,xmin=-0.5,xmax=10,taxa_font_size = 3,node_label_font_size = 3)
