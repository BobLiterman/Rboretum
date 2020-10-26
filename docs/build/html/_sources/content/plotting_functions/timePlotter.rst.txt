.. _timePlotter:

################
**timePlotter**
################
.. |br| raw:: html

   <br />

**timePlotter** takes as main arguments: 

  - A dataframe of node age information (ie. the output from **extractNodeAges()**)
  - A dataframe of node site support data (ie. the output from **getAlignmentSupport()**)

Taken together, **timePlotter** calculates the proportion of node support deriving from each alignment, and plots proportional support over time. |br| 

  - Signicance of time effects is determined using linear models via the function **time_lm.R**
  - Results are interpreted at a user-set alpha (default: 0.05) and Bonferroni correction is applied based on the number of genes/subsets analyzed
  - Genes/subsets with a signficant change in node support over time will have filled geoms

=======================
Function and Arguments
=======================

**Usage**:
::

  timePlotter <- function(node_age_df,tree_support,plot_datasets,all_sites_col,lm_alpha,wrap,wrap_scales,return_stats){

===========================      ===================================================================================================================================================================================================================================
 Argument                         Description
===========================      ===================================================================================================================================================================================================================================
**node_age_df**				            Output from extractNodeAges(phylo) or extractNodeAges(multiPhylo,return_summary = 'mean' or 'median')
**tree_support**                  Output from getAlignmentSupport with information about all clades in node_age_df and 2+ alignment support columns
**plot_datasets**                 OPTIONAL: Character vector specifying the datasets desired for the plot (Must be 2+ and present in 'tree_support')
**all_sites_col**                 OPTIONAL (Very special case): If there is is a column in the support data corresponding to all non-overlapping sites, specify the column name to be used when calculating proportions
**lm_alpha**                      OPTIONAL: Run linear models at this alpha level (must be > 0 and < 1; Alpha will be automatically Bonferroni corrected based on dataset count) [Default: Do not run linear models]
**wrap**                          OPTIONAL: (TRUE/FALSE); Plot datasets in separate facets [Default: FALSE; do not wrap and produce 1 plot]
**wrap_scales**                   OPTIONAL: ('free'/'fixed'): If wrapping, allow free or fixed scales on each facet [Default: 'fixed']
**return_stats**                  OPTIONAL: (TRUE/FALSE): In addition to returning the plot, also return the statistical output from linear models in a list
===========================      ===================================================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  library(Rboretum)
  sourceRboretum()

  # Read in ultrametric timetrees (branch lengths ~ time)
  timeTrees <- readRooted(c(rb_timeTree1_path,rb_timeTree2_path,rb_timeTree3_path),root_taxa = c('Species_C','Species_H')) %>% treeNamer()

  # Extract dates from 3 timetrees and get the mean node ages
  tree_dates <- extractNodeAges(timeTrees,return_summary = 'mean')

  # Read in signal from 5 multiple sequence alignments
  mySignals <- getAlignmentSignal(rb_alignment_dir,species_info = timeTrees,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))
  mySupports <- getAlignmentSupport(mySignals,timeTrees,include_root = TRUE,dataset_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

  # Plot percent breakdown of split support by alignment over time, looking for time effect at alpha = 0.05 (pre-Bonferroni correction)
  # Ie. For each node in the tree, what percentage of support comes from each gene/subset?
  timePlotter(node_age_df = tree_dates,tree_support = mySupports,lm_alpha = 0.05,return_stats = TRUE)
  
    $Stats
          Slope     StdErr    t_value   p_value Dataset BF_Sig    Adj_R_Sq Intercept
  1  0.11845685 0.08182761  1.4476392 0.1783318  Gene_A      N  0.09058284  6.644187
  2  0.07014702 0.24570274  0.2854955 0.7810906  Gene_B      N -0.09110665  35.72811
  3 -0.18168713 0.11408523 -1.5925562 0.1423441  Gene_C      N   0.1225436  26.54905
  4 -0.03640601 0.20012941 -0.1819123 0.8592861  Gene_D      N -0.09637188  22.50165
  5  0.02948928 0.11397543  0.2587337 0.8010884  Gene_E      N -0.09268523  8.576999

  # Note: In this case, there is not a significant effect of time, so no geoms are filled in 

.. image:: ../images/timePlotter.png 
  :width: 600