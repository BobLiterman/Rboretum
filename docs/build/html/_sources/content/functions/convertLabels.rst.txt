.. _convertLabels:

##################
**convertLabels**
##################

*convertLabels* uses a dataframe/tibble of name equivalencies and can use it to convert:

  - Tree tip names
  - Names in the column of a dataframe/tibble (e.g. output from getTreeSupport or getTreeClades)
  - Semicolon-delimited clade groups

=======================
Function and Arguments
=======================

**Usage**:
::

  convertLabels <- function(to_convert,name_df,from,to)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**to_convert**				            A phylo, multiPhylo, character vector, or list of IDs to convert
**name_df**                       A dataframe/tibble with column names that has the name equivalencies
**from**                          Column name of current IDs [Default: First column]
**to**                            Column name of desired IDs [Default: Second column]
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/convertLabels.R
  
  library(Rboretum)
  
  # Set test data directory
  sourceRboretum()
  

  