.. _vectorSemi:

###############
**vectorSemi**
###############

**vectorSemi** is equivalent to:
::

  paste(to_semi,collapse=";")

**vectorSemi** will convert:

  - A multi-element character vector into a single, semicolon-delimted character
  - A list of character vectors into a list of semicolon-delimted characters

=======================
Function and Arguments
=======================

**Usage**:
::

  vectorSemi(to_semi)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**to_semi**				                Charcter vector, or a list of character vectors
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/vectorSemi.R

  library(Rboretum)

  character_vector <- c('a','b','c')
  vectorSemi(character_vector)
  [1] "a;b;c"

  character_list <- list(c('a','b','c'),c('d','e','f'))
  vectorSemi(character_list)
  [[1]]
  [1] "a;b;c"

  [[2]]
  [1] "d;e;f"
  
