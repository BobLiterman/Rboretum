.. _semiVector:

###############
**semiVector**
###############

Rboretum uses semicolon (;) separated characters for many functions. 

*semiVector* converts: 
- A single-element semicolon-separated character into a character vector split by ";"
- A multi-element semicolon-separated character into list of character vectors split by ";"

=======================
Function and Arguments
=======================

**Usage**:
::

  semiVector(string_to_split)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**string_to_split**				        Semicolon delimited character (can contain 1 or more elements)
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/semiVector.R
  
  library(Rboretum)
  
  semicolon_single_char <- 'a;b;c;d;e;f'
  semiVector(semicolon_single_char)
  [1] "a" "b" "c" "d" "e" "f"
  
  semicolon_multi_char <- c('a;b;c','d;e;f')
  semiVector(semicolon_multi_char)
  [[1]]
  [1] "a" "b" "c"

  [[2]]
  [1] "d" "e" "f"

  