.. _semiSorter:

###############
**semiSorter**
###############

**semiSorter** can accept: 

  - A character vector of unjoined elements
  - A semicolon-delimited character of elements
  - A character vector of semicolon-delimited elements

In any case, **semiSorter** returns the same elements, sorted via *naturalsort*, and joined by semicolons (";"). 

  - If a vector of elements is passed, a single character element is returned.
  - If a set of semicolon-delimted elements are given, a character vector is returned containing each set, sorted

=======================
Function and Arguments
=======================

**Usage**:
::

  semiSorter(string_to_sort)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**string_to_sort**				        Character vector of elements to join, or one or more semicolon-delimited sets
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/semiSorter.R

  library(Rboretum)

  character_vector <- c('c','b','a')
  semiSorter(character_vector)
  [1] "a;b;c"

  semicolon_character <- 'c;b;a'
  semiSorter(semicolon_character)
  [1] "a;b;c"

  semicolon_characters <- c('c;b;a','f;e;d')
  semiSorter(semicolon_characters)
  [1] "a;b;c" "d;e;f"
  
