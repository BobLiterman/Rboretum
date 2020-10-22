.. _tableCount:

###############
**tableCount**
###############

**tableCount** queries an R table for a specific item, and returns the count of that item. If the item does not occur in the table, tableCount() returns 0.

=======================
Function and Arguments
=======================

**Usage**:
::

  tableCount(search_table,name)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**search_table**				          R Table (Result of 'table()' call)
**name**                          Name of item you want the count for
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/tableCount.R

  library(Rboretum)

  myVector <- c('a','b','b','c','c','c')
  myTable <- table(myVector)

  myTable[['a']]
  [1] 1

  tableCount(search_table = myTable,name = 'a')
  [1] 1

  myTable[['b']]
  [1] 2

  tableCount(myTable,'b')
  [1] 2

  myTable[['c']]
  [1] 3

  tableCount(myTable,'c')
  [1] 3

  myTable[['d']]
  Error in myTable[["d"]] : subscript out of bounds

  tableCount(myTable,'d')
  [1] 0
