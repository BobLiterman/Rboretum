########################
**Installing Rboretum**
########################

==============
Prerequisites
==============

1. R (Version 3.6.1+)
2. Python (Version 3.6+) with Biopython
3. A configured version of *reticulate* pointing to the Python3.6+ with Biopython installed

=============
Installation
=============

**Rboretum** can be installed through devtools::
  
    devtools::install_github('BobLiterman/Rboretum')

.. note::

  As of the time of this README, there is a bug in *ape* v.5.4 that breaks functionality in is.monophyletic(), so **Rboretum** requires v.5.3.
  
