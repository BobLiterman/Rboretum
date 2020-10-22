.. _getSpeciesComposition:

##########################
**getSpeciesComposition**
##########################

*getSpeciesComposition* takes one or more alignment files and returns information about the sequence content by species, including base counts, 'N' counts, and gap counts

.. warning::
  
  **Note:** Rboretum currently does not support handling of degenerate bases, and classifies them broadly as N/missing. Support for degenerate nucleotides is in progress and set for a future release. 


=======================
Function and Arguments
=======================

**Usage**:
::

  getSpeciesComposition <- function(alignment_path,species_info,alignment_name,prefix,suffix)

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**alignment_path**				        An absolute or relative path to an alignment file(s) or a directory containing multiple alignments
**species_info**                  A list of species to consider, provided as a phylo, multiPhlyo, character vector or labels or a semicolon-delimited string [Default: Shared species among all alignments]
**alignment_name**                A character vector containing desired alignment IDs [Default: Derive name from filename]
**prefix**                        If **alignment_path** points to a directory, only query files starting with **prefix** (e.g. 'Alignment') [Default: Use all files in directory]
**suffix**                        If **alignment_path** points to a directory, only query files ending wtih **suffix** (e.g. '.nex') [Default: Use all files in directory]
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getSpeciesComposition.R

  library(Rboretum)
  sourceRboretum()
  
  # Set path to alignment data
  myAlignmentFile <- rb_align1_path
  myAlignmentDir <- rb_alignment_dir
  
  # Get species composition information for a single alignment
  getSpeciesComposition(alignment_path = myAlignmentFile)

         Taxon Total_Bases Total_Degenerate Total_N Total_Gaps Percent_GC Percent_Degenerate   Percent_N Percent_Gap Alignment_Name
  1  Species_A        1489                0      10         52  0.4002686                  0 0.006447453  0.03352676     Gene_1.phy
  2  Species_B        1418                0      70         63  0.4055007                  0 0.045132173  0.04061896     Gene_1.phy
  3  Species_C        1457                0      43         51  0.4015100                  0 0.027724049  0.03288201     Gene_1.phy
  4  Species_D        1467                0      26         58  0.3980913                  0 0.016763378  0.03739523     Gene_1.phy
  5  Species_E        1463                0      37         51  0.4005468                  0 0.023855577  0.03288201     Gene_1.phy
  6  Species_F        1496                0       8         47  0.4050802                  0 0.005157963  0.03030303     Gene_1.phy
  7  Species_G        1446                0      54         51  0.3865837                  0 0.034816248  0.03288201     Gene_1.phy
  8  Species_H        1464                0      38         49  0.3989071                  0 0.024500322  0.03159252     Gene_1.phy
  9  Species_I        1454                0      49         48  0.3899587                  0 0.031592521  0.03094778     Gene_1.phy
  10 Species_J        1470                0      37         44  0.3918367                  0 0.023855577  0.02836879     Gene_1.phy
  11 Species_K        1458                0      48         45  0.4080933                  0 0.030947776  0.02901354     Gene_1.phy
  12 Species_L        1458                0      44         49  0.4039781                  0 0.028368794  0.03159252     Gene_1.phy
  13 Species_M        1433                0      70         48  0.3963712                  0 0.045132173  0.03094778     Gene_1.phy
  14 Species_N        1431                0      69         51  0.3878407                  0 0.044487427  0.03288201     Gene_1.phy
  15 Species_O        1482                0      20         49  0.3994602                  0 0.012894907  0.03159252     Gene_1.phy
  
  # Get species composition information from all .phy files in a directory, providing new names
  getSpeciesComposition(alignment_path = myAlignmentDir,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))
  
         Taxon Total_Bases Total_Degenerate Total_N Total_Gaps Percent_GC Percent_Degenerate   Percent_N Percent_Gap Alignment_Name
  1  Species_A        1489                0      10         52  0.4002686                  0 0.006447453  0.03352676         Gene_A
  2  Species_B        1418                0      70         63  0.4055007                  0 0.045132173  0.04061896         Gene_A
  3  Species_C        1457                0      43         51  0.4015100                  0 0.027724049  0.03288201         Gene_A
  4  Species_D        1467                0      26         58  0.3980913                  0 0.016763378  0.03739523         Gene_A
  5  Species_E        1463                0      37         51  0.4005468                  0 0.023855577  0.03288201         Gene_A
  6  Species_F        1496                0       8         47  0.4050802                  0 0.005157963  0.03030303         Gene_A
  7  Species_G        1446                0      54         51  0.3865837                  0 0.034816248  0.03288201         Gene_A
  8  Species_H        1464                0      38         49  0.3989071                  0 0.024500322  0.03159252         Gene_A
  9  Species_I        1454                0      49         48  0.3899587                  0 0.031592521  0.03094778         Gene_A
  10 Species_J        1470                0      37         44  0.3918367                  0 0.023855577  0.02836879         Gene_A
  11 Species_K        1458                0      48         45  0.4080933                  0 0.030947776  0.02901354         Gene_A
  12 Species_L        1458                0      44         49  0.4039781                  0 0.028368794  0.03159252         Gene_A
  13 Species_M        1433                0      70         48  0.3963712                  0 0.045132173  0.03094778         Gene_A
  14 Species_N        1431                0      69         51  0.3878407                  0 0.044487427  0.03288201         Gene_A
  15 Species_O        1482                0      20         49  0.3994602                  0 0.012894907  0.03159252         Gene_A
  16 Species_A        2468                0      18        318  0.5178282                  0 0.006419401  0.11340942         Gene_B
  17 Species_B        2390                0     124        290  0.5196653                  0 0.044222539  0.10342368         Gene_B
        .
        .
        .
        
  # Get species composition information from all .phy files in a directory, providing new names, considering only Species A - E
  getSpeciesComposition(alignment_path = myAlignmentDir,species_info = 'Species_A;Species_B;Species_C;Species_D;Species_E',suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

         Taxon Total_Bases Total_Degenerate Total_N Total_Gaps Percent_GC Percent_Degenerate   Percent_N Percent_Gap Alignment_Name
  1  Species_A        1489                0      10         52  0.4002686                  0 0.006447453  0.03352676         Gene_A
  2  Species_B        1418                0      70         63  0.4055007                  0 0.045132173  0.04061896         Gene_A
  3  Species_C        1457                0      43         51  0.4015100                  0 0.027724049  0.03288201         Gene_A
  4  Species_D        1467                0      26         58  0.3980913                  0 0.016763378  0.03739523         Gene_A
  5  Species_E        1463                0      37         51  0.4005468                  0 0.023855577  0.03288201         Gene_A
  6  Species_A        2468                0      18        318  0.5178282                  0 0.006419401  0.11340942         Gene_B
  7  Species_B        2390                0     124        290  0.5196653                  0 0.044222539  0.10342368         Gene_B
  8  Species_C        2438                0      80        286  0.5258409                  0 0.028530670  0.10199715         Gene_B
  9  Species_D        2459                0      49        296  0.5233835                  0 0.017475036  0.10556348         Gene_B
  10 Species_E        2437                0      67        300  0.5190808                  0 0.023894437  0.10699001         Gene_B
  11 Species_A         735                0       6        290  0.5891156                  0 0.005819593  0.28128031         Gene_C
  12 Species_B         707                0      45        279  0.5714286                  0 0.043646945  0.27061106         Gene_C
  13 Species_C         728                0      29        274  0.5920330                  0 0.028128031  0.26576140         Gene_C
  14 Species_D         739                0      18        274  0.5818674                  0 0.017458778  0.26576140         Gene_C
  15 Species_E         726                0      23        282  0.5950413                  0 0.022308438  0.27352085         Gene_C
  16 Species_A        1992                0      14        213  0.4121486                  0 0.006309148  0.09598918         Gene_D
  17 Species_B        1912                0     100        207  0.4173640                  0 0.045065345  0.09328526         Gene_D
  18 Species_C        1942                0      64        213  0.4150360                  0 0.028841821  0.09598918         Gene_D
  19 Species_D        1964                0      38        217  0.4205703                  0 0.017124831  0.09779180         Gene_D
  20 Species_E        1961                0      53        205  0.4140745                  0 0.023884633  0.09238396         Gene_D
  21 Species_A        1491                0       9          0  0.4875922                  0 0.006000000  0.00000000         Gene_E
  22 Species_B        1437                0      63          0  0.4829506                  0 0.042000000  0.00000000         Gene_E
  23 Species_C        1458                0      42          0  0.4890261                  0 0.028000000  0.00000000         Gene_E
  24 Species_D        1475                0      25          0  0.4894915                  0 0.016666667  0.00000000         Gene_E
  25 Species_E        1465                0      35          0  0.4771331                  0 0.023333333  0.00000000         Gene_E

  # Get species composition from dummy alignment
  getSpeciesComposition(alignment_path = rb_dummy_align_path)
  
         Taxon Total_Bases Total_Degenerate Total_N Total_Gaps Percent_GC Percent_Degenerate  Percent_N Percent_Gap     Alignment_Name
  1  Species_A          12                0       1          0 0.00000000         0.00000000 0.07692308  0.00000000 Dummy_Alignment.fa
  2  Species_B          12                0       1          0 0.00000000         0.00000000 0.07692308  0.00000000 Dummy_Alignment.fa
  3  Species_C          11                0       2          0 0.00000000         0.00000000 0.15384615  0.00000000 Dummy_Alignment.fa
  4  Species_D          11                0       2          0 0.00000000         0.00000000 0.15384615  0.00000000 Dummy_Alignment.fa
  5  Species_E          11                0       2          0 0.09090909         0.00000000 0.15384615  0.00000000 Dummy_Alignment.fa
  6  Species_F          10                1       2          0 0.10000000         0.07692308 0.15384615  0.00000000 Dummy_Alignment.fa
  7  Species_G           9                1       2          1 0.33333333         0.07692308 0.15384615  0.07692308 Dummy_Alignment.fa
  8  Species_H          10                1       2          0 0.40000000         0.07692308 0.15384615  0.00000000 Dummy_Alignment.fa
  9  Species_I           9                1       2          1 0.33333333         0.07692308 0.15384615  0.07692308 Dummy_Alignment.fa
  10 Species_J           7                1       2          3 0.42857143         0.07692308 0.15384615  0.23076923 Dummy_Alignment.fa

**Dummy Alignment**

.. image:: ../images/Dummy_Align.png
  :width: 600
