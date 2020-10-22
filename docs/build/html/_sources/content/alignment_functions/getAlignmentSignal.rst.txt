.. _getAlignmentSignal:

########################
**getAlignmentSignal**
########################

**getAlignmentSignal** takes one or more alignment files and returns site-by-site information about:

  - Site variation patterns. Possible site patterns include:
  
    - **non_base**: Sites where less than three species have a called, fixed base (A/C/T/G/-)
    - **invariant**: Sites made up of a single allele (and possibly missing data)
    - **singleton**: Otherwise invariant sites where one or more individual taxa have a unique, separate A/T/C/G/- allele
    - **biallelic**: Sites with two possible alleles + possible missing or singleton taxa
    - **triallelic**: Sites with three possible alleles + possible missing or singleton taxa
    - **quadallelic**: Sites with four possible alleles + possible missing or singleton taxa
    - **pentallelic**: Sites with five possible alleles (A + C + T + G + '-') + possible missing or singleton taxa
    
  - Which taxa have indels (if gap data is not treated as missing)
  - Which taxa have missing data
  - Parsimony-based splits (i.e. which taxonomic groupings are supported by alleles at each site)
  
The output from **getAlignmentSignal** can be used to extract support values for different potential relationships via **getAlignmentSupport**. 

.. warning::
  
  **Note:** Rboretum currently does not support handling of degenerate bases, and classifies them broadly as N/missing. Support for degenerate nucleotides is in progress and set for a future release. 

=======================
Function and Arguments
=======================

**Usage**:
::

  getAlignmentSignal <- function(alignment_path,species_info,use_gaps,alignment_name,prefix,suffix,existing_signal){

===========================      ===============================================================================================================================================================================================================
 Argument                         Description
===========================      ===============================================================================================================================================================================================================
**alignment_path**				        An absolute or relative path to an alignment file(s) or a directory containing multiple alignments
**species_info**                  A list of species to consider, provided as a phylo, multiPhlyo, character vector or labels or a semicolon-delimited string [Default: Shared species among all alignments]
**use_gaps**                      If FALSE, consider gaps (-) in alignments as missing data. [Default: TRUE, treat gaps as indel characters]
**alignment_name**                A character vector containing desired alignment IDs [Default: Derive name from filename]
**prefix**                        If **alignment_path** points to a directory, only query files starting with **prefix** (e.g. 'Alignment') [Default: Use all files in directory]
**suffix**                        If **alignment_path** points to a directory, only query files ending wtih **suffix** (e.g. '.nex') [Default: Use all files in directory]
**existing_signal**               OPTIONAL: Append these results to the output from getAlignmentSignal() run with the same species_info and different alignment(s)
===========================      ===============================================================================================================================================================================================================

==============
Example Usage
==============

.. code-block:: r
  
  # Script: Rboretum/docs/content/Doc_Function_Scripts/getAlignmentSignal.R

  library(Rboretum)
  sourceRboretum()
  
  # Set path to alignment data
  myAlignmentFile <- rb_align1_path
  mySpecies <- getAlignmentSpecies(myAlignmentFile)
  mySubspecies <- semiVector(mySpecies)[1:5]
  myAlignmentDir <- rb_alignment_dir
  
  # Get alignment signal information for a single alignment
  getAlignmentSignal(alignment_path = myAlignmentFile)

  # A tibble: 1,551 x 14
     Alignment_Positi~ All_Site_Bases Site_Pattern Non_Base_Taxa Non_Base_Count Singleton_Taxa Singleton_Count Gap_Taxa Split_1                           Split_2                           Split_3 Split_4 Split_5 Alignment_Name
                 <dbl> <chr>          <chr>        <chr>                  <dbl> <chr>                    <dbl> <chr>    <chr>                             <chr>                             <chr>     <dbl>   <dbl> <chr>         
   1                 1 C              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   2                 2 C              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   3                 3 C;N;T          biallelic    Species_D                  1 NA                           0 NA       Species_A;Species_F               Species_B;Species_C;Species_E;Sp~ NA           NA      NA Gene_1.phy    
   4                 4 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   5                 5 C;T            biallelic    NA                         0 NA                           0 NA       Species_C;Species_H               Species_A;Species_B;Species_D;Sp~ NA           NA      NA Gene_1.phy    
   6                 6 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   7                 7 C;G            singleton    NA                         0 Species_H                    1 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   8                 8 A;G            biallelic    NA                         0 NA                           0 NA       Species_A;Species_B;Species_D;Sp~ Species_C;Species_H               NA           NA      NA Gene_1.phy    
   9                 9 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
  10                10 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
  # ... with 1,541 more rows
  
  # Get alignment signal information from all .phylip files in a directory, providing new names, consider gaps as missing data
  getAlignmentSignal(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

  # A tibble: 1,551 x 14
     Alignment_Positi~ All_Site_Bases Site_Pattern Non_Base_Taxa Non_Base_Count Singleton_Taxa Singleton_Count Gap_Taxa Split_1                           Split_2                           Split_3 Split_4 Split_5 Alignment_Name
                 <dbl> <chr>          <chr>        <chr>                  <dbl> <chr>                    <dbl> <chr>    <chr>                             <chr>                             <chr>     <dbl>   <dbl> <chr>         
   1                 1 C              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   2                 2 C              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   3                 3 C;N;T          biallelic    Species_D                  1 NA                           0 NA       Species_A;Species_F               Species_B;Species_C;Species_E;Sp~ NA           NA      NA Gene_1.phy    
   4                 4 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   5                 5 C;T            biallelic    NA                         0 NA                           0 NA       Species_C;Species_H               Species_A;Species_B;Species_D;Sp~ NA           NA      NA Gene_1.phy    
   6                 6 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   7                 7 C;G            singleton    NA                         0 Species_H                    1 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
   8                 8 A;G            biallelic    NA                         0 NA                           0 NA       Species_A;Species_B;Species_D;Sp~ Species_C;Species_H               NA           NA      NA Gene_1.phy    
   9                 9 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
  10                10 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA           NA      NA Gene_1.phy    
  # ... with 1,541 more rows
  
  # Get alignment signal information from all .phylip files in a directory, providing new names, consider gaps as missing data
  getAlignmentSignal(alignment_path = myAlignmentDir,species_info = mySubspecies,use_gaps = FALSE,suffix = ".phy",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

  # A tibble: 9,105 x 14
     Alignment_Positi~ All_Site_Bases Site_Pattern Non_Base_Taxa Non_Base_Count Singleton_Taxa Singleton_Count Gap_Taxa Split_1                           Split_2                           Split_3 Split_4 Split_5 Alignment_Name
                 <dbl> <chr>          <chr>        <chr>                  <dbl> <chr>                    <dbl> <chr>    <chr>                             <chr>                             <chr>   <chr>     <dbl> <chr>         
   1                 1 C              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA      NA           NA Gene_A        
   2                 2 C              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA      NA           NA Gene_A        
   3                 3 C;N;T          biallelic    Species_D                  1 NA                           0 NA       Species_A;Species_F               Species_B;Species_C;Species_E;Sp~ NA      NA           NA Gene_A        
   4                 4 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA      NA           NA Gene_A        
   5                 5 C;T            biallelic    NA                         0 NA                           0 NA       Species_C;Species_H               Species_A;Species_B;Species_D;Sp~ NA      NA           NA Gene_A        
   6                 6 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA      NA           NA Gene_A        
   7                 7 C;G            singleton    NA                         0 Species_H                    1 NA       NA                                NA                                NA      NA           NA Gene_A        
   8                 8 A;G            biallelic    NA                         0 NA                           0 NA       Species_A;Species_B;Species_D;Sp~ Species_C;Species_H               NA      NA           NA Gene_A        
   9                 9 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA      NA           NA Gene_A        
  10                10 T              invariant    NA                         0 NA                           0 NA       NA                                NA                                NA      NA           NA Gene_A        
  # ... with 9,095 more rows
  
  # Get alignment signal from dummy alignment, with and without gap support
  getAlignmentSignal(alignment_path = rb_dummy_align_path)

  # A tibble: 13 x 14
     Alignment_Positi~ All_Site_Bases Site_Pattern Non_Base_Taxa                   Non_Base_Count Singleton_Taxa      Singleton_Count Gap_Taxa    Split_1         Split_2         Split_3     Split_4     Split_5   Alignment_Name
                 <dbl> <chr>          <chr>        <chr>                                    <dbl> <chr>                         <dbl> <chr>       <chr>           <chr>           <chr>       <chr>       <chr>     <chr>         
   1                 1 N              non_base     Species_A;Species_B;Species_C;~             10 NA                               NA NA          NA              NA              NA          NA          NA        Dummy_Alignme~
   2                 2 A;N            non_base     Species_C;Species_D;Species_E;~              8 NA                               NA NA          NA              NA              NA          NA          NA        Dummy_Alignme~
   3                 3 A              invariant    NA                                           0 NA                                0 NA          NA              NA              NA          NA          NA        Dummy_Alignme~
   4                 4 A;T            singleton    NA                                           0 Species_J                         1 NA          NA              NA              NA          NA          NA        Dummy_Alignme~
   5                 5 -;A            singleton    NA                                           0 Species_J                         1 Species_J   NA              NA              NA          NA          NA        Dummy_Alignme~
   6                 6 -;A;C;G;T      singleton    NA                                           0 Species_G;Species_~               4 Species_G   NA              NA              NA          NA          NA        Dummy_Alignme~
   7                 7 A;T            biallelic    NA                                           0 NA                                0 NA          Species_A;Spec~ Species_F;Spec~ NA          NA          NA        Dummy_Alignme~
   8                 8 A;C;T          triallelic   NA                                           0 NA                                0 NA          Species_A;Spec~ Species_G;Spec~ Species_D;~ NA          NA        Dummy_Alignme~
   9                 9 A;C;G;T        quadallelic  NA                                           0 NA                                0 NA          Species_A;Spec~ Species_G;Spec~ Species_I;~ Species_D;~ NA        Dummy_Alignme~
  10                10 -;A;C;G;T      pentallelic  NA                                           0 NA                                0 Species_I;~ Species_I;Spec~ Species_A;Spec~ Species_E;~ Species_G;~ Species_~ Dummy_Alignme~
  11                11 -;A;T          biallelic    NA                                           0 Species_J                         1 Species_J   Species_A;Spec~ Species_F;Spec~ NA          NA          NA        Dummy_Alignme~
  12                12 A;G;T          biallelic    NA                                           0 Species_J                         1 NA          Species_A;Spec~ Species_F;Spec~ NA          NA          NA        Dummy_Alignme~
  13                13 A;K;M;R;S;W    invariant    Species_F;Species_G;Species_H;~              5 NA                                0 NA          NA              NA              NA          NA          NA        Dummy_Alignme~
  
  getAlignmentSignal(alignment_path = rb_dummy_align_path,use_gaps = FALSE)

  # A tibble: 13 x 14
     Alignment_Positi~ All_Site_Bases Site_Pattern Non_Base_Taxa                    Non_Base_Count Singleton_Taxa     Singleton_Count Gap_Taxa    Split_1          Split_2          Split_3     Split_4     Split_5 Alignment_Name
                 <dbl> <chr>          <chr>        <chr>                                     <dbl> <chr>                        <dbl> <chr>       <chr>            <chr>            <chr>       <chr>         <dbl> <chr>         
   1                 1 N              non_base     Species_A;Species_B;Species_C;S~             10 NA                              NA NA          NA               NA               NA          NA               NA Dummy_Alignme~
   2                 2 A;N            non_base     Species_C;Species_D;Species_E;S~              8 NA                              NA NA          NA               NA               NA          NA               NA Dummy_Alignme~
   3                 3 A              invariant    NA                                            0 NA                               0 NA          NA               NA               NA          NA               NA Dummy_Alignme~
   4                 4 A;T            singleton    NA                                            0 Species_J                        1 NA          NA               NA               NA          NA               NA Dummy_Alignme~
   5                 5 -;A            invariant    Species_J                                     1 NA                               0 Species_J   NA               NA               NA          NA               NA Dummy_Alignme~
   6                 6 -;A;C;G;T      singleton    Species_G                                     1 Species_H;Species~               3 Species_G   NA               NA               NA          NA               NA Dummy_Alignme~
   7                 7 A;T            biallelic    NA                                            0 NA                               0 NA          Species_A;Speci~ Species_F;Speci~ NA          NA               NA Dummy_Alignme~
   8                 8 A;C;T          triallelic   NA                                            0 NA                               0 NA          Species_A;Speci~ Species_G;Speci~ Species_D;~ NA               NA Dummy_Alignme~
   9                 9 A;C;G;T        quadallelic  NA                                            0 NA                               0 NA          Species_A;Speci~ Species_G;Speci~ Species_I;~ Species_D;~      NA Dummy_Alignme~
  10                10 -;A;C;G;T      quadallelic  Species_I;Species_J                           2 NA                               0 Species_I;~ Species_A;Speci~ Species_E;Speci~ Species_G;~ Species_C;~      NA Dummy_Alignme~
  11                11 -;A;T          biallelic    Species_J                                     1 NA                               0 Species_J   Species_A;Speci~ Species_F;Speci~ NA          NA               NA Dummy_Alignme~
  12                12 A;G;T          biallelic    NA                                            0 Species_J                        1 NA          Species_A;Speci~ Species_F;Speci~ NA          NA               NA Dummy_Alignme~
  13                13 A;K;M;R;S;W    invariant    Species_F;Species_G;Species_H;S~              5 NA                               0 NA          NA               NA               NA          NA               NA Dummy_Alignme~

  # Postion 2 is 'non_base' because < 3 species have a called base
  # Note: Sites 5, 6, 10, and 11 have species with gap positions. 
  # Treating gaps as missing data sets all gap taxa to missing taxa in the bottom dataframe, and also changes the reported site patterns for rows 5 + 10

**Dummy Alignment**
  
.. image:: ../images/Dummy_Align.png
  :width: 600
