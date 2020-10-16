.. _getAlignmentSignal:

########################
**getAlignmentSignal**
########################

*getAlignmentSignal* takes one or more alignment files and returns site-by-site information about:

  - Site variation patterns (i.e. is the site invariant, singleton, biallelic, triallelic, etc.)
  - Indels (if gap data is not treated as missing)
  - Missing data
  - Parsimony-based splits (i.e. which taxonomic groupings are supported by alleles at each site)
  
The output from *getAlignmentSignal* can be used to extract support values for different relationships via **getTreeSupport**

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

  # Set test data directory
  sourceRboretum()

  # Set path to alignment data
  myAlignmentFile <- rb_align1_path
  mySpecies <- getAlignmentSpecies(myAlignmentFile)
  myAlignmentDir <- rb_alignment_dir
  
  # Get alignment signal information for a single alignment
  getAlignmentSignal(alignment_path = myAlignmentFile)
  
  # A tibble: 2,000 x 14
     Alignment_Positi~ All_Site_Bases Site_Pattern Non_Base_Taxa Non_Base_Count Singleton_Taxa Singleton_Count Gap_Taxa Split_1              Split_2                    Split_3             Split_4         Split_5 Alignment_Name
                 <dbl> <chr>          <chr>                <dbl>          <dbl> <chr>                    <dbl>    <dbl> <chr>                <chr>                      <chr>               <chr>             <dbl> <chr>         
   1                 1 a;c;t          triallelic              NA              0 NA                           0       NA Species_A;Species_B~ Species_J;Species_K;Speci~ Species_C;Species_~ NA                   NA Gene_1.phylip 
   2                 2 a;c;t          triallelic              NA              0 NA                           0       NA Species_A;Species_B  Species_C;Species_D;Speci~ Species_G;Species_J NA                   NA Gene_1.phylip 
   3                 3 a;c            biallelic               NA              0 NA                           0       NA Species_G;Species_N~ Species_A;Species_B;Speci~ NA                  NA                   NA Gene_1.phylip 
   4                 4 a;c;g;t        triallelic              NA              0 Species_J                    1       NA Species_A;Species_B~ Species_C;Species_D;Speci~ Species_H;Species_~ NA                   NA Gene_1.phylip 
   5                 5 a;c;g;t        triallelic              NA              0 Species_L                    1       NA Species_M;Species_N~ Species_A;Species_B;Speci~ Species_C;Species_~ NA                   NA Gene_1.phylip 
   6                 6 a;c;t          triallelic              NA              0 NA                           0       NA Species_C;Species_D~ Species_J;Species_M;Speci~ Species_A;Species_B NA                   NA Gene_1.phylip 
   7                 7 a;c;g          biallelic               NA              0 Species_J                    1       NA Species_A;Species_B~ Species_C;Species_D;Speci~ NA                  NA                   NA Gene_1.phylip 
   8                 8 a;g            biallelic               NA              0 NA                           0       NA Species_H;Species_I  Species_A;Species_B;Speci~ NA                  NA                   NA Gene_1.phylip 
   9                 9 a;c;g;t        quadallelic             NA              0 NA                           0       NA Species_G;Species_H~ Species_M;Species_N;Speci~ Species_E;Species_F Species_A;Spec~      NA Gene_1.phylip 
  10                10 a;g;t          triallelic              NA              0 NA                           0       NA Species_C;Species_D~ Species_H;Species_I;Speci~ Species_A;Species_B NA                   NA Gene_1.phylip 
  # ... with 1,990 more rows
  
  
  # Get alignment signal information from all .phylip files in a directory, providing new names, consider gaps as missing data
  getAlignmentSignal(alignment_path = myAlignmentDir,species_info = mySpecies,use_gaps = FALSE,suffix = ".phylip",alignment_name = c('Gene_A','Gene_B','Gene_C','Gene_D','Gene_E'))

  # A tibble: 11,500 x 14
     Alignment_Positi~ All_Site_Bases Site_Pattern Non_Base_Taxa Non_Base_Count Singleton_Taxa Singleton_Count Gap_Taxa Split_1              Split_2                    Split_3             Split_4         Split_5 Alignment_Name
                 <dbl> <chr>          <chr>                <dbl>          <dbl> <chr>                    <dbl>    <dbl> <chr>                <chr>                      <chr>               <chr>             <dbl> <chr>         
   1                 1 a;c;t          triallelic              NA              0 NA                           0       NA Species_A;Species_B~ Species_J;Species_K;Speci~ Species_C;Species_~ NA                   NA Gene_A        
   2                 2 a;c;t          triallelic              NA              0 NA                           0       NA Species_A;Species_B  Species_C;Species_D;Speci~ Species_G;Species_J NA                   NA Gene_A        
   3                 3 a;c            biallelic               NA              0 NA                           0       NA Species_G;Species_N~ Species_A;Species_B;Speci~ NA                  NA                   NA Gene_A        
   4                 4 a;c;g;t        triallelic              NA              0 Species_J                    1       NA Species_A;Species_B~ Species_C;Species_D;Speci~ Species_H;Species_~ NA                   NA Gene_A        
   5                 5 a;c;g;t        triallelic              NA              0 Species_L                    1       NA Species_M;Species_N~ Species_A;Species_B;Speci~ Species_C;Species_~ NA                   NA Gene_A        
   6                 6 a;c;t          triallelic              NA              0 NA                           0       NA Species_C;Species_D~ Species_J;Species_M;Speci~ Species_A;Species_B NA                   NA Gene_A        
   7                 7 a;c;g          biallelic               NA              0 Species_J                    1       NA Species_A;Species_B~ Species_C;Species_D;Speci~ NA                  NA                   NA Gene_A        
   8                 8 a;g            biallelic               NA              0 NA                           0       NA Species_H;Species_I  Species_A;Species_B;Speci~ NA                  NA                   NA Gene_A        
   9                 9 a;c;g;t        quadallelic             NA              0 NA                           0       NA Species_G;Species_H~ Species_M;Species_N;Speci~ Species_E;Species_F Species_A;Spec~      NA Gene_A        
  10                10 a;g;t          triallelic              NA              0 NA                           0       NA Species_C;Species_D~ Species_H;Species_I;Speci~ Species_A;Species_B NA                   NA Gene_A        
  # ... with 11,490 more rows

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
