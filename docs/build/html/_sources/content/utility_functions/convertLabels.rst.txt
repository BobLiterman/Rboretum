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

  rb_name_file
  [1] "<PACKAGE_DIR>/inst/extdata/Name_Conversion_Table.tsv"

  # Read in table of name equivalencies

  name_information <- read_tsv(rb_name_file)
  head(name_information)

  # A tibble: 6 x 3
  Alignment_IDs Number_IDs Common_Names         
  <chr>         <chr>      <chr>                
  1 Species_A     Species_1  Fuzzy Wollop         
  2 Species_B     Species_2  Green-footed Squeazel
  3 Species_C     Species_3  Lefty Lucy           
  4 Species_D     Species_4  Cheeto Bandito       
  5 Species_E     Species_5  Louse-by-the-Sea     
  6 Species_F     Species_6  Damp Trickster       

  # Read in tree with 'Alignment_IDs' as tip labels

  alignmentTree <- readRooted(rb_tree1_path,root_taxa = c('Species_C','Species_H'))
  alignmentTree$tip.label

  [1] "Species_H" "Species_C" "Species_K" "Species_L" "Species_A" "Species_F" "Species_E" "Species_M" "Species_J" "Species_I"
  [11] "Species_N" "Species_G" "Species_B" "Species_O" "Species_D"

  # Convert 'Alignment_IDs' [First column in name table] to 'Number_IDs' [Second column in name table]

  numberTree <- convertLabels(to_convert = alignmentTree,name_df = name_information)
  numberTree$tip.label

  [1] "Species_8"  "Species_3"  "Species_11" "Species_12" "Species_1"  "Species_6"  "Species_5"  "Species_13" "Species_10"
  [10] "Species_9"  "Species_14" "Species_7"  "Species_2"  "Species_15" "Species_4" 

  # Convert 'Number_IDs' to 'Common_Names'

  commonTree <- convertLabels(to_convert = numberTree,name_df = name_information,from='Number_IDs',to='Common_Names')
  commonTree$tip.label

  [1] "Clammy Warbeetle"      "Lefty Lucy"            "Steve's Fwonkton"      "Fat Jellysnorp"        "Fuzzy Wollop"         
  [6] "Damp Trickster"        "Louse-by-the-Sea"      "Three-clawed Sampson"  "Deborah's Downer"      "Six-legged Snake"     
  [11] "Matzo Monster"         "Wonky Algo"            "Green-footed Squeazel" "Languishing Trunkle"   "Cheeto Bandito"       

  # Convert a set of clades

  alignment_clades <- getTreeClades(alignmentTree)
  alignment_clades

  [1] "Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N;Species_O"
  [2] "Species_A;Species_E;Species_F;Species_G;Species_I;Species_J;Species_K;Species_L;Species_M;Species_N"                              
  [3] "Species_A;Species_E;Species_F;Species_K;Species_L"                                                                                
  [4] "Species_A;Species_F"                                                                                                              
  [5] "Species_A;Species_F;Species_K;Species_L"                                                                                          
  [6] "Species_B;Species_D;Species_O"                                                                                                    
  [7] "Species_B;Species_O"                                                                                                              
  [8] "Species_C;Species_H"                                                                                                              
  [9] "Species_G;Species_I;Species_J;Species_M;Species_N"                                                                                
  [10] "Species_G;Species_I;Species_N"                                                                                                    
  [11] "Species_G;Species_N"                                                                                                              
  [12] "Species_J;Species_M"                                                                                                              
  [13] "Species_K;Species_L"                                                                                                              

  number_clades <- convertLabels(to_convert = alignment_clades,name_df = name_information,from='Alignment_IDs',to='Number_IDs')
  number_clades

  [1] "Species_1;Species_10;Species_11;Species_12;Species_13;Species_14;Species_15;Species_2;Species_4;Species_5;Species_6;Species_7;Species_9"
  [2] "Species_1;Species_10;Species_11;Species_12;Species_13;Species_14;Species_5;Species_6;Species_7;Species_9"                               
  [3] "Species_1;Species_11;Species_12;Species_5;Species_6"                                                                                    
  [4] "Species_1;Species_6"                                                                                                                    
  [5] "Species_1;Species_11;Species_12;Species_6"                                                                                              
  [6] "Species_15;Species_2;Species_4"                                                                                                         
  [7] "Species_15;Species_2"                                                                                                                   
  [8] "Species_3;Species_8"                                                                                                                    
  [9] "Species_10;Species_13;Species_14;Species_7;Species_9"                                                                                   
  [10] "Species_14;Species_7;Species_9"                                                                                                         
  [11] "Species_14;Species_7"                                                                                                                   
  [12] "Species_10;Species_13"                                                                                                                  
  [13] "Species_11;Species_12"                                                                                                                  

  # Convert columns in a dataframe

  alignment_splits <- getTreeSplits(alignmentTree)
  alignment_splits

  # A tibble: 12 x 4
  Clade                                           Mirror_Clade                                                     Split_Node Root 
  <chr>                                           <chr>                                                                 <int> <lgl>
  1 Species_A;Species_E;Species_F;Species_G;Specie~ Species_B;Species_C;Species_D;Species_H;Species_O                        18 FALSE
  2 Species_A;Species_E;Species_F;Species_K;Specie~ Species_B;Species_C;Species_D;Species_G;Species_H;Species_I;Spe~         19 FALSE
  3 Species_A;Species_F;Species_K;Species_L         Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Spe~         20 FALSE
  4 Species_K;Species_L                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Spe~         21 FALSE
  5 Species_A;Species_F                             Species_B;Species_C;Species_D;Species_E;Species_G;Species_H;Spe~         22 FALSE
  6 Species_G;Species_I;Species_J;Species_M;Specie~ Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Spe~         23 FALSE
  7 Species_J;Species_M                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Spe~         24 FALSE
  8 Species_G;Species_I;Species_N                   Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Spe~         25 FALSE
  9 Species_G;Species_N                             Species_A;Species_B;Species_C;Species_D;Species_E;Species_F;Spe~         26 FALSE
  10 Species_B;Species_D;Species_O                   Species_A;Species_C;Species_E;Species_F;Species_G;Species_H;Spe~         27 FALSE
  11 Species_B;Species_O                             Species_A;Species_C;Species_D;Species_E;Species_F;Species_G;Spe~         28 FALSE
  12 Species_C;Species_H                             Species_A;Species_B;Species_D;Species_E;Species_F;Species_G;Spe~         16 TRUE 

  number_splits <- alignment_splits %>%
  +   rowwise() %>%
  +   mutate(Clade=convertLabels(Clade,name_information),
  +          Mirror_Clade=convertLabels(Mirror_Clade,name_information))

  number_splits

  # A tibble: 12 x 4
  # Rowwise: 
  Clade                                             Mirror_Clade                                                   Split_Node Root 
  <chr>                                             <chr>                                                               <int> <lgl>
  1 Species_1;Species_10;Species_11;Species_12;Speci~ Species_15;Species_2;Species_3;Species_4;Species_8                     18 FALSE
  2 Species_1;Species_11;Species_12;Species_5;Specie~ Species_10;Species_13;Species_14;Species_15;Species_2;Species~         19 FALSE
  3 Species_1;Species_11;Species_12;Species_6         Species_10;Species_13;Species_14;Species_15;Species_2;Species~         20 FALSE
  4 Species_11;Species_12                             Species_1;Species_10;Species_13;Species_14;Species_15;Species~         21 FALSE
  5 Species_1;Species_6                               Species_10;Species_11;Species_12;Species_13;Species_14;Specie~         22 FALSE
  6 Species_10;Species_13;Species_14;Species_7;Speci~ Species_1;Species_11;Species_12;Species_15;Species_2;Species_~         23 FALSE
  7 Species_10;Species_13                             Species_1;Species_11;Species_12;Species_14;Species_15;Species~         24 FALSE
  8 Species_14;Species_7;Species_9                    Species_1;Species_10;Species_11;Species_12;Species_13;Species~         25 FALSE
  9 Species_14;Species_7                              Species_1;Species_10;Species_11;Species_12;Species_13;Species~         26 FALSE
  10 Species_15;Species_2;Species_4                    Species_1;Species_10;Species_11;Species_12;Species_13;Species~         27 FALSE
  11 Species_15;Species_2                              Species_1;Species_10;Species_11;Species_12;Species_13;Species~         28 FALSE
  12 Species_3;Species_8                               Species_1;Species_10;Species_11;Species_12;Species_13;Species~         16 TRUE

  # Convert a list of IDs

  root_taxa <- alignment_splits %>% filter(Root) %>% pull(Clade) %>% semiVector()
  nonroot_taxa <- alignment_splits %>% filter(!Root) %>% pull(Clade) %>% semiVector() %>% unlist() %>% unique()

  taxa_list <- list('Root'=root_taxa,'Non_Root'=nonroot_taxa)
  taxa_list

  $Root
  [1] "Species_C" "Species_H"

  $Non_Root
  [1] "Species_A" "Species_E" "Species_F" "Species_G" "Species_I" "Species_J" "Species_K" "Species_L" "Species_M" "Species_N"
  [11] "Species_B" "Species_D" "Species_O"

  number_list <- convertLabels(to_convert = taxa_list ,name_df = name_information,from = 'Alignment_IDs',to='Number_IDs')
  number_list

  $Root
  [1] "Species_3" "Species_8"

  $Non_Root
  [1] "Species_1"  "Species_5"  "Species_6"  "Species_7"  "Species_9"  "Species_10" "Species_11" "Species_12" "Species_13"
  [10] "Species_14" "Species_2"  "Species_4"  "Species_15"
  
