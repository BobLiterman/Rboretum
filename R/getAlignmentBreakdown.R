#' Rboretum Alignment Breaker Downer
#'
#' Given the path(s) to multiple alignments and an optional list of taxa, this script returns a full breakdown of alignment statistics (including %GC, %N, %Invariant, etc.)
#' @param alignment_path Where to find alignment files. Options include:
#' \itemize{
#'   \item A character vector of one or more alignment file paths  (relative or absolute)
#'   \item A path to a single directory containing all alignment files (relative or absolute)
#' }
#' @param species_info OPTIONAL: List of taxa to analyze [Default: Process alignment using all taxa if solo, all shared in multiple]. Can be provided as:
#' \itemize{
#'   \item phylo object from which species will be extracted; or
#'   \item multiPhylo object from which common species will be extracted; or
#'   \item Character vector with 3+ taxon IDs
#'   \item Semicolon-separated list of taxon IDs
#' }
#' @param use_gaps OPTIONAL: If FALSE, treat gaps as missing data (like N) [Default: TRUE, treat gaps as indel data]
#' @param alignment_name OPTIONAL: Chacter vector of names for each alignment. If missing or incomplete, the base filename is used
#' @param prefix OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file prefixes (e.g. all alignment files start with "Mafft_")
#' @param suffix OPTIONAL: If 'alignment_path' is a directory, provide a character vector of file suffixes (e.g. all alignment files end with ".phy")
#' @return Dataframe containing a breakdown of alignment statistics
#' @export

getAlignmentBreakdown <- function(alignment_path,species_info,use_gaps,alignment_name,prefix,suffix){
  
  # Ensure that a path and root taxa are provided as character vectors
  if(missing(alignment_path)){
    stop("No alignment file or directories indicated with 'alignment_path'")
  } else if(!is.character(alignment_path)){
    stop("'alignment_path' should be a character vector of file paths or the path to an alignment directory.")
  }
  
  # Set missing values
  if(missing(species_info)){
    species_info <- substitute()
  }
  
  if(missing(use_gaps)){
    use_gaps <- substitute()
  }
  
  if(missing(alignment_name)){
    alignment_name <- substitute()
  }
  
  if(missing(prefix)){
    prefix <- substitute()
  }
  
  if(missing(suffix)){
    suffix <- substitute()  
  }
  
  # Get alignment composition data
  align_comp_df <- getAlignmentComposition(alignment_path,species_info,alignment_name,prefix,suffix)
  
  # Get species composition data
  species_comp_df <- getSpeciesComposition(alignment_path,species_info,alignment_name,suffix)
  
  # Get alignment patterns
  pattern_df <- getAlignmentPatterns(alignment_path,species_info,use_gaps,alignment_name,prefix,suffix)

  # Create breakdown_df
  breakdown_df <- align_comp_df %>% 
    rename(Alignment_Percent_GC = 'Percent_GC')
  
  alignment_names <- unique(breakdown_df$Alignment_Name)
  
  # Add GC Skew
  gc_skew <- species_comp_df %>% filter(Total_Bases>0) %>% # Ensure samples have at least one base
    group_by(Alignment_Name) %>%
    summarize(Species_GC_Mean = mean(Percent_GC),Species_GC_StdDev = sd(Percent_GC)) %>% ungroup()
  
  breakdown_df <- left_join(breakdown_df,gc_skew,by='Alignment_Name')
  
  # Get pattern percentanges (Non-Base, Invariant, Singleton, Parsimony-Informative, Biallelic, Triallelic, Quadallelic, Pentallelic)
  pattern_breakdown_df <- tibble(Alignment_Name=character(),Percent_Nonbase=numeric(),Percent_Invariant=numeric(),
                                 Percent_Singleton=numeric(),Percent_Parsimony_Informative=numeric(),Percent_Biallelic=numeric(),
                                 Percent_Triallelic=numeric(),Percent_Quadallelic=numeric(),Percent_Pentallelic=numeric())
  
  for(align in alignment_names){
    
    temp_df <- pattern_df %>%
      filter(Alignment_Name == align)
    
    # If alignment contains data...
    if(!nrow(temp_df)==0){
      
      alignment_length <- nrow(temp_df)
      
      
      pattern_table <- pull(temp_df,Site_Pattern) %>% table()
      
      pattern_breakdown_df <- pattern_breakdown_df %>% add_row(Alignment_Name = align,
                                                               Percent_Nonbase = tableCount(pattern_table,'non_base')/alignment_length,
                                                               Percent_Invariant = tableCount(pattern_table,'invariant')/alignment_length,
                                                               Percent_Singleton = tableCount(pattern_table,'singleton')/alignment_length,
                                                               Percent_Parsimony_Informative = nrow(temp_df[temp_df$Parsimony_Informative=="Yes",])/alignment_length,
                                                               Percent_Biallelic = tableCount(pattern_table,'biallelic')/alignment_length,
                                                               Percent_Triallelic = tableCount(pattern_table,'triallelic')/alignment_length,
                                                               Percent_Quadallelic = tableCount(pattern_table,'quadallelic')/alignment_length,
                                                               Percent_Pentallelic = tableCount(pattern_table,'pentallelic')/alignment_length)
    }
  }
  
  breakdown_df <- left_join(breakdown_df,pattern_breakdown_df,by='Alignment_Name')
  
  return(breakdown_df)
}