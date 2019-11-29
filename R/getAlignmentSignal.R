#' Rboretum Alignment Signal Fetcher
#'
#' Given the path to an alignment and list of taxa, this script returns site patterns for each site in the alignment (after pruning if needed)
#' @param alignment_path Path to alignment file (absolute or relative)
#' @param species_info Set of 3+ species that can be passed as either:
#' \itemize{
#'   \item phylo object from which species will be extracted; or
#'   \item Character vector of desired taxa
#' }
#' @param use_gaps
#' \itemize{
#'   \item FALSE (Default: Treat all gaps (-) in alignment as missing data)
#'   \item TRUE (Treat all gaps in alignment as potentially informative signal)
#' }
#' @param alignment_name Name for alignment. If missing, the base filename is used
#' @return Dataframe containing split pattern for each site in the alignment, relative to the given set of taxa
#' @export

getAlignmentSignal <- function(alignment_path,species_info,use_gaps,alignment_name){

  # Set whether gaps are treated as missing data or indels
  if(missing(use_gaps)){
    use_gaps <- FALSE
  } else if(!is.logical(use_gaps)){
    use_gaps <- FALSE
  }
  
  if(has_error(file_path_as_absolute(alignment_path))){
    stop("'alignment_path' does not point to a valid file.")
  }

  if(missing(alignment_name)){
    # Set and check path to alignment file
    if(!file.exists(file_path_as_absolute(alignment_path))){
      stop("Argument 1 (Path to alignment) does not point to an existing file.")
    } else{
      alignment_path <- file_path_as_absolute(alignment_path)
      alignment_name <- file_path_sans_ext(basename(alignment_path))}
    } else{
      if(!file.exists(file_path_as_absolute(alignment_path))){
        stop("Argument 1 (Path to alignment) does not point to an existing file.")
    } else{
      alignment_path <- file_path_as_absolute(alignment_path)
      alignment_name <- as.character(alignment_name)}
    }

  # For Python passing...
  if(use_gaps){
    info_gap <- '1'
  } else{ info_gap <- '0'}

  # Get list of species from tree or from character vector
  if(Rboretum::isPhylo(species_info)){
    spp_list = paste(sort(species_info$tip.label),collapse = ";")
  } else if(typeof(species_info)=="character" && length(species_info) > 3){
    spp_list = paste(sort(species_info),collapse = ";")
  } else { stop("'species_info' is not a phylo object or character vector 3+ species IDs") }

  split_support <- getSplitSupport(alignment_path,info_gap,spp_list)

  if(!is.data.frame(split_support)){
    print(alignment_path)
    print("Error in parsing alignment above, returning empty dataframe. Check alignment to ensure chosen taxa are present and identically named.")
    return(tibble(Alignment_Name=character(),Alignment_Position=integer(),Site_Pattern=character(),Gap=logical(),Singleton=logical(),Singleton_Taxa=character(),Non_Base_Taxa=character(),Non_Base_Count=integer(),Split_1=character(),Split_2=character(),Split_3=character(),Split_4=character(),Split_5=character()))
  }

  split_support <- split_support %>%
    mutate(Alignment_Name = alignment_name) %>% # Add alignment name
    rename(Alignment_Position = Zeroed_Site_Position) %>%
    mutate(Alignment_Position = as.integer(Alignment_Position) + 1) %>% # Convert position to 1-based
    mutate(Gap = ifelse(Site_Pattern == 'pentallelic' | str_detect(Site_Pattern,'gap_'),TRUE,FALSE)) %>%
    mutate(Site_Pattern = str_replace(Site_Pattern,'gap_','')) %>%
    mutate(Singleton = ifelse(is.na(Singleton_Taxa),FALSE,TRUE)) %>%
    select(Alignment_Name,Alignment_Position,Site_Pattern,Gap,Singleton,Singleton_Taxa,Non_Base_Taxa,Non_Base_Count,Split_1,Split_2,Split_3,Split_4,Split_5)

  return(split_support)
}
