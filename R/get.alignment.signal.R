#' Rboretum Alignment Signal Fetcher
#'
#' Given the path to an alignment and list of taxa, this script returns site patterns for each site in the alignment (after pruning if needed)
#' @param alignment_path Path to alignment file (absolute or relative)
#' @param species_info Can be EITHER (1) A phylo object from which species will be extracted; or (2) a character vector of desired taxa (> 3 species required)
#' @param informative_gaps If TRUE, gaps in the alignment (-) are treated as potentially informative indels. If FALSE, gaps are considered missing data and not used (Default: FALSE)
#' @param alignment_name Name for alignment. If missing, the base filename is used
#' @return Dataframe containing split pattern for each site in the alignment, relative to the given set of taxa
#' @export
#' @examples
#' myAlignPath <- '/path/to/align.phy'
#'
#' mySpecies <- treeToCheck # Phylo object with species of interest
#' OR
#' mySpecies <- c('Spp1','Spp2','Spp3'...)
#'
#' # Gaps are indels
#' get.alignment.signal(myAlignPath,mySpecies,informative_gaps = TRUE)
#'
#' # Gaps are missing data
#' get.alignment.signal(myAlignPath,mySpecies,informative_gaps = FALSE)
#'
#' # Add alignment name
#' get.alignment.signal(myAlignPath,mySpecies,informative_gaps = FALSE, alignment_name = 'gene_XYZ')
#'
get.alignment.signal <- function(alignment_path,species_info,informative_gaps,alignment_name){

  # Set whether gaps are treated as missing data or indels
  if(missing(informative_gaps)){
    informative_gaps <- FALSE
  } else if(!is.logical(informative_gaps)){
    informative_gaps <- FALSE
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
  if(informative_gaps){
    info_gap <- '1'
  } else{ info_gap <- '0'}

  # Get list of species from tree or from character vector
  if(typeof(species_info)=="character" && length(species_info) > 3){
    spp_list = paste(sort(species_info),collapse = ";")
  } else if(!has_error(attributes(species_info)$class)){
      if(attributes(species_info)$class == "phylo"){
        spp_list = paste(sort(species_info$tip.label),collapse = ";")
      } else{ stop("Argument 3 (Species information) is neither a phylo object, nor a character vector of species longer than 3.")}
      } else{ stop("Argument 3 (Species information) is neither a phylo object, nor a character vector of species longer than 3.")} 

  split_support <- getSplitSupport(alignment_path,info_gap,spp_list)

  if(!is.data.frame(split_support)){
    stop("Error in parsing alingment. Check alignment to ensure chosen taxa are present and identically named.")
  }

  split_support <- split_support %>%
    mutate(Alignment_Name = alignment_name) %>% # Add alignment name
    rename(Alignment_Position = Zeroed_Site_Position) %>%
    mutate(Alignment_Position = as.integer(Alignment_Position) + 1) %>% # Convert position to 1-based
    select(Alignment_Name,Alignment_Position,Site_Pattern,Non_Base_Taxa,Non_Base_Count,Singleton_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5)

  return(split_support)
}
