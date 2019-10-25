#' Get Topology Support Counts from Alignment Data
#'
#' Given the path to an alignment and list of taxa, this script returns site patterns for each site in the alignment (after pruning if needed)
#' @param alignment_path Path to alignment file (absolute or relative)
#' @param species_info Can be EITHER (1) A phylo object from which species will be extracted; or (2) a vector of desired taxa (> 3 species required)
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
#' getAlignmentSignal(myAlignPath,mySpecies,informative_gaps = TRUE)
#'
#' # Gaps are missing data
#' getAlignmentSignal(myAlignPath,mySpecies,informative_gaps = FALSE)
#'
#' # Add alignment name
#' getAlignmentSignal(myAlignPath,mySpecies,informative_gaps = FALSE, alignment_name = 'gene_XYZ')
#'
getAlignmentSignal <- function(alignment_path,species_info,informative_gaps,alignment_name){

  # Set whether gaps are treated as missing data or indels
  if(missing(informative_gaps)){
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

  if(informative_gaps != TRUE){
    informative_gaps <- FALSE
  }

  # For Python passing...
  if(informative_gaps){
    info_gap <- '1'
  } else{ info_gap <- '0'}

  # Get list of species from tree or from character vector
  if(typeof(species_info)=="character" && length(species_info) > 3){
    spp_list = paste(sort(species_info),collapse = ";")
  } else{
    if(attributes(species_info)$class == "phylo"){
      spp_list = paste(sort(species_info$tip.label),collapse = ";")
    } else{ stop("Argument 3 (Species information) is neither a phylo object, nor a character vector of species longer than 3.") }
  }

  split_support <- getSplitSupport(alignment_path,info_gap,spp_list)

  if(!is.data.frame(split_support)){
    stop("Error in parsing alingment. Check alignment to ensure chosen taxa are present and identically named.")
  }

  site_count <- nrow(split_support)
  non_base_count <- split_support %>% filter(Site_Pattern=='non_base') %>% nrow()
  invariant_count <- split_support %>% filter(Site_Pattern %in% c('invariant','gap_invariant')) %>% nrow()
  singleton_count <- split_support %>% filter(Site_Pattern %in% c('singleton','gap_singleton')) %>% nrow()
  biallelic_count <- split_support %>% filter(Site_Pattern %in% c('biallelic','gap_biallelic')) %>% nrow()
  triallelic_count <- split_support %>% filter(Site_Pattern %in% c('triallelic','gap_triallelic')) %>% nrow()
  quadallelic_count <- split_support %>% filter(Site_Pattern %in% c('quadallelic','gap_quadallelic')) %>% nrow()
  pentallelic_count <- split_support %>% filter(Site_Pattern %in% c('pentallelic')) %>% nrow()

  # Print output information
  print(paste("Alignment File:",alignment_path))
  print(paste("Species Requested:",spp_list))

  print(paste("Alignment contains",site_count,"sites."))

  if(info_gap == "0"){
    print(paste("Alignment contains",non_base_count,"non-base sites. (Note: Gaps were interpreted as missing/non-base data.)"))
  }
  if(info_gap == "1"){
    print(paste("Alignment contains",non_base_count,"non-base sites. (Note: Gaps were interpreted as informative indels.)"))
  }

  print(paste("Alignment contains",invariant_count,"invariant sites."))
  print(paste("Alignment contains",singleton_count,"singleton sites."))
  print(paste("Alignment contains",biallelic_count,"biallelic sites."))
  print(paste("Alignment contains",triallelic_count,"triallelic sites."))
  print(paste("Alignment contains",quadallelic_count,"quadallelic sites."))
  print(paste("Alignment contains",pentallelic_count,"pentallelic sites."))

  split_support <- split_support %>%
    mutate(Alignment_Name = alignment_name) %>% # Add alignment name
    rename(Alignment_Position = Zeroed_Site_Position) %>%
    mutate(Alignment_Position = as.integer(Alignment_Position) + 1) %>% # Convert position to 1-based
    select(Alignment_Name,Alignment_Position,Site_Pattern,Non_Base_Taxa,Singleton_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5)
  return(split_support)
}
