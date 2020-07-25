#' Rboretum Alignment Species Fetcher
#'
#' Given the path to an alignment, this function returns the taxon labels as a sorted, semicolon-separated character
#' @param alignment_path Relative or absolute path to alignment file
#' @return Sorted, semicolon-separated list of taxa
#' @export

getAlignmentSpecies <- function(alignment_path){
  
  # Ensure that a path is provided as a character
  if(missing(alignment_path)){
    stop("getAlignmentSpecies requires an 'alignment_path' argument...")
  } else if(!is.character(alignment_path)){
    stop("'alignment_path' should be a character path to an alignment file...")
  } else if(length(alignment_path) != 1){
    stop("getAlignmentSpecies can only accept a single path as an argument...")
  }

  # Ensure file exists
  isFile <- file.exists(alignment_path) & !dir.exists(alignment_path)
  
  if(isFile){
    alignment_path <- file_path_as_absolute(alignment_path)
  }
  else{
    stop("'alignment_path' doesn't point to a valid file...")
  }
  
  alignment_species <- fetchAlignmentSpecies(alignment_path)
  
  # Catch errors
  if(length(alignment_species) == 1){
    stop("No species could be extracted from file at 'alignment_path'...")
  } else{
    alignment_species <- Rboretum::semiSorter(alignment_species)
  }
  
  return(alignment_species)
}