#' Rboretum Identical Taxon Checker
#'
#' This function takes a multiPhylo object and returns TRUE if all species are shared among all trees; otherwise, FALSE
#' @param trees multiPhylo object
#' @return TRUE if all trees have all species in common; else, FALSE
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3,...)
#' same.taxa(trees)
#'
same.taxa <- function(trees){

  # Check that input is multiphylo and has at least 2 trees
  if(!Rboretum::is.multiPhylo(trees)){
    stop("'trees' does not appear to be a valid multiPhylo object with 2+ trees")
  }

  # Check that all trees contain all species
  species_lists <- c()
  
  for(i in 1:length(trees)){
    species_lists <- c(species_lists,paste(sort(trees[[i]]$tip.label),collapse = "_"))
  }

  if(length(unique(species_lists))==1){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}
