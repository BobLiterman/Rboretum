#' Rboretum Three Common Taxa Checker
#'
#' This function takes a multiphylo object and returns TRUE if all trees share at least three common taxa (e.g. enough to make a tree); otherwise, FALSE
#' @param trees Multiphylo object
#' @return TRUE (all trees share >= 3 species) or FALSE
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#' checkSharedSpeciesMulti(trees)
#'
checkSharedSpecies <- function(trees){

  # Check that input is multiphylo and has at least 2 trees
  if(has_error(unlist(attributes(trees)$class))){ 
    stop("'trees' argument should be a multiPhylo object")
  } else if(!"multiPhylo" %in% unlist(attributes(trees)$class)){
    stop("'trees' argument should be a multiPhylo object")
  } else if(length(trees)<2){
    stop("At least two trees are required for comparison. For a single tree, use checkTreeTaxa()")
  }
  
  # Check that all trees contain at least three taxa
  species_check <- c()
  for(i in 1:length(trees)){
    temp_tree <- trees[[i]]
    species_check <- c(species_check,temp_tree$tip.label)
  }
  
  all_species <- sort(unique(species_check))
  shared_species <- all_species[table(species_check)==length(trees)]
  if(length(shared_species)>=3){
    return(TRUE)
  }
  else{return(FALSE)}
}
