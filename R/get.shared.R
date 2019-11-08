#' Rboretum Common Taxa Fetcher
#'
#' This function takes a multiPhylo object and returns a sorted list of taxa common to all trees
#' @param trees multiPhylo object
#' @return Sorted character vector of tip labels present in all trees
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#' get.shared(trees)
#'
get.shared <- function(trees){
  
  # Check that input is multiphylo and has at least 2 trees
  if(has_error(unlist(attributes(trees)$class))){ 
    stop("'trees' argument should be a multiPhylo object")
  } else if(!"multiPhylo" %in% unlist(attributes(trees)$class)){
    stop("'trees' argument should be a multiPhylo object")
  } else if(length(trees)<2){
    stop("At least two trees are required for comparison.")
  }
  
  tree_species <- c()
  for(i in 1:length(trees)){
    tree_species <- c(tree_species,trees[[i]]$tip.label)
  }
  
  all_species <- sort(unique(tree_species))
  shared_species <- all_species[table(tree_species)==length(trees)]
  
  if(length(shared_species)>0){
    return(shared_species)
  } else{ stop('No species common to all trees.')}
}
