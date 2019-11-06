#' Rboretum Tree Trimmer
#'
#' This ape wrapper returns a phylo object that has been pruned to include only specified taxa. Note: Function will STOP if tree is missing taxa.
#' @param tree Phylo object
#' @param taxa Character vector of desired tip labels to keep
#' @return Pruned phylo object
#' @export
#' @examples
#' getTrimmedTree(tree,taxa)
#'
getTrimmedTree <- function(tree,taxa){

  # Get tree species
  tree_species <- sort(tree$tip.label)

  # If tree species match subset list, return whole tree
  if(identical(tree_species,sort(taxa))){
    return(tree)
  }

  # Otherwise, check if all species from subset list are in tree. If so, return pruned tree.
  if(checkTreeTaxa(tree,taxa)){
    pruned_tree <- ape::drop.tip(tree,tree$tip.label[-match(taxa, tree$tip.label)])
    return(pruned_tree)
  } else{ stop("ERROR: Some taxa from reqested list not in tree.")}
}
