#' Get Trimmed Tree
#'
#' Given a tree and a list of taxa, this function returns a tree pruned to include only those taxa. Note: Function will STOP if all taxa are not in tree.
#' @param whole_tree Phylo object
#' @param taxa Vector containing desired taxa list
#' @return Pruned phylo object
#' @export
#' @examples
#' getTrimmedTree(whole_tree,taxa)
#'
getTrimmedTree <- function(whole_tree,taxa){

  # Get tree species
  tree_species <- sort(whole_tree$tip.label)

  # If tree species match subset list, return whole tree
  if(identical(tree_species,sort(taxa))){
    return(whole_tree)
  }

  # Otherwise, check if all species from subset list are in tree. If so, return pruned tree.
  if(checkTreeTaxa(whole_tree,taxa)){
    pruned_tree <- ape::drop.tip(whole_tree,whole_tree$tip.label[-match(taxa, whole_tree$tip.label)])
    return(pruned_tree)
  }

  else{
    stop("ERROR: Some taxa from reqested list not in tree.")
  }
}
