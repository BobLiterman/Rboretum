#' Get Split Node
#'
#' This function returns the node ID of the MRCA of a set of 'taxa' in 'tree'
#' @param tree Phylo object
#' @param taxa Vector containing taxa to check
#' @return Node ID of MRCA
#' @export
#' @examples
#' getSplitNode(tree,taxa_to_check)
#'

# Given a tree and a vector of taxa, return node of MRCA (if present and monophyletic)
getSplitNode <- function(tree,taxa){

  # Ensure speecies are in the tree
  if(has_error(Rboretum::checkMonophyletic(tree,taxa)) | !Rboretum::checkMonophyletic(tree,taxa)){
    stop('Taxa are not monophyletic, or not in the tree.')
  }

  # If so, return split node ID
  return(ape::getMRCA(tree,taxa))
}
