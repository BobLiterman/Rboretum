#' Rboretum Multiphylo Rooted Tree Reader
#'
#' This function is an ape wrapper, and returns a multiPhylo object, with each tree rooted as specified by the user.
#' @param tree_paths Vector of paths to tree files that can be read by ape::read.tree() or ape::read.nexus()
#' @param root_taxa Character vector containing outgroup species IDs (Must be in all trees and always monophyletic)
#' @param tree_names OPTIONAL: Character vector of names to assign to trees. Length must equal the number of trees.
#' @return A multiPhylo object, with each tree rooted at specified taxa
#' @export

batch_readRootedTree <- function(tree_paths,root_taxa,tree_names){
  
  if(length(tree_paths)<2){
    stop("'tree_paths' contains fewer than 2 items, use readRootedTree()")
  }
  
  trees <- purrr::map(.x = tree_paths,root_taxa=root_taxa, .f = readRootedTree)
  class(trees) <- "multiPhylo"
  
  if(missing(tree_names)){
    return(trees)
  } else if(length(tree_names) != length(trees)){
    print("Incorrect number of tree names provided. Returning unnamed trees.")
    return(trees)
  } else{
    names(trees) <- tree_names
    return(trees)
  }
}
