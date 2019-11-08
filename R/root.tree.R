#' Rboretum Tree Rooter
#'
#' This ape wrapper accepts a phylo object and a character vector of new root taxa and returns a rooted phylo object, if possible
#' @param tree Phylo object
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree and monophyletic)
#' @return A phylo object rooted at specified taxa
#' @export
#' @examples
#' root.tree(birdTree,c('Alligator','Turtle'))
root.tree <- function(tree,root_taxa){
  
  if(has_error(ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!(all(root_taxa %in% tree$tip.label))){
    stop("Root taxa not found in tree.")
  } else if(has_error(ape::root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
    stop("Ape cannot root tree on these taxa.")}
  
  return(ape::root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))
}
