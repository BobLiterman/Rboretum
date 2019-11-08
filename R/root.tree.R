#' Rboretum Tree Re-Rooter
#'
#' This ape wrapper accepts a phylo object and a character vector of new root taxa and returns a rooted phylo object, if possible
#' @usage newRootTaxa <- c('NewRoot_1','NewRoot_2')
#' rerootedTree  <- root.tree(myTree,newRootTaxa)
#' @param tree Phylo object
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree and monophyletic)
#' @return A phylo object rooted at specified taxa
#' @export
#' @examples
#' root.tree(birdTree,c('Alligator','Turtle'))

root.tree <- function(tree,root_taxa){

  # Read in tree and fetch species
  tree_species <- tree$tip.label

  # Ensure root species in tree
  if(!(all(root_taxa %in% tree_species))){
    stop("Root taxa not found in tree.")
  }

  # If root species in tree, find MRCA and root at that node
  else if(!has_error(ape::root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
    return(ape::root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))
  }
  else{
    stop("Ape cannot root tree on these taxa.")
  }
}
