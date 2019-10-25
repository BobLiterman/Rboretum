#' Re-Root Phylogenetic Tree
#'
#' This function returns a rooted tree, assuming the tree contains the outgroup taxa, and they are monophyletic.
#' @param tree Path to tree file
#' @param root_taxa Vector containing outgroup species IDs (Must be in tree and monophyletic)
#' @return A phylo object rooted at specified taxa
#' @export
#' @examples
#' rerootTree(tree,taxa)

rerootTree <- function(tree,root_taxa){

  # Read in tree and fetch species
  tree_species <- tree$tip.label

  # Ensure root species in tree
  if(!(all(root_taxa %in% tree_species))){
    stop("Root taxa not found in tree.")
  }

  # If root species in tree, find MRCA and root at that node
  else if(!has_error(root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
    return(ape::root.phylo(tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))
  }
  else{
    stop("Ape cannot root tree on these taxa.")
  }
}
