#' Read Rooted Phylogenetic Tree from Path
#'
#' This function returns a rooted tree, assuming the tree contains the outgroup taxa, and they are monophyletic.
#' @param tree_path Path to tree file
#' @param root_taxa Vector containing outgroup species IDs (Must be in tree and monophyletic)
#' @return A phylo object rooted at specified taxa
#' @export
#' @examples
#' readRootedPhylo('/path/to/tetrapod_tree.nwk',c('Fish_1','Fish_2'))
#' OR
#' tree_path <- '/path/to/tetrapod_tree.nwk'
#' outgroup_species <- c('Fish_1','Fish_2')
#' readRootedPhylo(tree_path,outgroup_species)

readRootedPhylo <- function(tree_path,root_taxa){

  # Read in tree and fetch species
  raw_tree <- ape::read.tree(tree_path)
  tree_species <- raw_tree$tip.label

  # Ensure root species in tree, and are monophyletic for rooting
  if(!(all(root_taxa %in% tree_species))){
    stop("Root taxa not found in tree.")
  }

  else if(!has_error(root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
    return(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))
  }
  else{
    stop("Ape cannot root tree on these taxa.")
  }
}
