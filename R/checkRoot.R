#' Rboretum Root Checker
#'
#' This function returns TRUE if  'taxa' are part of a root split
#' @param tree Phylo object
#' @param taxa Character vector containing taxa to check
#' @return TRUE if 'taxa' are either side of the root, else, FALSE
#' @export
#' @examples
#' checkRoot(tree,taxa_to_check)
#'

checkRoot <- function(tree,taxa){
  if(!(all(taxa %in% tree$tip.label))){
    stop("Taxa missing from tree.")
  }

  # Get tree species
  tree_species <- sort(tree$tip.label)
  taxa <- sort(taxa)

  # Get mirror clade (all species - focal group)
  mirror_taxa <- sort(dplyr::setdiff(tree_species, taxa))

  if(ape::is.monophyletic(tree,taxa) & ape::is.monophyletic(tree,mirror_taxa)){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
