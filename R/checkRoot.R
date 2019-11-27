#' Rboretum Root Checker
#'
#' This function returns TRUE if  'taxa' are part of a root split
#' @param tree Rooted phylo object
#' @param taxa Character vector containing taxa to check
#' @return TRUE if 'taxa' are either side of the root, else, FALSE
#' @export

checkRoot <- function(tree,taxa){

  if(has_error(ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for checkRoot")
  } else if(!(all(taxa %in% tree$tip.label))){
    stop("Specified taxa missing from tree.")}

  # Get mirror clade (all species - focal group)
  mirror_taxa <- sort(dplyr::setdiff(tree$tip.label, taxa))

  if(ape::is.monophyletic(tree,taxa) & ape::is.monophyletic(tree,mirror_taxa)){
    return(TRUE)
  } else{ return(FALSE) }
}
