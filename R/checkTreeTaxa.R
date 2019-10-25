#' Check Tree Taxa
#'
#' This function returns TRUE if all 'taxa' are present in 'tree', FALSE otherwise
#' @param tree Phylo object
#' @param taxa Vector containing taxa to check
#' @return TRUE if 'taxa' in 'tree', else, FALSE
#' @export
#' @examples
#' checkTreeTaxa(tree,taxa_to_find)
#'

checkTreeTaxa <- function(tree,taxa){
  tree_species  <- tree$tip.label
  if(all(taxa %in% tree_species)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}
