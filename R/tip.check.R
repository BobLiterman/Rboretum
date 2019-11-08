#' Rboreturm Phylo Taxon Checker
#'
#' This function returns TRUE if all 'taxa' are present in 'tree', FALSE otherwise
#' @param tree Phylo object
#' @param taxa Character vector containing taxa to query
#' @return TRUE if all 'taxa' in 'tree', else, FALSE
#' @export
#' @examples
#' 
#' taxa_to_find <- c('Spp1','Spp2','Spp3')
#' tip.check(myTree,taxa_to_find)
#'

tip.check <- function(tree,taxa){
  tree_species  <- tree$tip.label
  if(all(taxa %in% tree_species)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}
