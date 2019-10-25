#' Get Shared Species From Two Trees
#'
#' Returns sorted list of species that occur in 'tree_1' and 'tree_2'
#' @param tree_1 Phylo object
#' @param tree_2 Phylo object
#' @return Vector of species shared between 'tree_1' and 'tree_2
#' @export
#' @examples
#' getSharedSpecies(tree_1,tree_2)
#'

getSharedSpecies <- function(tree_1,tree_2){

  species_1 <- sort(tree_1$tip.label)
  species_2 <- sort(tree_2$tip.label)

  shared_species <- Reduce(intersect, list(species_1,species_2))

  if(length(shared_species ) > 0){
    return(sort(shared_species))
  }
  else{
    stop("Trees share no species.")
  }
}
