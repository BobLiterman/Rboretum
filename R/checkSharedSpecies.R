#' Check for Shared Species
#'
#' This function returns TRUE if 'tree_1' and 'tree_2' have 3 or more species in commmon
#' @param tree_1 Phylo object
#' @param tree_2 Phylo object
#' @return TRUE if 'tree_1' and 'tree_2' have 3 or more  species in commmon; else, FALSE
#' @export
#' @examples
#' checkSharedSpecies(tree_1,tree_2)
#'

# Returns TRUE if tree_1 and tree_2 share at least three species
checkSharedSpecies <- function(tree_1,tree_2){

  species_1 <- sort(tree_1$tip.label)
  species_2 <- sort(tree_2$tip.label)

  shared_species <- Reduce(intersect, list(species_1,species_2))

  if(length(shared_species) >= 3){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
