#' Check if Trees Share All Species
#'
#' This function returns TRUE if 'tree_1' and 'tree_2' have all species in commmon
#' @param tree_1 Phylo object
#' @param tree_2 Phylo object
#' @return TRUE if 'tree_1' and 'tree_2' have all species in commmon; else, FALSE
#' @export
#' @examples
#' checkIdenticalSpecies(tree_1,tree_2)
#'

# Returns TRUE if tree_1 and tree_2 share at least three species
checkIdenticalSpecies <- function(tree_1,tree_2){

  species_1 <- sort(tree_1$tip.label)
  species_2 <- sort(tree_2$tip.label)

  if(all(species_1 %in% species_2) & all(species_2 %in% species_1)){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
