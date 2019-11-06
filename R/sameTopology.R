#' Check for Identical Tree Topology
#'
#' This function returns TRUE if 'tree_1' and 'tree_2' have the same topology after pruning to identical species lists
#' @param tree_1 Rooted phylo object
#' @param tree_2 Rooted phylo object
#' @return TRUE if 'tree_1' and 'tree_2' have the same topology after pruning to identical species lists; else, FALSE
#' @export
#' @examples
#' sameTopology(tree_1,tree_2)
#'

# Returns TRUE if tree_1 and tree_2 have identical topology (after pruning if necessary)
sameTopology <- function(tree_1,tree_2){

  shared_species <- Rboretum::getSharedSpecies(c(tree_1,tree_))

  # If trees have identical species lists, compare topologies
  if(identical(species_1,species_2)){

    if(length(species_1)<3){
      stop("Trees contain fewer than three species.")
    }

    else{
      if(ape::all.equal.phylo(tree_1,tree_2,use.edge.length = FALSE)){
        return(TRUE)}
      else{
        return(FALSE)}
    }
  }

  # If trees have different species lists, check >3 and prune
  else{
    if(has_error(Rboretum::checkSharedSpecies(tree_1,tree_2)) | !Rboretum::checkSharedSpecies(tree_1,tree_2)){
      stop("Trees contain fewer than three shared species.")
    } else{
      shared_species <- Rboretum::getSharedSpecies(tree_1,tree_2)
    }

    pruned_tree_1 <- Rboretum::getTrimmedTree(tree_1,shared_species)
    pruned_tree_2 <- Rboretum::getTrimmedTree(tree_2,shared_species)

    if(ape::all.equal.phylo(pruned_tree_1,pruned_tree_2,use.edge.length = FALSE)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}
