#' Rboretum Two-Tree Comparability Checker
#'
#' This function returns TRUE if 'tree_1' and 'tree_2' have 3 or more species in commmon, and a unique topology
#' @param tree_1 Phylo object
#' @param tree_2 Phylo object
#' @return TRUE if 'tree_1' and 'tree_2' have 3 or more species in commmon and a unique topology; else, FALSE
#' @export
#' @examples
#' checkComparable(tree_1,tree_2)
#'

checkComparable <- function(tree_1,tree_2){
  if(!has_error(Rboretum::sameTopology(tree_1,tree_2))){
    return(FALSE)
  } else if(!Rboretum::sameTopology(tree_1,tree_2)){
    return(TRUE)
  } else{ return(FALSE) }
}
