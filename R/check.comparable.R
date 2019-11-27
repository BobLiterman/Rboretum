#' Rboretum Two-Tree Comparability Checker
#'
#' This function returns TRUE if 'tree_1' and 'tree_2' have 3 or more species in common, and a unique topology
#' @param tree_1 phylo object
#' @param tree_2 phylo object
#' @return TRUE if 'tree_1' and 'tree_2' have 3 or more species in common and a unique topology; else, FALSE
#' @export
#' @examples
#' check.comparable(tree_1,tree_2)
#'

check.comparable <- function(tree_1,tree_2){
  if(has_error(Rboretum::sameTopology(c(tree_1,tree_2)))){
    return(FALSE)
  } else if(!Rboretum::sameTopology(c(tree_1,tree_2))){
    return(TRUE)
  } else{ return(FALSE) }
}
