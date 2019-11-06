#' Rboretum Clade Fetcher
#'
#' This function takes a tree and returns a vector containing each monophyletic group
#' @param tree Phylo object
#' @return Character vector of monophyletic clades
#' @export
#' @examples
#' getAllClades(tree)
#'

getAllClades <- function(tree){
  
  splits <- Rboretum::getSplits(tree)

}
