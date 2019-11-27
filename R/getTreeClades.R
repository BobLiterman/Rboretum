#' Rboretum Clade Fetcher
#'
#' This function takes a tree and returns a sorted character vector containing each monophyletic group (not including the root clades)
#' @param tree Rooted phylo object
#' @return Character vector of semicolon-separated monophyletic clades
#' @export
#'

getTreeClades <- function(tree){
  
  if(has_error(ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for getTreeClades")}
  
  clades <- Rboretum::getTreeSplits(tree) %>%
    filter(!is.na(Split_Node)) %>%
    pull(Clade) %>%
    as.character() %>%
    sort()

  return(clades)
}
