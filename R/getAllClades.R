#' Rboretum Clade Fetcher
#'
#' This function takes a tree and returns a sorted character vector containing each monophyletic group
#' @param tree Phylo object
#' @return Character vector of semicolon-separated monophyletic clades
#' @export
#' @examples
#' getAllClades(tree)
#'

getAllClades <- function(tree){
  
  if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for getAllClades")
  }
  
  splits <- Rboretum::getSplits(tree)
  
  root_split <- splits %>% filter(is.na(Split_Node))
  
  clades <- splits %>% 
    filter(!is.na(Split_Node)) %>%
    pull(Clade) %>%
    as.character()
  
  clades <- c(clades,as.character(root_split$Clade),as.character(root_split$Mirror_Clade)) %>%
    sort()
  return(clades)
}
