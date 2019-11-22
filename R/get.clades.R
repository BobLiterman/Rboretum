#' Rboretum Clade Fetcher
#'
#' This function takes a tree and returns a sorted character vector containing each monophyletic group
#' @param tree Rooted phylo object
#' @return Character vector of semicolon-separated monophyletic clades
#' @export
#' @examples
#' get.clades(tree)
#'

get.clades <- function(tree){
  
  if(has_error(ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for get.clades")}
  
  splits <- Rboretum::get.splits(tree)
  
  root_split <- splits %>% filter(is.na(Split_Node))
  
  clades <- splits %>% 
    filter(!is.na(Split_Node)) %>%
    pull(Clade) %>%
    as.character()
  
  clades <- c(clades,as.character(root_split$Clade),as.character(root_split$Mirror_Clade)) %>%
    sort()
  
  # Don't return single taxon clades
  clades <- clades[str_detect(clades,";")]
  
  return(clades)
}
