#' Rboretum MultiPhylo Autonamer
#'
#' This function adds/replaces names of multiPhylo trees with dummy names (Tree_1, Tree_2, etc.)
#' @param trees multiPhylo object
#' @return A named multiPhylo object
#' @export

treeNamer <- function(trees){
  
  # Check if 'trees' is a valid multiPhylo object
  if(!Rboretum::isMultiPhylo(trees)){
    stop("'trees' does not appear to be a valid multiPhylo object")
  }
  
  tree_numbers <- 1:length(trees)
  tree_names <- purrr::map(.x = tree_numbers,.f=function(x){paste(c('Tree',x),collapse = '_')}) %>% unlist()
  names(trees) <- tree_names
  return(trees)
}