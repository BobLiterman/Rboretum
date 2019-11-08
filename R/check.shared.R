#' Rboretum Three Common Taxa Checker
#'
#' This function takes a multiPhylo object and returns TRUE if all trees share at least three common taxa (e.g. enough to make a tree); otherwise, FALSE
#' @param trees multiPhylo object
#' @return TRUE (all trees share >= 3 species) or FALSE
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#' check.shared(trees)
#'
check.shared <- function(trees){

  # Check that input is multiphylo and has at least 2 trees
  if(has_error(unlist(attributes(trees)$class))){ 
    stop("'trees' argument should be a multiPhylo object")
  } else if(!"multiPhylo" %in% unlist(attributes(trees)$class)){
    stop("'trees' argument should be a multiPhylo object")
  } else if(length(trees)<2){
    stop("At least two trees are required for comparison.")
  } 
  
  if(has_error(Rboretum::get.shared(trees))){
    return(FALSE)
  } else if(length(Rboretum::get.shared(trees))<3){
    return(FALSE)
  } else if(length(Rboretum::get.shared(trees))>=3){
    return(TRUE)
  } else { stop("Something unexpected happened...")} 
}