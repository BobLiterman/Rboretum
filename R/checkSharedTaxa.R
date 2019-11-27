#' Rboretum Three Common Taxa Checker
#'
#' This function takes a multiPhylo object and returns TRUE if all trees share at least three common taxa (e.g. enough to make a tree); otherwise, FALSE
#' @param trees multiPhylo object
#' @return TRUE (all trees share >= 3 species) or FALSE
#' @export

checkSharedTaxa <- function(trees){

  # Check that input is multiphylo and has at least 2 trees
  if(!Rboretum::isMultiPhylo(trees)){
    stop("'trees' does not appear to be a multiPhylo object with 2+ trees.")
  } 
  
  if(has_error(Rboretum::getSharedTaxa(trees))){
    return(FALSE)
  } else if(length(Rboretum::getSharedTaxa(trees))<3){
    return(FALSE)
  } else if(length(Rboretum::getSharedTaxa(trees))>=3){
    return(TRUE)
  } else { stop("Something unexpected happened...")} 
}