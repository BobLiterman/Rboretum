#' Rboretum Common Taxa Fetcher
#'
#' This function takes a multiPhylo object and returns a sorted list of taxa common to all trees. Note: Function will STOP if <3 taxa are shared among trees.
#' @param trees multiPhylo object where >= 3 species are shared among all trees
#' @return Sorted character vector of tip labels present in all trees
#' @export

getSharedTaxa <- function(trees){
  
  if(!Rboretum::isMultiPhylo(trees,check_three_taxa = TRUE)){
    stop("'trees' does not appear to be a valid multiPhylo where all trees share at least three taxa.")
  } else{
    
    all_species <- purrr::map(.x=trees,.f=function(x){x$tip.label}) %>% unlist() %>% unique() %>% sort()
    tip_table  <- purrr::map(.x=trees,.f=function(x){x$tip.label}) %>% unlist() %>% table()
  
    shared_species <- all_species[tip_table==length(trees)]
  
    return(shared_species)
  }
}
