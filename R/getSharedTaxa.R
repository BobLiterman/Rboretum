#' Rboretum Common Taxa Fetcher
#'
#' This function takes a multiPhylo object and returns a sorted list of taxa common to all trees. Note: Function will STOP if <3 taxa are shared among trees.
#' @param trees multiPhylo object where >= 3 species are shared among all trees
#' @return Sorted character vector of tip labels present in all trees
#' @export

getSharedTaxa <- function(trees){
  
  # Check that input is multiphylo and has at least 2 trees
  if(!Rboretum::isMultiPhylo(trees, check_three_taxa = TRUE)){
    stop("'trees' does not appear to be a valid multiPhylo object with 2+ trees")
  }

  all_species <- purrr::map(.x=trees,.f=function(x){x$tip.label}) %>% unlist() %>% unique() %>% sort()
  shared_species <- all_species[table(purrr::map(.x=trees,.f=function(x){x$tip.label}) %>% unlist())==length(trees)]
  
  return(shared_species)
}
