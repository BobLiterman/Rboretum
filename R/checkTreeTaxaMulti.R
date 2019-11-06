#' Rboreturm Multiphlyo Taxon Checker
#'
#' This function returns TRUE if all 'taxa' are present in all 'trees', FALSE otherwise
#' @param trees Multiphylo object
#' @param taxa Character vector containing taxa to query
#' @return TRUE if all 'taxa' in all 'trees', else, FALSE
#' @export
#' @examples
#' 
#' myTrees <- c(treeA,treeB,treeC)
#' taxa_to_find <- c('Spp1','Spp2','Spp3')
#' checkTreeTaxa(myTrees,taxa_to_find)
#'

checkTreeTaxaMulti <- function(trees,taxa){
  
  # Check that input is multiphylo and has at least 2 trees
  if(has_error(unlist(attributes(trees)$class))){ 
    stop("'trees' argument should be a multiPhylo object")
  } else if(!"multiPhylo" %in% unlist(attributes(trees)$class)){
    stop("'trees' argument should be a multiPhylo object")
  } else if(length(trees)<2){
    stop("At least two trees are required for comparison. For a single tree, use checkTreeTaxa()")
  }
  
  # Check that all trees contain all species
  species_check <- c()
  for(i in 1:length(trees)){
    species_check <- c(species_check,Rboretum::checkTreeTaxa(trees[[i]],taxa))
  }

  if(all(species_check)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}
