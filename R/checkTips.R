#' Rboreturm Tip Checker
#'
#' This function returns TRUE if all 'taxa' are present in 'tree' (or a multiPhylo of trees); FALSE otherwise
#' @param tree phylo or multiPhylo object
#' @param taxa Character vector containing taxa to query
#' @return TRUE if all 'taxa' in 'tree'; else, FALSE
#' @export

checkTips <- function(tree,taxa){
  
  if(Rboretum::isPhylo(tree)){
    if(all(taxa %in% tree$tip.label)){
      return(TRUE)
    } else{ return(FALSE) }
  } else if(Rboretum::isMultiPhylo(tree)){
    species_check <- c()
      
    for(i in 1:length(tree)){
      species_check <- c(species_check,all(taxa %in% tree[[i]]$tip.label))
    }
      
    if(all(species_check)){
      return(TRUE)
    } else{ return(FALSE) }
  }
  else{
    stop("'tree' does not appear to be a valid phylo or multiPhylo object")
  }
}