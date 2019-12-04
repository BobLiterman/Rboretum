#' Rboreturm Tip Checker
#'
#' This function returns TRUE if all 'taxa' are present in 'tree' (or a multiPhylo of trees); FALSE otherwise
#' @param tree phylo or multiPhylo object
#' @param taxa Character vector containing taxa to query
#' @return TRUE if all 'taxa' in 'tree'; else, FALSE
#' @export

checkTips <- function(tree,taxa){
  
  if(!Rboretum::isPhylo(tree) & !Rboretum::isMultiPhylo(tree)){
    stop("'tree' does not appear be a valid phylo or multiPhylo object.")
  }
  
  if(Rboretum::isPhylo(tree)){
    
    if(all(taxa %in% tree$tip.label)){
      return(TRUE)
    } else{ 
      return(FALSE) 
    }
    
  } else if(Rboretum::isMultiPhylo(tree)){
    
    if(purrr::map(.x = tree,.f = function(x){all(taxa %in% x$tip.label)}) %>% unlist() %>% all()){
      return(TRUE)
    } else{ return(FALSE) }
  }
}