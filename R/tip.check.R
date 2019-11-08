#' Rboreturm Taxon Checker
#'
#' This function returns TRUE if all 'taxa' are present in 'tree' (or a multiphylo of trees); FALSE otherwise
#' @param tree Phylo or Multiphylo object
#' @param taxa Character vector containing taxa to query
#' @return TRUE if all 'taxa' in 'tree'; else, FALSE
#' @export
#' @examples
#' 
#' taxa_to_find <- c('Spp1','Spp2','Spp3')
#' tip.check(myTree,taxa_to_find)
#' 
#' trees_to_check <- c(Tree_1,Tree_2,Tree_3)
#' tip.check(trees_to_check,taxa_to_find)
#'

tip.check <- function(tree,taxa){
  
  # Check that input is multiphylo and has at least 2 trees
  if(has_error(unlist(attributes(trees)$class))){ 
    stop("'trees' argument should be a phylo or multiPhylo object")
  } else if("phylo" %in% unlist(attributes(trees)$class)){
    if(all(taxa %in% tree$tip.label)){
      return(TRUE)
    } else{ return(FALSE) }
  } else if("multiPhylo" %in% unlist(attributes(trees)$class)){
    if(length(tree)>=2){
      
      species_check <- c()
      
      for(i in 1:length(trees)){
        species_check <- c(species_check,all(taxa %in% tree[[i]]$tip.label))
      }
      
      if(all(species_check)){
        return(TRUE)
      }
      else{ return(FALSE) }
    } else{
      if(all(taxa %in% tree[[1]]$tip.label)){
        return(TRUE)
      } else{ return(FALSE) }
    }
  } else{ stop("'tree' argument should be a phylo or multiPhylo object") }
}
  
  
  
  
  
  
  

