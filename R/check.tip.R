#' Rboreturm Taxon Checker
#'
#' This function returns TRUE if all 'taxa' are present in 'tree' (or a multiPhylo of trees); FALSE otherwise
#' @param tree phylo or multiPhylo object
#' @param taxa Character vector containing taxa to query
#' @return TRUE if all 'taxa' in 'tree'; else, FALSE
#' @export
#' @examples
#' 
#' taxa_to_find <- c('Spp1','Spp2','Spp3')
#' check.tip(myTree,taxa_to_find)
#' 
#' trees_to_check <- c(Tree_1,Tree_2,Tree_3)
#' check.tip(trees_to_check,taxa_to_find)
#'

check.tip <- function(tree,taxa){
  
  if(has_error(unlist(attributes(tree)$class))){ 
    stop("'tree' argument should be a phylo or multiPhylo object")
  } else if(!"phylo" %in% unlist(attributes(tree)$class) & !"multiPhylo" %in% unlist(attributes(tree)$class)){
    stop("'tree' argument should be a phylo or multiPhylo object")}
  
  if("phylo" %in% unlist(attributes(tree)$class)){
    if(all(taxa %in% tree$tip.label)){
      return(TRUE)
    } else{ return(FALSE) }
  } 
  
  if("multiPhylo" %in% unlist(attributes(tree)$class)){
      species_check <- c()
      
      for(i in 1:length(tree)){
        species_check <- c(species_check,all(taxa %in% tree[[i]]$tip.label))
      }
      
      if(all(species_check)){
        return(TRUE)
      } else{ return(FALSE) }
  }
}