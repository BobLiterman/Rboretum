#' Rboretum Tree Trimmer
#'
#' This ape wrapper returns a phylo or multiPhylo object that has been pruned to include only specified taxa. Note: Function will STOP if any tree is missing taxa.
#' @param tree phylo or multiPhylo object
#' @param taxa Character vector of desired tip labels to keep
#' @return Pruned phylo or multiPhylo object
#' @export
#' @examples
#' myFavoriteSpecies  <- c('Spp1','Spp2','Spp3')
#' trim.tree(myBigTree,myFavoriteSpecies)
#' 
#' myHugeTrees <- c(BigTree_1,BigTree_2,BigTree_3)
#' trim.tree(myHugeTrees,myFavoriteSpecies)
trim.tree <- function(tree,taxa){
  
  if(has_error(unlist(attributes(tree)$class))){ 
    stop("'tree' argument should be a phylo or multiPhylo object")
  } else if(!"phylo" %in% unlist(attributes(tree)$class) & !"multiPhylo" %in% unlist(attributes(tree)$class)){
    stop("'tree' argument should be a phylo or multiPhylo object")}
  
  if(tip.check(tree,taxa)){
    
    if("phylo" %in% unlist(attributes(tree)$class)){
      return(ape::drop.tip(tree,tree$tip.label[-match(taxa, tree$tip.label)]))
    }
    
    if("multiPhylo" %in% unlist(attributes(tree)$class)){
      for(i in 1:length(tree)){
        tree[[i]] <- ape::drop.tip(tree[[i]],tree[[i]]$tip.label[-match(taxa, tree[[i]]$tip.label)])
      }
      return(tree)
    }
  } else{ stop("ERROR: One or more trees missing one or more taxa.") }
}

  
  
  
  
  
  
