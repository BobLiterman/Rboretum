#' Rboretum Node Label Stripper
#'
#' This function takes a phylo or multiPhylo of trees and removed node labels
#' @param tree A phylo or multiPhylo object
#' @return A phylo or multiPhylo object with no node labels
#' @export
#' 
stripNodeLabels <- function(tree){
  
  if(missing(tree)){
    stop("'stripNodeLabels' requires a phylo or multiPhylo object.")
  } else if(!Rboretum::isMultiPhylo(tree) & !Rboretum::isPhylo(tree)){
    stop("'stripNodeLabels' requires a phylo or multiPhylo object.")
  }
  
  if(Rboretum::isPhylo(tree)){
    tree$node.label <- NULL
    return(tree)
  }
  
  if(Rboretum::isMultiPhylo(tree)){
    tree_count <- length(tree)
    for(i in 1:tree_count){
      no_label_tree <- tree[[i]]
      no_label_tree$node.label <- NULL
      tree[[i]] <- no_label_tree
    }
    return(tree)
  }
}