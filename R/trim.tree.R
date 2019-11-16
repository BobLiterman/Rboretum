#' Rboretum Tree Trimmer
#'
#' This ape wrapper returns a phylo or multiPhylo object that has been pruned to include only specified taxa. Note: Function will STOP if any 'tree' is missing any 'taxa'.
#' @param tree phylo or multiPhylo object
#' @param taxa Character vector of desired tip labels to keep (or discard if remove=TRUE)
#' @param remove OPTIONAL: If TRUE, tip labels specified by 'taxa' are removed from all trees rather than retained [Default: FALSE, prune 'tree' down to 'taxa']
#' @return Pruned phylo or multiPhylo object
#' @export
#' @examples
#' myFavoriteSpecies  <- c('Spp1','Spp2','Spp3')
#' trim.tree(myBigTree,myFavoriteSpecies)
#' 
#' myLeastFavoriteSpecies  <- c('Spp4','Spp5','Spp6')
#' trim.tree(myBigTree,myLeastFavoriteSpecies,remove=TRUE)
#' 
#' myHugeTrees <- c(BigTree_1,BigTree_2,BigTree_3)
#' trim.tree(myHugeTrees,myFavoriteSpecies)
trim.tree <- function(tree,taxa,remove){
  
  if(missing(remove)){
    remove <- FALSE
  } else if(!is.logical(remove)){
    remove <- FALSE
  }
  
  if(has_error(unlist(attributes(tree)$class))){ 
    stop("'tree' argument should be a phylo or multiPhylo object")
  } else if(!"phylo" %in% unlist(attributes(tree)$class) & !"multiPhylo" %in% unlist(attributes(tree)$class)){
    stop("'tree' argument should be a phylo or multiPhylo object")}
  
  if(tip.check(tree,taxa)){
    
    if(!remove){
      
      if("phylo" %in% unlist(attributes(tree)$class)){
        return(ape::drop.tip(tree,tree$tip.label[-match(taxa, tree$tip.label)]))
      }
      
      if("multiPhylo" %in% unlist(attributes(tree)$class)){
        for(i in 1:length(tree)){
          tree[[i]] <- ape::drop.tip(tree[[i]],tree[[i]]$tip.label[-match(taxa, tree[[i]]$tip.label)])
        }
        return(tree)
      }
      
    } else{
      
      if("phylo" %in% unlist(attributes(tree)$class)){
        keep_taxa <- tree$tip.label[!tree$tip.label %in% taxa]
        return(ape::drop.tip(tree,tree$tip.label[-match(keep_taxa, tree$tip.label)]))
      }
      
      if("multiPhylo" %in% unlist(attributes(tree)$class)){
        for(i in 1:length(tree)){
          keep_taxa <- tree[[i]]$tip.label[!tree[[i]]$tip.label %in% taxa]
          tree[[i]] <- ape::drop.tip(tree[[i]],tree[[i]]$tip.label[-match(keep_taxa, tree[[i]]$tip.label)])
        }
        return(tree)
    }
  } 
  } else{ stop("ERROR: One or more trees missing one or more taxa.") }
}

  
  
  
  
  
  
