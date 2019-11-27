#' Rboretum Tree Trimmer
#'
#' This ape wrapper returns a phylo or multiPhylo object that has been pruned to include only specified taxa. Note: Function will STOP if any 'tree' is missing any 'taxa'.
#' @param tree phylo or multiPhylo object
#' @param taxa Character vector of desired tip labels to keep (or discard if remove=TRUE)
#' @param remove OPTIONAL: If TRUE, tip labels specified by 'taxa' are removed from all trees rather than retained [Default: FALSE, prune 'tree' down to 'taxa']
#' @return Pruned phylo or multiPhylo object
#' @export

treeTrimmer <- function(tree,taxa,remove){
  
  if(!Rboretum::isMultiPhylo(tree) & !Rboretum::isPhylo(tree)){
    stop("'tree' does not appear to be a valid phylo or multiPhylo object")
  }
  
  if(missing(remove)){
    remove <- FALSE
  } else if(!is.logical(remove)){
    remove <- FALSE
  }
  
  if(!remove){
    
    if(!checkTips(tree,taxa)){
      stop("Specified 'taxa' are missing from at least one tree.")
    }
    
    if(length(taxa) < 3){
      stop("Can't trim to fewer than three tips.")
    }
    
    if(Rboretum::isPhylo(tree)){
      return(ape::drop.tip(tree,tree$tip.label[-match(taxa, tree$tip.label)]))
    } 
    
    if(Rboretum::isMultiPhylo(tree)){
      
      if(!is.null(names(trees))){
        tree_names <- names(trees)
        namedTrees <- TRUE
      } else{ namedTrees <- FALSE }
      
      for(i in 1:length(tree)){
        tree[[i]] <- ape::drop.tip(tree[[i]],tree[[i]]$tip.label[-match(taxa, tree[[i]]$tip.label)])
      }
      
      if(namedTrees){ names(tree)  <- tree_names }
      
      return(tree)
    }
  }
  
  if(remove){
    
    if(Rboretum::isPhylo(tree)){
      
      keep_taxa <- tree$tip.label[!tree$tip.label %in% taxa]
      
      if(length(keep_taxa) < 3){
        stop("Can't trim to fewer than three tips.")
      } else{
          return(ape::drop.tip(tree,tree$tip.label[-match(keep_taxa, tree$tip.label)]))
      }
    }
    
    if(Rboretum::isMultiPhylo(tree)){
      
      if(!is.null(names(trees))){
        tree_names <- names(trees)
        namedTrees <- TRUE
      } else{namedTrees <- FALSE }
        
      for(i in 1:length(tree)){
        
        keep_taxa <- tree[[i]]$tip.label[!tree[[i]]$tip.label %in% taxa]

        if(length(keep_taxa) < 3){
          stop("Can't trim to fewer than three tips.")
        } else{
          tree[[i]] <- ape::drop.tip(tree[[i]],tree[[i]]$tip.label[-match(keep_taxa, tree[[i]]$tip.label)])
        } 
      }
        if(namedTrees){ names(tree)  <- tree_names }
        
        return(tree)
    }
  } 
}