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
    
    if(length(taxa) < 2){
      stop("Can't trim to fewer than two tips.")
    }
    
    if(Rboretum::isPhylo(tree)){
      return(ape::drop.tip(tree,tree$tip.label[-match(taxa, tree$tip.label)]))
    } 
    
    if(Rboretum::isMultiPhylo(tree)){
      
      if(!is.null(names(tree))){
        tree_names <- names(tree)
        namedTrees <- TRUE
      } else{ namedTrees <- FALSE }
      
      tree <- purrr::map(.x = tree,.f = function(x){ape::drop.tip(x,x$tip.label[-match(taxa, x$tip.label)])})
      class(tree) <- "multiPhylo"
      
      if(namedTrees){ names(tree)  <- tree_names }
      
      return(tree)
    }
  }
  
  if(remove){
    
    if(Rboretum::isPhylo(tree)){
      
      keep_taxa <- tree$tip.label[!tree$tip.label %in% taxa]
      
      if(length(keep_taxa) < 2){
        stop("Can't trim to fewer than two tips.")
      }
      
      return(ape::drop.tip(tree,tree$tip.label[match(taxa, tree$tip.label)]))
    } else if(Rboretum::isMultiPhylo(tree)){
      
      if(!is.null(names(tree))){
        tree_names <- names(tree)
        namedTrees <- TRUE
      } else{namedTrees <- FALSE }
        
      if(any(purrr::map(.x = tree,.f = function(x){length(x$tip.label[!x$tip.label %in% taxa])}) %>% unlist() < 2)){
        stop("Can't trim to fewer than two tips.")
      } else{
        tree <- purrr::map(.x = tree, .f = function(x){ape::drop.tip(x,x$tip.label[match(taxa, x$tip.label)])})
        class(tree) <- "multiPhylo"
        if(namedTrees){ names(tree)  <- tree_names }
        return(tree)
      }
    }
  } 
}