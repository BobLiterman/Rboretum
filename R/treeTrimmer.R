#' Rboretum Tree Trimmer
#'
#' This function returns a phylo or multiPhylo object that has been pruned to include only specified taxa. Note: Function will STOP if any 'tree' is missing any 'taxa' or if trimming will result in <2 final taxa
#' @param tree Tree(s) to trim. Options include:
#' \itemize{
#'   \item A phylo object
#'   \item A multiPhylo object
#' }
#' @param taxa Character vector of desired tip labels to keep (or discard if remove=TRUE)
#' @param remove OPTIONAL: If TRUE, tip labels specified by 'taxa' are removed from all trees rather than retained [Default: FALSE, prune 'tree' down to 'taxa']
#' @return Pruned phylo or multiPhylo object
#' @examples 
#' # myTree is a phylo object with 4 tips
#' desiredTaxa <- c('Species_1','Species_2','Species_3')
#' undesiredTaxa <- 'Species_4'
#' trimmedTree <- treeTrimmer(myTree,desiredTaxa) # Trim myTree down to taxa in 'desiredTaxa'
#' trimmedTree <- treeTrimmer(myTree,undesiredTaxa,remove=TRUE) # Remove 'undesiredTaxa' from myTree
#' @export

treeTrimmer <- function(tree,taxa,remove){
  
  # Check if 'tree' is a valid object
  if(!Rboretum::isMultiPhylo(tree) & !Rboretum::isPhylo(tree)){
    stop("'tree' does not appear to be a valid phylo or multiPhylo object")
  }
  
  # Remove 'taxa' or retain 'taxa'?
  if(missing(remove)){
    remove <- FALSE
  } else if(!is.logical(remove)){
    remove <- FALSE
  }
  
  if(!remove){ # If 'taxa' is the desired list of species to keep...
    
    if(!checkTips(tree,taxa)){ # Ensure 'taxa' exist in 'tree' or all trees in 'tree'
      stop("Specified 'taxa' are missing from at least one tree.")
    } else if(length(taxa) < 2){ # Ensure 'taxa' includes 2+ tip labels
      stop("Can't trim to fewer than two tips.")
    }
    
    # Prune single tree down to just those tip labels in 'taxa'
    if(Rboretum::isPhylo(tree)){
      return(ape::drop.tip(tree,tree$tip.label[-match(taxa, tree$tip.label)]))
    } 
    
    if(Rboretum::isMultiPhylo(tree)){ # If a multiPhylo is provided...
      
      # Grab tree names, if named
      if(!is.null(names(tree))){
        tree_names <- names(tree)
        namedTrees <- TRUE
      } else{ namedTrees <- FALSE }
      
      # Prune all trees down to just those tip labels in 'taxa' and rename if necessary.
      tree <- purrr::map(.x = tree,.f = function(x){ape::drop.tip(x,x$tip.label[-match(taxa, x$tip.label)])})
      class(tree) <- "multiPhylo"
      if(namedTrees){ names(tree)  <- tree_names }
      
      return(tree)
    }
  } else{ # If 'taxa' is the desired list of species to remove...
    
    if(Rboretum::isPhylo(tree)){ # If a single tree is provided...
      
      keep_taxa <- tree$tip.label[!tree$tip.label %in% taxa] # Get tip labels to keep (those not in 'taxa')
      
      if(length(keep_taxa) < 2){ # Ensure final tree includes 2+ tip labels
        stop("Can't trim to fewer than two tips.")
      }
      
      return(ape::drop.tip(tree,tree$tip.label[match(taxa, tree$tip.label)]))
      
    } else if(Rboretum::isMultiPhylo(tree)){ # If a multiPhylo is provided...
      
      # Grab tree names, if named
      if(!is.null(names(tree))){
        tree_names <- names(tree)
        namedTrees <- TRUE
      } else{namedTrees <- FALSE }
      
      # Check that all trees will  have 2+ tip labels after pruning  
      if(any(purrr::map(.x = tree,.f = function(x){length(x$tip.label[!x$tip.label %in% taxa])}) %>% unlist() < 2)){
        stop("Can't trim to fewer than two tips.")
      } else{
        
        tree <- purrr::map(.x = tree, .f = function(x){ape::drop.tip(x,x$tip.label[match(taxa, x$tip.label)])}) # Remove 'taxa' from each tree in 'tree'
        class(tree) <- "multiPhylo"
        if(namedTrees){ names(tree)  <- tree_names } # Append names if necessary
        return(tree)
      }
    }
  } 
}