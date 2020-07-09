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
#' @param keep_bs OPTIONAL: If TRUE, do not remove bootstrap labels when taxa are trimmed off of tree [Default: FALSE, if tree has taxa removed, also remove node labels]
#' @return Pruned phylo or multiPhylo object
#' @examples 
#' # myTree is a phylo object with 4 tips
#' desiredTaxa <- c('Species_1','Species_2','Species_3')
#' undesiredTaxa <- 'Species_4'
#' trimmedTree <- treeTrimmer(myTree,desiredTaxa) # Trim myTree down to taxa in 'desiredTaxa'
#' trimmedTree <- treeTrimmer(myTree,undesiredTaxa,remove=TRUE) # Remove 'undesiredTaxa' from myTree
#' @export

treeTrimmer <- function(tree,taxa,remove,keep_bs){
  
  # Check if 'tree' is a valid object
  if(!Rboretum::isMultiPhylo(tree) & !Rboretum::isPhylo(tree)){
    stop("'tree' does not appear to be a valid phylo or multiPhylo object")
  }
  
  # Get tree taxa and names for multiPhylo
  if(Rboretum::isPhylo(tree)){
    tree_taxa <- naturalsort(tree$tip.label)
    common_taxa <- tree_taxa
  } else{
    tree_taxa <- purrr::map(.x=tree,.f=function(x){x$tip.label}) %>% unlist() %>% unique() %>% naturalsort() # Get all unique tip labels among 'trees'
    common_taxa <- Rboretum::getSharedTaxa(tree)
    
    # If trees have names, fetch them
    if(Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree_names <- names(tree)
    } else{ # Otherwise,create dummy names
      tree <- treeNamer(tree)
      tree_names <- names(tree)
    }
  }
  
  # Remove or retain 'taxa'?
  if(missing(remove)){
    remove <- FALSE
  } else if(!is.logical(remove)){
    remove <- FALSE
  } else if(length(remove)!=1){
    remove <- FALSE
  }
  
  # Strip bootstraps after trim?
  if(missing(keep_bs)){
    keep_bs <- FALSE
  } else if(!is.logical(keep_bs)){
    keep_bs <- FALSE
  } else if(length(keep_bs)!=1){
    keep_bs <- FALSE
  }
  
  # Check taxa
  if(Rboretum::isPhylo(tree)){
    if(missing(taxa)){
      stop("'treeTrimmer' requires a character vector of taxa to remove or retain.")
    } else if(!is.character(taxa)){
      stop("'treeTrimmer' requires a character vector of taxa to remove or retain.")
    }
  } else{
    
    # If a multiPhylo is provided without taxa, treeTrimmer will trim all trees down to the common set of taxa
    if(missing(taxa)){
      taxa <- common_taxa
      remove <- FALSE
    } else if(!is.character(taxa)){
      stop("'treeTrimmer' requires a character vector of taxa to remove or retain.")
    }
  }

  # Ensure 'taxa' are in 'tree' if remove=FALSE
  if(!remove & !all(taxa %in% common_taxa)){
    stop("Specified 'taxa' are missing from at least one tree.")
  }
  
  # Ensure 2+ 'taxa' if remove=FALSE
  if(!remove & length(taxa)<2){
    stop("Can't trim to fewer than two tips.")
  }
  
  # Check that all trees will have 2+ tip labels after pruning if remove=TRUE
  if(remove){
    if(Rboretum::isPhylo(tree)){
      taxa_to_keep <- tree_taxa[!tree_taxa %in% taxa]
      if(length(taxa_to_keep)<2){
        stop("Can't trim to fewer than two tips.")
      }
    } else{
      if(any(purrr::map(.x = tree,.f = function(x){length(x$tip.label[!x$tip.label %in% taxa])}) %>% unlist() < 2)){
        stop("Can't trim to fewer than two tips.")
      }
    }
  }
  
  # Set taxa to remove
  if(remove){
    taxa_to_remove <- taxa
  } else{
    if(Rboretum::isPhylo(tree)){
      taxa_to_remove <- tree_taxa[!tree_taxa %in% taxa]
    } else{
      taxa_to_remove <- purrr::map(.x=tree,.f=function(x){x$tip.label[!x$tip.label %in% taxa]})
    }
  }
  
  # Trim trees and return
  if(Rboretum::isPhylo(tree)){
    
    # If the phylo already contains the desired taxa, return unchanged
    if(!any(taxa_to_remove %in% tree_taxa)){
      return(tree)
    } else{ # If the phylo needs to be pruned, return tree without bootstrap values
      
      if(!keep_bs){
        tree$node.label <- NULL
      }
      return(ape::drop.tip(tree,taxa_to_remove))
    }
  } else{ # Process multiPhylo
    
    return_tree <- tree
    tree_count <- length(return_tree)
    
    # Find trees that require pruning (i.e. that contain taxa)
    for(i in 1:tree_count){
      
      temp_tree <- tree[[i]]
      
      # If the phylo already contains the desired taxa, return unchanged
      if(!any(taxa_to_remove[[i]] %in% temp_tree$tip.label)){
        return_tree[[i]] <- temp_tree
      } else{ # If the phylo needs to be pruned, return tree without bootstrap values
        
        if(!keep_bs){
          temp_tree$node.label <- NULL
        }
        
        return_tree[[i]] <- ape::drop.tip(temp_tree,taxa_to_remove[[i]])
      }
    }
    class(return_tree) <- "multiPhylo"
    names(return_tree) <- tree_names
    return(return_tree)
  }
}