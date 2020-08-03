#' Rboreturm Tip Checker (X)
#'
#' This function returns TRUE if all 'taxa' are present in 'tree' (or a multiPhylo of trees); FALSE otherwise. If multiPhylo, all comparisons are made on the common set of taxa.
#' @param tree Tree(s) to check. Options include:
#' \itemize{
#'   \item A phylo object
#'   \item A multiPhylo object where all trees share 3+ taxa
#' }
#' @param taxa Character vector containing taxa to query
#' @param check_mono OPTIONAL: If TRUE, also query whether 'taxa' are monophyletic in tree(s) [Default: FALSE]
#' @param check_root OPTIONAL: If TRUE, also query whether 'taxa' are monophyletic and make up one of the two root splits in tree(s) [Default: FALSE]
#' @return TRUE if all 'taxa' in 'tree'; else, FALSE (barring other selected options)
#' @examples 
#' # Check if all labels in 'myTaxa' are present in phylo object 'myTree'
#' checkTips(myTree,myTaxa)
#' # Check if all labels in 'myTaxa' are present in all trees in multiPhylo object 'myTrees'
#' checkTips(myTrees,myTaxa)
#' # Check if clade 'myTaxa' are monophyletic in phylo object 'myTree'
#' checkTips(myTree,myTaxa,check_mono=TRUE)
#' # Check if clade 'myTaxa' are one of the two root splits in phylo object 'myTree' 
#' checkTips(myTree,myTaxa,check_root=TRUE)
#' @export

checkTips <- function(tree,taxa,check_mono,check_root){
  
  # Check if 'tree' is valid
  if(missing(tree)){
    stop("'checkTips' requires a 'tree' argument.")
  } else if(!Rboretum::isMultiPhylo(tree,check_three_taxa = TRUE) & !Rboretum::isPhylo(tree)){
    stop("'tree' does not appear be a valid phylo or multiPhylo object where all trees share 3+ taxa.")
  }
  
  # Check if 'taxa' is a character vector
  if(missing(taxa)){
    stop("No 'taxa' to check")
  } else if(!is.character(taxa)){
    stop("'taxa' should be a character vector of tip labels to check.")
  }
  
  # Check if 'taxa' is a semi-colon separated character of length == 1. If so, vectorize
  if(length(taxa) == 1 & Rboretum::semiChecker(taxa)){
    taxa <- semiVector(taxa)
  }
  
  # Check monophly?
  if(missing(check_mono)){
    check_mono <- FALSE
  } else if(!is.logical(check_mono)){
    check_mono <- FALSE
  }
  
  # Check root?
  if(missing(check_root)){
    check_root <- FALSE
  } else if(!is.logical(check_root)){
    check_root <- FALSE
  } else if(check_root & !Rboretum::isPhylo(tree,check_root=TRUE) & !Rboretum::isMultiPhylo(tree,check_root=TRUE)){
    stop("'check_root' is enabled but 'tree' contains unrooted trees.") # Can't check rootness of unrooted trees
  }

  if(Rboretum::isPhylo(tree)){ # If 'tree' is a single tree...

    if(all(taxa %in% tree$tip.label)){ # Check if all 'taxa' are in 'tree'
      
      tree_taxa <- tree$tip.label
      mirror_taxa <- sort(dplyr::setdiff(tree_taxa, taxa))
      
      if(check_mono & !ape::is.monophyletic(tree,taxa)){
        return(FALSE) # 'taxa' in 'tree', but not monophyletic --> FALSE
      } else if(check_root & !(ape::is.monophyletic(tree,taxa) & ape::is.monophyletic(tree,mirror_taxa))){
        return(FALSE) # 'taxa' in 'tree', but is not one of the two root splits --> FALSE
      } else{
        return(TRUE) # 'taxa' in  tree and neither check failed --> TRUE
      }
    } else{
      return(FALSE) # 'taxa' not in 'tree' -> FALSE
    }
  }
  
  if(Rboretum::isMultiPhylo(tree)){
    
    # Trim to common taxa set prior to testing
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa=TRUE)){
      tree <- Rboretum::treeTrimmer(tree) 
    }

    if(purrr::map(.x = tree,.f = function(x){all(taxa %in% x$tip.label)}) %>% unlist() %>% all()){ # Check if all 'taxa' are in all trees in 'tree'
      
      tree_taxa <- Rboretum::getSharedTaxa(tree)
      mirror_taxa <- sort(dplyr::setdiff(tree_taxa, taxa))
      
      if(check_mono & !purrr::map(.x = tree, .f=function(x){ape::is.monophyletic(x,taxa)}) %>% unlist() %>% all()){
        return(FALSE) # 'taxa' in 'tree', but not monophyletic in all trees --> FALSE
      } else if(check_root & !purrr::map(.x = tree, .f=function(x){ape::is.monophyletic(x,taxa) & ape::is.monophyletic(x,mirror_taxa)}) %>% unlist() %>% all()){
        return(FALSE) # 'taxa' in 'tree', but is not the root of all trees --> FALSE
      } else{
        return(TRUE) # 'taxa' in  tree and neither check failed --> TRUE
      }
    } else{
      return(FALSE) # 'taxa' not present in all trees in 'tree' -> FALSE
    }
  }
}