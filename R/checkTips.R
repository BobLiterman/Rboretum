#' Rboreturm Tip Checker
#'
#' This function returns TRUE if all 'taxa' are present in 'tree' (or a multiPhylo of trees); FALSE otherwise. If multiPhylo, monophyly + root comparisons are made on the common set of taxa.
#' @param tree Tree(s) to check. Options include:
#' \itemize{
#'   \item A phylo object
#'   \item A multiPhylo object
#' }
#' @param taxa Character vector containing taxa to query
#' @param check_mono OPTIONAL: If TRUE, also query whether 'taxa' are monophyletic in tree(s) [Default: FALSE]
#' @param check_root OPTIONAL: If TRUE, also query whether 'taxa' are monophyletic and make up one of the two root splits in tree(s) [Default: FALSE]
#' @return TRUE if all 'taxa' in 'tree'; else, FALSE
#' @examples 
#' checkTips(myTree,myTaxa) # Check if all labels in 'myTaxa' are present in phylo object 'myTree'
#' checkTips(myTrees,myTaxa) # Check if all labels in 'myTaxa' are present in all trees in multiPhylo object 'myTrees'
#' checkTips(myTree,myTaxa,check_mono=TRUE) # Check if clade 'myTaxa' are monophyletic in phylo object 'myTree'
#' checkTips(myTree,myTaxa,check_root=TRUE) # Check if clade 'myTaxa' are one of the two root splits in phylo object 'myTree'
#' @export

checkTips <- function(tree,taxa,check_mono,check_root){
  
  # Check if 'tree' is valid
  if(!Rboretum::isPhylo(tree) & !Rboretum::isMultiPhylo(tree)){
    stop("'tree' does not appear be a valid phylo or multiPhylo object.")
  }
  
  # Check if taxa is a character vector
  if(missing(taxa)){
    stop("No 'taxa' to check")
  } else if(!is.character(taxa)){
    stop("'taxa' should be a character vector of tip labels to check.")
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
    
    if(!Rboretum::isMultiPhylo(tree,check_three_taxa=TRUE)){
      if(check_mono | check_root){
        stop("Trees in 'tree' do not share at least three taxa, and cannot be assessed for monophyly or root.")
      }
    }

    if(purrr::map(.x = tree,.f = function(x){all(taxa %in% x$tip.label)}) %>% unlist() %>% all()){ # Check if all 'taxa' are in all trees in 'tree'
      
      tree_taxa <- Rboretum::getSharedTaxa(tree)
      
      if(!Rboretum::isMultiphylo(tree,check_all_taxa=TRUE)){
        tree <- Rboretum::treeTrimmer(tree,tree_taxa) # Trim to common taxa set prior to assessing monophyly or root
      }
      
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