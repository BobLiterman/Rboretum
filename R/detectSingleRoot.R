#' Rboreturm Single Root Checker
#'
#' This function returns TRUE if 'tree' is rooted on a single taxon, or FALSE if not
#' @param tree Tree to check. Must be a rooted phylo object
#' @param return_root OPTIONAL: If TRUE and tree is rooted on a single taxon, return that tip label; else return FALSE [Default: FALSE, don't return root label] 
#' @return TRUE if 'tree' is rooted on a single taxon; else, FALSE
#' @export
#' 

detectSingleRoot <- function(tree,return_root){
  
  # Check valid rooted phylo
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'detect_single_root' requires a rooted phylo object")
  }
  
  # Check return_root
  if(missing(return_root)){
    return_root <- FALSE
  } else if(!is.logical(return_root)){
    return_root <- FALSE
  } else if(length(return_root)!=1){
    return_root <- FALSE
  }
  
  tree_taxa <- tree$tip.label
  subtree_taxa <- purrr::map(.x=subtrees(tree)[2:length(ape::subtrees(tree))],.f=function(x){x$tip.label}) %>% unlist() %>% unique()
  
  # Single root trees will have one fewer tip labels in the subtrees relative to the tree
  if(length(tree_taxa) == length(subtree_taxa)){
    return(FALSE)
  }
  
  if(return_root){
    if(length(tree_taxa) - length(subtree_taxa) == 1){
      return(setdiff(tree_taxa,subtree_taxa))
    } else{
      warning("UNEXPECTED RESULT FROM detectSingleRoot")
      return(FALSE)
    }
  } else{
    if(length(tree_taxa) - length(subtree_taxa) == 1){
      return(TRUE)
    } else{
      warning("UNEXPECTED RESULT FROM detectSingleRoot")
      return(FALSE)
    }
  }
}