#' Rboretum Identical Topology Checker
#'
#' This function returns TRUE if all trees in multiPhylo have the same topology after pruning to identical species lists
#' @param trees multiPhylo object
#' @return TRUE if all trees have the same topology after pruning to identical species lists; else, FALSE
#' @export
#' @examples
#' myTrees <- c(Tree_1,Tree_2,Tree_3)
#' same.topology(trees)
#'
same.topology <- function(trees){
  if(!Rboretum::is.multiPhylo(trees)){
    stop("'trees' does not appear to be a valid multiPhylo object with 2+ trees")
  } else if(has_error(Rboretum::check.shared(trees))){
    stop("Trees do not appear to share three species in common.")
  } else if(!Rboretum::check.shared(trees)){
    stop("Trees do no share three species in common.")
  }
  
  if(!Rboretum::same.taxa(trees)){
    shared_speces <- Rboretum::get.shared(trees)
    trees <- Rboretum::treeTrimmer(trees,shared_speces)
  }
  
  tree_count <- length(trees)
  tree_list <- c()
  
  for(i in 1:(tree_count-1)){
    for(j in (i+1):tree_count){
      tree_list <- c(tree_list,ape::all.equal.phylo(trees[[i]],trees[[j]],use.edge.length = FALSE))
    }
  }
  
  if(all(tree_list)){
    return(TRUE)
  } else{ return(FALSE) }
}