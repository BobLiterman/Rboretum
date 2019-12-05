#' Rboretum Identical Topology Checker
#'
#' This function returns TRUE if all trees in multiPhylo have the same topology after pruning to identical species lists
#' @param trees multiPhylo object
#' @return TRUE if all trees have the same topology after pruning to identical species lists; else, FALSE
#' @export

checkSameTopology <- function(trees){
  
  if(!Rboretum::isMultiPhylo(trees,check_rooted = TRUE,check_three_taxa = TRUE)){
    stop("'trees' does not appear to be a valid multiPhylo object where all trees are rooted and share at least three taxa.")
  }
  
  if(!Rboretum::isMultiPhylo(trees,check_all_taxa = TRUE)){
    tree_taxa <- Rboretum::getSharedTaxa(trees)
    trees <- treeTrimmer(trees,tree_taxa)
  }
  
  tree_count <- length(trees)
  
  top_check <- c()
  
  for(i in 1:(tree_count-1)){
    for(j in 2:tree_count){
      top_check <- c(top_check,ape:all.equal.phylo(trees[[i]],trees[[j]],use.edge.length = FALSE))
    }
  }
  
  if(all(top_check)){
    return(TRUE)
  } else{
    return(FALSE)
  }
  
}