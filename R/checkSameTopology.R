#' Rboretum Identical Topology Checker
#'
#' This function returns TRUE if all trees in multiPhylo have the same topology after pruning to identical species lists
#' @param trees Rooted, multiPhylo object where all trees share at least three taxa
#' @param check_any OPTIONAL: If TRUE, check if ANY (rather than all) trees share a topology [Default: FALSE, check for unanimous topology]
#' @return TRUE if all trees have the same topology after pruning to identical species lists; else, FALSE
#' @examples 
#' checkSameTopology(myTrees) # Check if all trees in multiPhylo 'myTrees' share the same topology
#' checkSameTopology(myTrees,check_any=TRUE) # Check if any trees in multiPhylo 'myTrees' share a topology
#' @export

checkSameTopology <- function(trees,check_any){
  
  # Check if 'trees' is a valid multiPhylo where all trees share >= 3 taxa
  if(!Rboretum::isMultiPhylo(trees,check_rooted = TRUE,check_three_taxa = TRUE)){
    stop("'trees' does not appear to be a valid multiPhylo object where all trees are rooted and share at least three taxa.")
  }
  
  # Check for any duplicated topologies (versus all)?
  if(missing(check_any)){
    check_any <- FALSE
  } else if(!is.logical(check_any)){
    check_any <- FALSE
  }
  
  tree_taxa <- Rboretum::getSharedTaxa(trees)
  tree_count <- length(trees)
  
  # Trim to common taxa set if necessary
  if(!Rboretum::isMultiPhylo(trees,check_all_taxa = TRUE)){
    trees <- treeTrimmer(trees,tree_taxa)
  }

  # Compare all tree topologies
  top_check <- c()

  for(i in 1:(tree_count-1)){
    for(j in 2:tree_count){
      top_check <- c(top_check,ape::all.equal.phylo(trees[[i]],trees[[j]],use.edge.length = FALSE))
    }
  }

  # check_any= TRUE will return TRUE if any two trees share a topology
  if(check_any){
    if(any(top_check)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  } else{ # check_any= FALSE [Default] will return TRUE if all trees share the same topology
    if(all(top_check)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}