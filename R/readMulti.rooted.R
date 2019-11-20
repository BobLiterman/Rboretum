#' Rboretum Multiphylo Rooted Tree Reader
#'
#' This function is an ape wrapper, and returns a multiPhylo object, with each tree rooted as specified by the user.
#' @param tree_paths Vector of paths to tree files that can be read by ape::read.tree() or ape::read.nexus()
#' @param root_taxa Character vector containing outgroup species IDs (Must be in all trees and always monophyletic)
#' @return A multiPhylo object, with each tree rooted at specified taxa
#' @export
#' @examples
#' myTrees <- c('/path/to/tree1','/path/to/tree2')
#' myRootTaxa <- c('Spp1','Spp2')
#' myMulti <- readMulti.rooted(myTrees,myRootTaxa)

readMulti.rooted <- function(tree_paths,root_taxa){
  if(length(tree_paths)<2){
    stop("'tree_paths' contains fewer than 2 items, use read.rooted()")
  } else{ 
    tree_count <- length(tree_paths)
    }

  trees <- c(ape::rtree(10),ape::rtree(10))

  for(i in 1:tree_count){
    if(has_error(read.rooted(tree_paths[[i]],root_taxa))){
      print(tree_paths[[i]])
      stop("The above file above cannot be read or rooted using read.rooted()")
    } else{
      trees <- c(trees,read.rooted(tree_paths[[i]],root_taxa))
    }
  }
  
  trees <- trees[3:2+tree_count]
}
