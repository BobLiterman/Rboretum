#' Get Split Bootstrap
#'
#' This function returns the boostrap support for the MRCA of 'taxa' in 'tree'
#' @param tree Phylo object
#' @param taxa Vector containing taxa to check
#' @return Bootstrap support label for MRCA of 'taxa' in 'tree'
#' @export
#' @examples
#' getSplitBootstrap(tree,taxa)
#'

getSplitBootstrap <- function(tree,taxa){

  # Ensure speecies are in the tree
  if(has_error(Rboretum::checkMonophyletic(tree,taxa)) | !Rboretum::checkMonophyletic(tree,taxa)){
    stop('Taxa are not monophyletic, or not in the tree.')
  }

  # If specified taxa are a root split, find boostrap support for one of the two associated root nodes (only one should have a value if rooted from an unrooted tree)
  if(Rboretum::checkRoot(tree,taxa)){
    tree_species <- sort(tree$tip.label)

    root_1 <- sort(taxa)
    root_tree_1 <- Rboretum::getTrimmedTree(tree,root_1)
    root_1_BS <- root_tree_1$node.label[1]

    root_2 <- sort(dplyr::setdiff(tree_species, taxa))
    root_tree_2 <- Rboretum::getTrimmedTree(tree,root_2)
    root_2_BS <- root_tree_2$node.label[1]

    if(!is.na(as.numeric(root_1_BS))){
      return(as.numeric(root_1_BS))
    }

    else if(!is.na(as.numeric(root_2_BS))){
      return(as.numeric(root_2_BS))
    }

    else{
      return(NA)
    }
  }

  # If not a root node, return boostrap node label
  else{
    pruned.tree <- Rboretum::getTrimmedTree(tree,taxa)
    clade_bs <- pruned.tree$node.label[1]
    if(!is.na(as.numeric(clade_bs))){
      return(as.numeric(clade_bs))
    }
    else{
      return(NA)
    }
  }
}
