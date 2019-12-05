#' Rboretum Rooted Tree Reader Worker
#'
#' This function reads in and roots a phylo object
#' @param to_root_worker Path to tree file
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree and monophyletic)
#' @return A phylo object, rooted at specified taxa
#' @export

readRooted_Worker <- function(to_root_worker,root_taxa){
  
  # Ensure that a path and root taxa are provided as character vectors
  if(missing(to_root_worker)){
    return(NA) # No tree file or directories indicated with 'to_root_worker'
  } else if(!is.character(to_root_worker)){
    return(NA) # 'to_root_worker' should be a character path to a tree file.
  } else if(missing(root_taxa)){
    return(NA) # No root taxa provided
  } else if(!is.character(root_taxa)){
    return(NA) # 'root_taxa' should be a character vector of tip labels
  }

  if(file.exists(to_root_worker) & !dir.exists(to_root_worker)){
    
    # Read in tree and fetch species
    if(!has_error(silent=TRUE,expr=ape::read.tree(to_root_worker))){
      raw_tree <- ape::read.tree(to_root_worker)
      tree_species <- raw_tree$tip.label
    } else if(!has_error(silent=TRUE,expr=ape::read.nexus(to_root_worker))){
      raw_tree <- ape::read.nexus(to_root_worker)
      tree_species <- raw_tree$tip.label
    } else{ 
      return(NA) # 'to_root_worker' does not point to file that can be read in by ape::read.tree() or ape::read.nexus()
    }
    
    # Ensure at least one root taxon is present in the tree
    if(!any(root_taxa %in% tree_species)){
      return(NA) # 'root_taxa' completely absent from tree
    } else{
      root_taxa <- root_taxa[root_taxa %in% tree_species]
    }
    
    mirror_clade <- tree_species[!tree_species %in% root_taxa]
    
    # If tree is already rooted, check root and return, or unroot
    if(ape::is.rooted(raw_tree)){
      if(Rboretum::checkTips(raw_tree,root_taxa,check_root=TRUE)){
        return(raw_tree)
      } else{
        raw_tree <- ape::unroot.phylo(raw_tree)
      }
    }
    
    # Ensure root species in tree, and are monophyletic for rooting
    if(!Rboretum::checkTips(raw_tree,root_taxa,check_mono=TRUE)){
      return(NA) # 'root_taxa' are not monophyletic in 'to_root_worker', and cannot be used as a root
    } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))){
      return(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = TRUE))
    } else if(!has_error(silent=TRUE,expr=ape::root.phylo(raw_tree,outgroup = mirror_clade,edgelabel = TRUE,resolve.root = TRUE))){
      return(ape::root.phylo(raw_tree,outgroup = mirror_clade,edgelabel = TRUE,resolve.root = TRUE))
    } else{ 
      return(NA) # Ape cannot root tree on these taxa
    }
  
  } else{
    return(NA) # 'to_root_worker' does not point to a valid file
  }
}