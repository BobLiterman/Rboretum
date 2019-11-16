#' Rboretum Rooted Tree Reader
#'
#' This function is an ape wrapper, and returns a tree rooted as specified by the user.
#' @usage treePath <- '/path/to/my/tree.nwk'
#' rootTaxa <- c('Root_species_1','Root_species_2')
#' rootedTree <- read.rooted(tree_path,root_taxa)
#' @param tree_path Path to tree file that can be read by ape::read.tree() or ape::read.nexus()
#' @param root_taxa Character vector containing outgroup species IDs (Must be in tree and monophyletic)
#' @param resolve_root OPTIONAL: Set ape::resolve.root [Default: TRUE]
#' @return A phylo object, rooted at specified taxa
#' @export
#' @examples
#' myRootedTetrapods <- read.rooted('/path/to/tetrapod_tree.nwk',c('Fish_1','Fish_2'))

read.rooted <- function(tree_path,root_taxa,resolve_root){
  
  if(missing(resolve_root)){
    resolve_root <- TRUE
  } else if(!is.logical(resolve_root)){
    resolve_root <- TRUE
  }

  # Read in tree and fetch species
  if(!has_error(ape::read.tree(tree_path))){
    raw_tree <- ape::read.tree(tree_path)
    tree_species <- raw_tree$tip.label
  } else if(!has_error(ape::read.nexus(tree_path))){
    raw_tree <- ape::read.nexus(tree_path)
    tree_species <- raw_tree$tip.label
  } else{ stop("Path does not point to file that can be read in by ape::read.tree() or ape::read.nexus()") }

  # Ensure root species in tree, and are monophyletic for rooting
  if(!(all(root_taxa %in% tree_species))){
    stop("Root taxa not found in tree.")
  }

  if(ape::is.rooted(raw_tree)){
    if(Rboretum::check.root(raw_tree,root_taxa)){
      return(raw_tree)
    } else if(!has_error(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = resolve_root))){
      return(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = resolve_root))
    } else{ stop("Ape cannot root tree on these taxa.") }
  } else if(!has_error(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = resolve_root))){
    return(ape::root.phylo(raw_tree,outgroup = root_taxa,edgelabel = TRUE,resolve.root = resolve_root))
  } else{ stop("Ape cannot root tree on these taxa.") }
}
