#' Rboretum Root Fetcher
#'
#' This function takes a rooted phylo object, and returns the representative root split as two semicolon-delimited characters
#' @param tree A rooted phylo object
#' @return A two-element character vector containing semicolon-delimited root clades
#' @export

getRoot <- function(tree){
  
  # Ensure tree is passed
  if(missing(tree)){ stop("'getRoot' requires a phylo object passed via the 'tree' argument...")} 

  # Ensure tree is rooted
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE)){ stop("'getRoot' requires a rooted phylo object passed via the 'tree' argument...")} 
  
  # Get tip labels and count
  tree_tips <- naturalsort(tree$tip.label)
  tip_count <- length(tree_tips)
  
  # Get subtrees (except for subtree with all taxa)
  all_subtrees <- ape::subtrees(tree)
  subtree_taxa_counts <- purrr::map(.x=all_subtrees,.f=function(x){length(x$tip.label)}) %>% unlist()
  clade_subtrees <- all_subtrees[subtree_taxa_counts != tip_count]
  
  # Get all taxa + clades included in clade subtrees
  subtree_taxa <- purrr::map(.x=clade_subtrees,.f=function(x){x$tip.label}) %>% unlist() %>% unique() %>% naturalsort()
  subtree_taxa_count <- length(subtree_taxa)
  subtree_clades <- purrr::map(.x=clade_subtrees,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()
  
  # The clade subtrees of trees rooted on a single taxon are missing that taxon (tip_count - subtree_taxa_count = 1)
  # The clade subtrees of trees rooted on a 2+ taxon should contain all tree tips (tip_count - subtree_taxa_count = 0)
  
  if(tip_count - subtree_taxa_count > 1){ stop("Clade subtrees are missing more than 1 tip relative to all tip labels...Unexpected error in ape::subtrees()")} 

  if(tip_count - subtree_taxa_count == 1){
    return_root <- c(tree_tips[!tree_tips %in% subtree_taxa],semiSorter(tree_tips[tree_tips %in% subtree_taxa]))
  } else if(tip_count - length(subtree_taxa) == 0){
    mirror_clades <- purrr::map(.x=subtree_clades,.f=function(x){semiSorter(tree_tips[!tree_tips %in% semiVector(x)])}) %>% unlist()
    root_split <- mirror_clades[mirror_clades %in% subtree_clades]
    
    # Return root split as two-element vector, where element 1 has fewer taxa than element 2 (if same, return sorted)
    if(str_count(root_split[1],";") < str_count(root_split[2],";")){
      return_root <- c(root_split[1],root_split[2])
    } else if(str_count(root_split[1],";") > str_count(root_split[2],";")){
      return_root <- c(root_split[2],root_split[1])
    } else if(str_count(root_split[1],";") == str_count(root_split[2],";")){
      return_root <- naturalsort(root_split)
    }
  }
  return(return_root)
}