#' Rboretum Tree Splitter
#'
#' This function breaks down a rooted phylo into its respective set of splits
#' @param tree Rooted phylo object
#' @return Dataframe: Column 1: Semi-colon separated monophyletic clade; Column 2: 'Mirror Clade'; Column 3: Phylo Node ID (NA for root split)
#' @export

getTreeSplits_Worker <- function(tree){
  
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'getTreeSplits_Worker' requires a rooted tree.")
  }
  
  # Strip bootsrap values if present
  tree <- stripNodeLabels(tree)
  
  # Get species + subtree
  tree_species <- naturalsort(tree$tip.label)
  tree_subtree <- subtrees(tree)
  subtree_length <- length(tree_subtree)
  
  # Root detection for clade listing
  if(detectSingleRoot(tree)){
    tree_subtree <- subtrees(tree)[2:subtree_length] # Remove top subtree (whole tree)
  } else{
    tree_subtree <- subtrees(tree)[3:subtree_length] # Remove top two subtrees (whole tree + other side of root)
  }
  
  # Get clades, mirror clades, and node list
  tree_clades <- purrr::map(.x=tree_subtree,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()
  mirror_clades <- purrr::map(.x=tree_clades,.f=function(x){semiSorter(setdiff(tree_species,semiVector(x)))}) %>% unlist()
  node_list <- purrr::map(.x=tree_subtree,.f=function(x){x$node.label[1]}) %>% unlist()
  
  split_df <- tibble("Clade"=as.character(tree_clades),"Mirror_Clade"=as.character(mirror_clades),"Split_Node"=as.integer(node_list)) %>%
    rowwise() %>%
    mutate(Split_Node = ifelse(ape::is.monophyletic(tree,semiVector(Clade)) & ape::is.monophyletic(tree,semiVector(Mirror_Clade)),NA,Split_Node)) %>%
    ungroup()

  return(split_df)
}