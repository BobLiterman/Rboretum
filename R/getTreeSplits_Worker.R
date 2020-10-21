#' Rboretum Tree Splitter
#'
#' This function breaks down a rooted phylo into its respective set of splits
#' @param tree Rooted phylo object
#' @return Dataframe: Column 1: Semi-colon separated monophyletic clade; Column 2: 'Mirror Clade'; Column 3: Phylo Node ID
#' @export

getTreeSplits_Worker <- function(tree){
  
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'getTreeSplits_Worker' requires a rooted tree.")
  }
  
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
  mirror_clades <- purrr::map(.x=subtree_clades,.f=function(x){semiSorter(tree_tips[!tree_tips %in% semiVector(x)])}) %>% unlist()
  
  # Get root clades and remove from larger list
  root_clades <- getTreeRoot(tree)
  root_node <- getMRCA(tree,tree$tip.label)
  
  # Get non-root clades
  non_root_clades <- subtree_clades[!subtree_clades %in% root_clades]
  non_root_mirror <- mirror_clades[!subtree_clades %in% root_clades]
  
  # Get node IDs
  non_root_mrca <- purrr::map(.x=non_root_clades,.f=function(x){getMRCA(tree,semiVector(x))}) %>% unlist()
  
  tree_split_df <- tibble(
    Clade = non_root_clades,
    Mirror_Clade = non_root_mirror,
    Split_Node = non_root_mrca,
    Root = FALSE) %>%
    add_row(Clade = root_clades[1],
            Mirror_Clade = root_clades[2],
            Split_Node = root_node,
            Root=TRUE)
  
  return(tree_split_df)
}