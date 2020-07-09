#' Rboretum Bootstrap Combiner
#'
#' This function takes a multiPhylo of trees that all share a common topolgy, and returns a single Phylo object with combined node labels (e.g. bootstrap values)
#' @param trees A multiPhylo object where all trees share all taxa and a common topology
#' @return A phylo object with concatenated bootstrap values/node labels
#' @export
#' 
combineBootstraps <- function(trees){
  
  # Assess validity of multiPhylo
  if(missing(trees)){
    stop("'combineBootstraps' requires a rooted multiPhylo where all trees share all taxa, and share a common topology.")
  } else if(!Rboretum::isMultiPhylo(trees,check_all_taxa=TRUE,check_all_equal=TRUE,check_rooted = TRUE)){
    stop("'combineBootstraps' requires a rooted multiPhylo where all trees share all taxa, and share a common topology.")
  }
  
  tree_count <- length(trees)  
  
  # Name trees if necessary
  if(!Rboretum::isMultiPhylo(trees,check_named = TRUE)){
    trees <- Rboretum::treeNamer(trees)
  }
  
  tree_names <- names(trees)

  # Strip out first tree to get the base topology
  base_tree <- trees[[1]]

  # Create bootstrap table
  bootstrap_tibble <- tibble(Clade=Rboretum::getTreeClades(base_tree,include_root = TRUE))
  
  # For each tree, fetch the bootstraps as identified by the clade they support (this avoids mismatching node IDs)
  for(i in 1:tree_count){
    
    tree <- trees[[i]]
    tree_name <- tree_names[[i]]
    
    subtree <- ape::subtrees(base_tree)
    subtree_length <- length(subtree)
    
    # If tree lacks node labels, add "-"
    if(is.null(tree$node.label)){
      tree$node.label <- rep("-",subtree_length)
    }
    
    clade_subtree <- subtree[2:subtree_length]
    
    clades <- purrr::map(.x=clade_subtree,.f=function(x){Rboretum::semiSorter(x$tip.label)}) %>% unlist()
    node_labels <- purrr::map(.x=clade_subtree,.f=function(x){x$node.label[[1]]}) %>% unlist()
    tree_labels <- tibble(Clade=clades,Label=node_labels)
    
    bootstrap_tibble <- bootstrap_tibble %>% left_join(tree_labels,by='Clade')
  }
  
  names(bootstrap_tibble) <- c('Clade',tree_names)
  
  # Create a named list from the tibble
  bootstrap_list <- list()
  
  for(i in 1:nrow(bootstrap_tibble)){
    clade <- bootstrap_tibble$Clade[[i]]
    node_labels <- select(bootstrap_tibble,2:(tree_count+1)) %>% slice(i) %>% unlist()
    
    # Replace "" with "-"
    node_labels[node_labels==""]  <- "-"
    
    # Replace NA with "-"
    node_labels[is.na(node_labels)]  <- "-"
    
    bootstrap_list[[clade]] <- ifelse(all(node_labels=="-"),"-",paste(node_labels,collapse = "/"))
  }
  
  # Add new labels to phylo for return
  base_tree_clades <- purrr::map(.x=subtrees(base_tree),.f=function(x){Rboretum::semiSorter(x$tip.label)}) %>% unlist()
  base_tree$node.label <- c('Root',purrr::map(.x=base_tree_clades,.f=function(x){bootstrap_list[[x]]}) %>% unlist())
  
  return(base_tree)
}