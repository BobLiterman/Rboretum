#' Rboretum Tree Splitter
#'
#' This function breaks down a rooted phylo or multiPhylo object into its respective set of splits
#' @param tree Tree(s) to split. Options include:
#' \itemize{
#'   \item A rooted phylo object
#'   \item A rooted, named multiPhylo object where all trees share 3+ taxa
#' }
#' @return Dataframe (or a list of dataframes for each unique topology): Column 1: Semi-colon separated monophyletic clade; Column 2: 'Mirror Clade'; Column 3: Phylo Node ID (NA for root split)
#' @export

getTreeSplits <- function(tree){
  
  if(Rboretum::isPhylo(tree,check_rooted = TRUE)){ # If a phylo object is provided, get splits...
    split_df <- getTreeSplits_Worker(tree)
    return(split_df)
  } else if(Rboretum::isMultiPhylo(tree,check_named = TRUE,check_rooted = TRUE,check_three_taxa = TRUE,check_all_equal = TRUE)){ # If a valid multiPhylo is provided, but all trees have the same topology, return splits for the first tree
    
    common_taxa <- Rboretum::getSharedTaxa(tree)
    first_tree <- treeTrimmer(tree[[1]],common_taxa)
    split_df <- getTreeSplits_Worker(first_tree)
    return(split_df)
    
  } else if(Rboretum::isMultiPhylo(tree,check_named = TRUE,check_rooted = TRUE,check_three_taxa = TRUE,check_some_equal = TRUE)){ # If a valid multiPhylo is provided, but some trees have the same topology, return splits for the unique trees
    
    tree <- getUniqueTopologies(tree,print_table = TRUE)
    tree_names <- names(tree)
    split_df <- purrr::map(.x = tree,.f = function(x){getTreeSplits_Worker(x)})
    names(split_df) <- tree_names
    return(split_df)
    
  } else if(Rboretum::isMultiPhylo(tree,check_named = TRUE,check_rooted = TRUE,check_three_taxa = TRUE,check_all_unique = TRUE)){ # If a valid multiPhylo is provided, and all trees have a unique topology, return splits for all trees
    
    tree_names <- names(tree)
    split_df <- purrr::map(.x = tree,.f = function(x){getTreeSplits_Worker(x)})
    names(split_df) <- tree_names
    return(split_df)
    
  } else{
    stop("'tree' does not appear to be a rooted phylo, or a named, rooted multiPhylo object where all trees share 3+ taxa.")
  }
}