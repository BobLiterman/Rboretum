#' Rboretum Tree Splitter
#'
#' This function breaks down a rooted phylo or multiPhylo object into its respective set of splits
#' @param tree Tree(s) to split. Options include:
#' \itemize{
#'   \item A rooted phylo object
#'   \item A rooted, named multiPhylo object where all trees share 3+ taxa
#' }
#' @return Dataframe (or a list of dataframes for each topology): Column 1: Semi-colon separated monophyletic clade; Column 2: 'Mirror Clade'; Column 3: Phylo Node ID (NA for root split)
#' @export

getTreeSplits <- function(tree){
  
  if(Rboretum::isPhylo(tree,check_rooted = TRUE)){ # If a phylo object is provided, get splits...
    split_df <- getTreeSplits_Worker(tree)
    return(split_df)
  } else if(Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_three_taxa = TRUE)){ # If a multiPhylo, trim to common species and get splits...
    
    # Name trees if unnamed
    if(!Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree <- Rboretum::treeNamer(tree)
    }
    
    # Trim to commmon taxa if necessary
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree <- treeTrimmer(tree)
    }
    
    # If all topologies in multiPhylo match, return a single dataframe containing the splits
    if(Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
      tree <- tree[[1]]
      split_df <- getTreeSplits_Worker(tree)
      return(split_df)
    } else{ # If multiple topologies exist in the multiPhylo, return splits for each tree
      tree_names <- names(tree)
      split_df <- purrr::map(.x = tree,.f = function(x){getTreeSplits_Worker(x)}) %>% `names<-`(tree_names)
      return(split_df)
    }
  } else{
    stop("'tree' does not appear to be a rooted phylo, or a named, rooted multiPhylo object where all trees share 3+ taxa.")
  }
}