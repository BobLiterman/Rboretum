#' Rboretum MultiPhylo Summarizer
#'
#' This function breaks down a rooted multiPhylo object into a number of summary values
#' @param trees A rooted multiPhylo object where all trees share 3+ taxa
#' @export

summarizeMultiPhylo <- function(trees){
  
  if(!Rboretum::isMultiPhylo(trees,check_rooted = TRUE,check_three_taxa = TRUE)){ stop("'trees' does not appear to be a valid, rooted multiPhylo object where all trees share 3+ taxa.") }
  
  # Ensure that all trees are not the same topology
  if(Rboretum::isMultiPhylo(trees,check_all_equal = TRUE)){ stop("All trees in 'trees' have the same topologly after trimming to the common taxa. Nothing to compare.") }
  
  # Add tree names if unnamed
  if(!Rboretum::isMultiPhylo(trees,check_named = TRUE)){
    print("Trees unnamed, adding dummy names...")
    print('Command: treeNamer(trees)')
    trees <- treeNamer(trees)
  }
  
  tree_count <- length(trees)
  tree_names <- names(trees)
  
  # Trim trees to common species if necessary
  common_taxa <- Rboretum::getSharedTaxa(trees)
  
  print(paste(c('Read in',tree_count,'trees, that all share',length(common_taxa),'tip labels...'),collapse = ' '))
  
  if(!Rboretum::isMultiPhylo(trees,check_all_taxa = TRUE)){
    print('Trees do not have the same taxa, trimming to common tip labels...')
    print('Command: treeTrimmer(trees,getSharedTaxa(trees)')
    trees <- Rboretum::treeTrimmer(trees,common_taxa)
  } else{
    print('All trees contain identical tip labels...')
  }
  
  all_clades <- getTreeClades(trees,return_counts = TRUE)
  
  print(paste(c('Among all trees, and discounting root splits, there are',nrow(all_clades),'unique monophyletic clades...'),collapse = ' '))
  print('Command: getTreeClades(trees,return_counts = TRUE)')
  
  # Check if all trees are unique
  if(!Rboretum::isMultiPhylo(trees,check_all_unique = TRUE)){
    raw_tree_count <- tree_count
    tree_table <- Rboretum::getUniqueTopologies(trees,return_table = TRUE)
    trees <- Rboretum::getUniqueTopologies(trees)
    tree_count <- length(trees)
    tree_names <- names(trees)
    print(paste(c('Of',raw_tree_count,'raw trees, there were',tree_count,'unique topologies...'),collapse = ' '))
    print('Command: getUniqueTopologies(trees,print_table = TRUE)')
    print(tree_table)
  } else{
    print('All trees have a unique topology...')
  }
}