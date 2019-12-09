#' Rboretum Clade Fetcher
#'
#' This function breaks down a rooted phylo or multiPhylo object into its respective monophyletic clades (root clades excluded unless noted), sorted alphabetically.
#' @param tree Tree(s) to split. Options include:
#' \itemize{
#'   \item A rooted phylo object
#'   \item A rooted, named multiPhylo object where all trees share 3+ taxa
#' }
#' @param include_root OPTIONAL: If TRUE, return two root clades [Default: FALSE, exlclude root clades]
#' @param print_counts OPTIONAL: If TRUE, print a summary table of unique splits and how many/which trees contain them [Default: FALSE, no printing]
#' @param return_counts OPTIONAL: If TRUE, instead of returning a vector of all clades, return a summary table of unique splits and how many/which trees contain them [Default: FALSE, return clades as sorted vector]
#' @param return_shared OPTIONAL: If TRUE, instead of returning a vector of all clades, return a vector of only those clades present in all trees [Default: FALSE, return all clades]
#' @return Either:
#' \itemize{
#'   \item A character vector of monophyletic clades; or,
#'   \item A dataframe containing information about clades and their occurrence among trees
#' }
#' @export

getTreeClades <- function(tree,include_root,print_counts,return_counts,return_shared){
  
  # Check if tree is valid
  if(!Rboretum::isMultiPhylo(tree,check_rooted = TRUE, check_named = TRUE, check_three_taxa = TRUE) & !Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'tree' does not appear to be a valid rooted phylo object, or a named, rooted multiPhylo")
  }

  # Return information about the two root splits?
  if(missing(include_root)){
    include_root <- FALSE
  } else if(!is.logical(include_root)){
    include_root <- FALSE
  }
  
  # Print summary table?
  if(missing(print_counts)){
    print_counts <- FALSE
  } else if(!is.logical(print_counts)){
    print_counts <- FALSE
  }
    
  # Return summary table instead of clade vector?
  if(missing(return_counts)){
    return_counts <- FALSE
  } else if(!is.logical(return_counts)){
    return_counts <- FALSE
  }
  
  # Return shared clades instead of all clades?
  if(missing(return_shared)){
    return_shared <- FALSE
  } else if(!is.logical(return_shared)){
    return_shared <- FALSE
  }
  
  if(return_shared & return_counts){
    stop("Choose either 'return_shared' or 'return_counts'...")
  }
  
  # Get tree splits
  tree_splits <- Rboretum::getTreeSplits(tree)
  
  if(Rboretum::isPhylo(tree)){
    
    # Get non-root clades  
    tree_clades <- tree_splits %>% filter(!is.na(Split_Node)) %>%
      pull(Clade) %>% as.character() %>% sort()
    
    if(include_root){ 
      
      # Extract root clades if requested
      root_clades <- tree_splits %>% filter(is.na(Split_Node)) %>% 
        select(Clade,Mirror_Clade) %>% slice(1) %>% 
        unlist() %>% as.character()
      
      tree_clades <- c(tree_clades,root_clades) %>% sort()
      
      return(tree_clades)
      
    } else{
      
      return(tree_clades)
    }
    
  } else if(Rboretum::isMultiPhylo(tree)){
    
    tree_count <- length(tree)
    
    # Trim to common taxa, if necessary
    shared_taxa <- Rboretum::getSharedTaxa(tree)
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree <- Rboretum::treeTrimmer(tree,shared_taxa)
    }
    
    # Get tree clades
    
    listed_tree_clades <- purrr::map(.x = tree ,.f = function(x){ Rboretum::getTreeSplits_Worker(x) %>% filter(!is.na(Split_Node)) %>% pull(Clade) %>% as.character()})
    listed_root_clades <- purrr::map(.x = tree ,.f = function(x){ Rboretum::getTreeSplits_Worker(x) %>% filter(is.na(Split_Node)) %>% select(Clade,Mirror_Clade) %>% slice(1) %>%  unlist() %>% as.character()})
    
    if(include_root){
      listed_tree_clades <- purrr::map(.x=1:tree_count,.f=function(x){c(listed_tree_clades[[x]],listed_root_clades[[x]])})
      names(listed_tree_clades) <- names(listed_root_clades)
    }
    
    clade_table <- table(unlist(listed_tree_clades))
    
    trees_with_clade <- purrr::map(.x=names(clade_table),.f=function(x){paste(sort(names(list.search(listed_tree_clades,expr = x %in% .))),collapse = ";")}) %>% as.character()
    
    clade_sorter <- data.frame(Clade=as.character(names(clade_table)),Count=as.integer(clade_table),Trees=as.character(trees_with_clade)) %>% 
      arrange(desc(Count),Clade)
    
    if(print_counts){
      print(clade_sorter)
    }
    
    if(return_shared){
      shared_clades <- clade_sorter %>%
        filter(Count == tree_count)
      
      if(nrow(shared_clades)>0){
        shared_clades <- shared_clades %>%
          pull(Clade) %>% as.character() %>% sort()
        return(shared_clades)
      } else{
        stop("Trees in 'trees' share no clades...")
      }
    } else if(return_counts){
      return(clade_sorter)
    } else{
      sorted_clades <- pull(clade_sorter,Clade) %>% as.character() %>% sort()
      return(sorted_clades)
    }
    
  } else{ stop("Unknown error.") }
}