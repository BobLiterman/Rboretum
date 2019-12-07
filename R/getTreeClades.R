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
#' @return Dataframe (or a list of dataframes for each unique topology): Column 1: Semi-colon separated monophyletic clade; Column 2: 'Mirror Clade'; Column 3: Phylo Node ID (NA for root split)
#' @export

getTreeClades <- function(tree,include_root,print_counts,return_counts){
  
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
    
    # Get tree clades
    listed_tree_clades <- purrr::map(.x = tree_splits ,.f = function(x){  pull(x,Clade) %>% as.character()})
    tree_clade_vec <- as.character(unlist(tree_clades))

    if(include_root){
      
      # Extract root clades if requested
      tree_names <- names(listed_tree_clades)
      
      listed_root_clades <- purrr::map(.x = tree_splits ,.f = function(x){ x %>% filter(is.na(Split_Node)) %>% select(Clade,Mirror_Clade) %>% slice(1) %>%  unlist() %>% as.character()})
      root_clade_vec <- as.character(unlist(root_clades))

      listed_tree_clades <- purrr::map(.x=tree_names,.f=function(x){c(listed_tree_clades[[x]],listed_root_clades[[x]])})
      names(listed_tree_clades) <- tree_names
      tree_clade_vec <- c(tree_clade_vec,root_clade_vec)
      
    }
    
    clade_table <- table(tree_clade_vec)
    
    trees_with_clade <- purrr::map(.x=names(clade_table),.f=function(x){paste(sort(names(list.search(listed_tree_clades,expr = x %in% .))),collapse = ";")}) %>% as.character()
    
    clade_sorter <- data.frame(Clade=as.character(names(clade_table)),Count=as.integer(clade_table),Trees=as.character(trees_with_clade)) %>% 
      arrange(desc(Count),Clade)
    
    if(print_counts){
      print(clade_sorter)
    }
    
    if(return_counts){
      return(clade_sorter)
    } else{
      sorted_clades <- pull(clade_sorter,Clade) %>% as.character() %>% sort()
      return(sorted_clades)
    }
    
  } else{ stop("Unknown error.") }
}