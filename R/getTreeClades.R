#' Rboretum Clade Fetcher
#'
#' This function breaks down a rooted phylo or multiPhylo object into its respective monophyletic groups, sorted alphanumerically.
#' @param tree Tree(s) to split. Options include:
#' \itemize{
#'   \item A rooted phylo object
#'   \item A rooted multiPhylo object where all trees share 3+ taxa
#' }
#' @param include_root OPTIONAL: If FALSE, don't include root clades as part of the return (Automatic if tree is rooted on a single taxon) [Default: TRUE, include root clades]
#' @param print_counts OPTIONAL: If TRUE, print a summary table of unique splits and how many/which trees contain them [Default: FALSE, no printing]
#' @param return_counts OPTIONAL: If TRUE, instead of returning a vector of clades, return a summary table of unique splits and how many/which trees contain them [Default: FALSE, return clades as sorted vector]
#' @param return_shared OPTIONAL: If TRUE, instead of returning a vector of clades, return a vector of only those clades present in all trees [Default: FALSE, return all clades]
#' @return Either:
#' \itemize{
#'   \item A character vector of monophyletic groups; or,
#'   \item A dataframe containing information about clades and their occurrence among trees
#' }
#' @export

getTreeClades <- function(tree,include_root,print_counts,return_counts,return_shared){
  
  # Check if tree is valid
  if(missing(tree)){
    stop("getTreeClades requires a valid rooted phylo object, or a rooted multiPhylo where all trees share 3+ taxa.")
  } else if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE, check_three_taxa = TRUE)){
    stop("getTreeClades requires a valid rooted phylo object, or a rooted multiPhylo where all trees share 3+ taxa.")
  }

  # Return information about the two root splits?
  if(missing(include_root)){
    include_root <- TRUE
  } else if(!is.logical(include_root)){
    include_root <- TRUE
  } else if(length(include_root)!=1){
    stop("include_root should be TRUE or FALSE")
  }
  
  # Print summary table?
  if(missing(print_counts)){
    print_counts <- FALSE
  } else if(!is.logical(print_counts)){
    print_counts <- FALSE
  } else if(length(include_root)!=1){
    stop("print_counts should be TRUE or FALSE")
  }
    
  # Return summary table instead of clade vector?
  if(missing(return_counts)){
    return_counts <- FALSE
  } else if(!is.logical(return_counts)){
    return_counts <- FALSE
  } else if(length(return_counts)!=1){
    stop("return_counts should be TRUE or FALSE")
  }
  
  # Return shared clades instead of all clades?
  if(missing(return_shared)){
    return_shared <- FALSE
  } else if(!is.logical(return_shared)){
    return_shared <- FALSE
  } else if(length(return_shared)!=1){
    stop("return_shared should be TRUE or FALSE")
  }
  
  if(return_shared & return_counts){
    stop("Choose either 'return_shared' or 'return_counts'...")
  }
  
  if(Rboretum::isPhylo(tree)){
    
    # Get tree splits
    tree_splits <- Rboretum::getTreeSplits(tree)
    
    # Get non-root clades  
    tree_clades <- tree_splits %>% filter(!is.na(Split_Node)) %>%
      pull(Clade) %>% as.character()

    if(include_root & !Rboretum::detectSingleRoot(tree)){
      
      # Extract root clades if requested, unless tree is rooted on a single taxon
      root_clades <- tree_splits %>% filter(is.na(Split_Node)) %>% 
        select(Clade,Mirror_Clade) %>% slice(1) %>% 
        unlist() %>% as.character()
      
      tree_clades <- c(tree_clades,root_clades)
    }
    return(naturalsort(tree_clades))
    
  } else if(Rboretum::isMultiPhylo(tree)){ # If a multiPhylo is provided...

    tree_count <- length(tree)
    
    # Add names, if necessary
    if(!Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree <- treeNamer(tree)
    }
    
    # Trim to common taxa, if necessary
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree <- Rboretum::treeTrimmer(tree)
    }
    
    # If all topologies in multiPhylo match, return a single vector of clades
    if(Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
      print("All trees supplied to getTreeClades share a common topology...returning results from common topology...")
      tree <- tree[[1]]
      tree_splits <- Rboretum::getTreeSplits(tree)
      
      # Get non-root clades  
      tree_clades <- tree_splits %>% filter(!is.na(Split_Node)) %>%
        pull(Clade) %>% as.character()
      
      if(include_root & !Rboretum::detectSingleRoot(tree)){
        
        # Extract root clades if requested, unless tree is rooted on a single taxon
        root_clades <- tree_splits %>% filter(is.na(Split_Node)) %>% 
          select(Clade,Mirror_Clade) %>% slice(1) %>% 
          unlist() %>% as.character()
        
        tree_clades <- c(tree_clades,root_clades)
      }
      
      return(naturalsort(tree_clades))
      
    } else{
      tree_names <- names(tree)
      
      # Detect trees rooted on a single taxon
      root_check <- purrr::map(.x=tree,.f=function(x){Rboretum::detectSingleRoot(x,return_root = FALSE)}) %>% unlist()
      
      # Get non-root clades
      listed_tree_clades <- purrr::map(.x = tree ,.f = function(x){ Rboretum::getTreeSplits_Worker(x) %>% filter(!is.na(Split_Node)) %>% pull(Clade) %>% as.character()})
      clade_names <- names(listed_tree_clades)
      
      # Get root clades
      listed_root_clades <- purrr::map(.x = tree ,.f = function(x){ Rboretum::getTreeSplits_Worker(x) %>% filter(is.na(Split_Node)) %>% select(Clade,Mirror_Clade) %>% slice(1) %>% unlist() %>% as.character()})
      
      # Extract root clades if requested, except for trees rooted on a single taxon
      if(include_root){
          tree_pos <- 1:tree_count
          root_index <- tree_pos[!root_check]
          for(i in root_index){
            listed_tree_clades[[i]] <- naturalsort(c(listed_tree_clades[[i]],listed_root_clades[[i]]))
          }
          names(listed_tree_clades) <- clade_names
        }
        
        clade_table <- table(unlist(listed_tree_clades))
        
        trees_with_clade <- purrr::map(.x=names(clade_table),.f=function(x){paste(naturalsort(names(list.search(listed_tree_clades,expr = x %in% .))),collapse = ";")}) %>% as.character()
        
        clade_sorter <- tibble(Clade=as.character(names(clade_table)),Count=as.integer(clade_table),Trees=as.character(trees_with_clade)) %>% 
          arrange(desc(Count),Clade)
        
        if(print_counts){
          print(clade_sorter)
        }
        
        if(return_shared){
          shared_clades <- clade_sorter %>%
            filter(Count == tree_count)
          
          if(nrow(shared_clades)>0){
            shared_clades <- shared_clades %>%
              pull(Clade) %>% as.character() %>% naturalsort()
            return(shared_clades)
          } else{
            stop("Trees in 'trees' share no clades...")
          }
        } else if(return_counts){
          return(clade_sorter)
        } else{
          sorted_clades <- pull(clade_sorter,Clade) %>% as.character() %>% naturalsort()
          return(sorted_clades)
        }
    }
  }
}