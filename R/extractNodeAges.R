#' Rboretum Node Age Extractor
#'
#' This function takes a tree(s) where branches are scaled to time and returns the estimated node ages, or a summary of node ages across trees
#' @param tree Tree(s) to extract date information from. Must be ultrametric. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees share 3+ taxa and support a single topology 
#' }
#' @param return_summary OPTIONAL: If a multiPhylo is provided, function will return the mean, median, or both set of node ages across datasets.
#' \itemize{
#'   \item mean: Return mean node age for each node
#'   \item median: Return median node age for each node
#'   \item both: Return mean, median, and summary statistics for each node
#' }
#' @return A dataframe containing (a) node/date information for each tree or (b) clade-level ages summarized across trees in a multiPhlyo if return_summary is set to 'mean',',median', or 'both'
#' @export

extractNodeAges <- function(tree,return_summary){
  
  # Ensure tree is valid
  if(missing(tree)){
    stop("extractNodeAges requires a ultrametric, rooted phylo object, or a rooted set of ultrametric multiPhylo trees that all support a common topology.")
  }
  
  # If phylo, ensure ultrametric
  if(Rboretum::isPhylo(tree)){
    if(!ape::is.ultrametric.phylo(tree)){
      stop("extractNodeAges requires an ultrametric phylo object...")
    }
  } else if(Rboretum::isMultiPhylo(tree,check_three_taxa = TRUE,check_rooted = TRUE)){ # If multiPhylo, ensure ultrametric and single topology
    
    # Add dummy names if necessary
    if(!Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree <- Rboretum::treeNamer(tree)
    }
    
    # Trim multiPhylo to common taxa if necessary
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree <- Rboretum::treeTrimmer(tree)
    }
    
    # Ensure a single topology after trimming
    if(!Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
      stop("extractNodeAges requires a rooted phylo object or a multiPhylo object where all trees share 3+ taxa and a common topology...")  
    }
    
    # Ensure all trees are ultrametric
    if(!purrr::map(.x=tree,.f=function(x){ape::is.ultrametric.phylo(x)}) %>% unlist() %>% all()){
      stop("extractNodeAges requires an ultrametric multiPhylo object...")
    }
  } else{
    stop("extractNodeAges requires a rooted phylo object or a multiPhylo object where all trees share 3+ taxa")
  }
  
  # Return summary data?
  if(missing(return_summary)){
    return_summary <- FALSE
  } else if(!is.character(return_summary)){
    stop("'return_summary' should be either 'mean','median',or 'both'")
  } else if(!return_summary %in% c('mean','median','both')){
    stop("'return_summary' should be either 'mean','median',or 'both'")
  } else if(return_summary == 'mean'){
    summary_col <- 'mean'
    return_summary <- TRUE
  } else if(return_summary == 'median'){
    summary_col <- 'median'
    return_summary <- TRUE
  } else{
    summary_col <- 'both'
    return_summary <- TRUE
  }
  
  # Can't summarize a single tree
  if(Rboretum::isPhylo(tree)){
    return_summary  <- FALSE
    tree <- c(tree,tree) # Dummy tree for processing
    tree_count <- 1
  } else{
    
    tree_count <- length(tree)
  }
  
  tree_df_list <- list()
  
  for(i in 1:tree_count){
    
    # If trees have node labels, node IDs can't be used to pull branching times
    no_bs_tree <- tree[[i]]
    no_bs_tree$node.label <- NULL
    no_bs_subtree <- subtrees(no_bs_tree)
    no_bs_branching_times <- branching.times(no_bs_tree)

    tree_clades <- purrr::map(.x=no_bs_subtree,.f=function(x){semiSorter(x$tip.label)}) %>% unlist()
    tree_nodes <- purrr::map(.x=tree_clades,.f=function(x){ape::getMRCA(no_bs_tree,semiVector(x))}) %>% unlist()
    node_ages <- purrr::map(.x=tree_nodes,.f=function(x){no_bs_branching_times[[as.character(x)]] %>% as.numeric()}) %>% unlist()
    
    tree_date_df <- tibble(Clade=as.character(tree_clades),Node_Age=as.numeric(node_ages))
    
    if(tree_count == 1){
      return(tree_date_df)
    } else{
      tree_df_list[[i]] <- tree_date_df %>% mutate(Tree_Name=names(tree)[i])
    }
  }
  
  tree_date_df <- tree_df_list %>% bind_rows()
  
  if(return_summary){
    
    tree_date_df <- tree_date_df %>% 
      select(-Tree_Name) %>% 
      group_by(Clade) %>% 
      summarise(Mean_Node_Age=mean(Node_Age),Median_Node_Age=median(Node_Age),StdDev_Node_Age=sd(Node_Age),MAD_Node_Age=mad(Node_Age)) %>%
      rowwise() %>%
      mutate(CI_95_Low = Mean_Node_Age - ((qnorm(0.975)*StdDev_Node_Age)/sqrt(tree_count)),
             CI_95_High = Mean_Node_Age + ((qnorm(0.975)*StdDev_Node_Age)/sqrt(tree_count))) %>%
      ungroup()
  
    if(summary_col == 'mean'){
      tree_date_df <- tree_date_df %>% 
        select(Clade,Mean_Node_Age,StdDev_Node_Age,CI_95_Low,CI_95_High) %>% `names<-`(c('Clade','Mean_Node_Age','StdDev_Node_Age','CI_95_Low','CI_95_High'))
    } else if(summary_col == 'median'){
      tree_date_df <- tree_date_df %>% 
        select(Clade,Median_Node_Age,MAD_Node_Age) %>% `names<-`(c('Clade','Median_Node_Age','MAD_Node_Age'))
    }
  }
  return(tree_date_df)
}