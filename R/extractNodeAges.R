#' Rboretum Node Age Extractor
#'
#' This function takes a tree(s) where branches are scaled to time and returns the estimated node ages, or a summary of node ages across trees
#' @param tree Tree(s) to extract date information from. Must be ultrametric. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees share 3+ taxa and support a single topology 
#' }
#' @param return_summary OPTIONAL: If a multiPhylo is provided, function will retrurn the mean, median, or both set of node ages across datasets.
#' \itemize{
#'   \item mean: Return mean node age for each node
#'   \item median: Return median node age for each node
#'   \item both: Return mean, median, and summary statstics for each node
#' }
#' @return A dataframe containing (a) node/date information for each tree or (b) clade-level ages summarized across trees in a multiPhlyo if return_summary is set to 'mean',',median', or 'both'
#' @export

extractNodeAges <- function(tree,return_summary){
  
  # Ensure tree is valid
  if(missing(tree)){
    stop("extractNodeAges requires a ultrametric, rooted phylo object, or a rooted set of ultrametric multiPhylo trees that all support a common topology.")
  } else if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_three_taxa = TRUE,check_all_equal = TRUE)){
    stop("extractNodeAges requires a rooted phylo object or a multiPhylo object where all trees share 3+ taxa and a common topology...")  
  }
  
  if(Rboretum::isPhylo(tree)){
    if(!ape::is.ultrametric.phylo(tree)){
      stop("extractNodeAges requires an ultrametric phylo object...")
    }
  }
  
  # If multiPhylo, ensure a single topology and ultrametric
  if(Rboretum::isMultiPhylo(tree)){
    if(Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_three_taxa = TRUE) & !Rboretum::isMultiPhylo(tree)){
      stop("Topologies vary among trees in supplied multiPhylo. getNodeAges only accepts a single tree, or a multiPhylo where all trees share a single topology...")
    } else if(!purrr::map(.x=tree,.f=function(x){ape::is.ultrametric.phylo(x)}) %>% unlist() %>% all()){
      stop("extractNodeAges requires an ultrametric multiPhylo object...")
    }
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
    
    # Trim multiPhylo to common taxa if necessary
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree <- Rboretum::treeTrimmer(tree)
    }
    
    # Add dummy names if necessary
    if(!Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree <- Rboretum::treeNamer(tree)
    }
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
    if(summary_col == 'both'){
      tree_date_df <- tree_date_df %>% 
        select(-Tree_Name) %>% 
        group_by(Clade) %>% 
        summarise(Mean_Node_Age=mean(Node_Age),Median_Node_Age=median(Node_Age),StdDev_Node_Age=sd(Node_Age),MAD_Node_Age=mad(Node_Age))
    } else if(summary_col == 'mean'){
      tree_date_df <- tree_date_df %>% 
        select(-Tree_Name) %>% 
        group_by(Clade) %>% 
        summarise(Mean_Node_Age=mean(Node_Age)) %>%
        select(Clade,Mean_Node_Age) %>% `names<-`(c('Clade','Node_Age'))
    } else if(summary_col == 'median'){
      tree_date_df <- tree_date_df %>% 
        select(-Tree_Name) %>% 
        group_by(Clade) %>% 
        summarise(Median_Node_Age=median(Node_Age)) %>%
        select(Clade,Median_Node_Age) %>% `names<-`(c('Clade','Node_Age'))
    }
  }
  return(tree_date_df)
}