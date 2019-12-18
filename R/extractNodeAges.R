#' Rboretum Node Age Extractor
#'
#' This function takes a tree where branches are scaled to time and returns the estimated node ages, or a summary of node ages acrosss trees
#' @param tree Tree(s) to extract date information from. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees share 3+ taxa
#' }
#' @param return_summary OPTIONAL: If TRUE and a multiPhylo is provided, return a summary of node ages across trees with mean and median values per clade
#' @export

extractNodeAges <- function(tree,return_summary){
  
  # Ensure tree is valid
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_three_taxa = TRUE)){
    stop("extractNodeAges requires a rooted phylo or multiPhylo object where all trees share 3+ taxa")  
  }
  
  # Return summary data?
  if(missing(return_summary)){
    return_summary <- FALSE
  } else if(!is.logical(return_summary)){
    return_summary  <- FALSE
  }
  
  if(Rboretum::isPhylo(tree)){
    tree <- c(tree,tree)
    tree_count <- 1
    
  } else{
    
    tree_count <- length(tree)
    
    # Trim multiPhylo to common taxa if necessary
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree_taxa <- Rboretum::getSharedTaxa(tree)
      tree <- Rboretum::treeTrimmer(tree,tree_taxa)
    }
    
    # Add dummy names if necessary
    if(!Rboretum::isMultiPhylo(tree,check_named = TRUE)){
      tree <- Rboretum::treeNamer(tree)
    }
  }
  
  tree_df_list <- list()
  
  for(i in 1:tree_count){
    
    tree_clades <- Rboretum::getTreeClades(tree[[i]],include_root = TRUE)
    
    tree_nodes <- purrr::map(.x=tree_clades,.f=function(x){ape::getMRCA(tree[[i]],semiVector(x))}) %>% unlist() %>% as.character()
    node_ages <- purrr::map(.x=tree_nodes,.f=function(x){ape::branching.times(tree[[i]])[x] %>% as.numeric()}) %>% unlist()
    
    tree_date_df <- data.frame(Clade=as.character(tree_clades),Node=as.integer(tree_nodes),Node_Age=as.numeric(node_ages),stringsAsFactors = FALSE)
    
    if(tree_count == 1){
      return(tree_date_df)
    } else{
      tree_df_list[[i]] <- tree_date_df %>% mutate(Tree_Name=names(tree)[i])
    }
  }
  
  tree_date_df <- tree_df_list %>% bind_rows()
  
  if(return_summary){
    
    # Warn if topologies differ
    if(!Rboretum::isMultiPhylo(check_all_equal=TRUE)){
      warning("WARNING: Trees have different topologies, which may impact node age summaries")
    }
    
    tree_date_df <- tree_date_df %>% 
      select(-Node,-Tree_Name) %>% 
      group_by(Clade) %>% 
      summarise(Mean_Node_Age=mean(Node_Age),Median_Node_Age=median(Node_Age),StdDev_Node_Age=sd(Node_Age),MAD_Node_Age=mad(Node_Age))
  }
  
  return(tree_date_df)
}