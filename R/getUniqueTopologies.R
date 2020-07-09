#' Rboretum Unique Topology Fetcher
#' 
#' This function takes a multiPhylo where all trees share 3 or more taxa, and returns the unique topologies after pruning to a common set of taxa
#' @param trees Named, rooted multiPhylo object where all trees share at least three taxa
#' @param tree_names: OPTIONAL: If TRUE, name unique trees based off of individual tree names [Default: FALSE, return trees named Topology_1, Topology_2, etc.]
#' @param print_table OPTIONAL: If TRUE, return trees, and print summary table
#' @param return_table OPTIONAL: If TRUE, return summary table rather than multiPhylo
#' @return multiPhylo containing unique topologies
#' @export

getUniqueTopologies <- function(trees,tree_names,print_table,return_table){
  
  if(!Rboretum::isMultiPhylo(trees,check_rooted = TRUE,check_three_taxa = TRUE)){
    stop("'trees' must be a rooted multiPhylo object where all trees share at least three taxa.")
  }
  
  # Get tree count
  raw_tree_count <- length(trees)
  
  # Check return/print status
  if(missing(print_table)){
    print_table <- FALSE
  } else if(!is.logical(print_table)){
    print_table <- FALSE
  } else if(length(print_table)!=1){
    print_table <- FALSE
  }
  
  if(missing(return_table)){
    return_table <- FALSE
  } else if(!is.logical(return_table)){
    return_table <- FALSE
  } else if(length(return_table)!=1){
    return_table <- FALSE
  }
  
  if(missing(tree_names)){
    tree_names <- FALSE
  } else if(!is.logical(tree_names)){
    tree_names <- FALSE
  } else if(length(tree_names)!=1){
    tree_names <- FALSE
  }
  
  # Fetch/Add names as needed
  if(!Rboretum::isMultiPhylo(trees,check_named = TRUE)){
    trees <- Rboretum::treeNamer(trees)
    raw_tree_names <- names(trees)
  } else{
    raw_tree_names <- names(trees)
  }
  
  # Reduce multiPhylo to shared taxa
  if(!Rboretum::isMultiPhylo(trees,check_all_taxa = TRUE)){
    trees <- Rboretum::treeTrimmer(trees)
  }
  
  # Get unique topologies
  unique_trees <- ape::unique.multiPhylo(trees)
  unique_count <- length(unique_trees)
  
  # If a single unique topology, combine bootstraps before return
  if(unique_count == 1){
    return_tree <- Rboretum::combineBootstraps(trees)
    
    if(print_table | return_table){
      print("Single topology detected. No table/summary to print...")
    }
    
    return(return_tree)
  } else{
    
    # If multiple topologies exist, tally and return
    unique_tree_names <- paste0("Topology_", 1:unique_count)
    names(unique_trees) <- unique_tree_names
    
    # Create topology list
    unique_tree_list <- vector("list", unique_count)
    names(unique_tree_list) <- unique_tree_names
    
    # Get tree matches
    for(tree_num in 1:raw_tree_count){
      
      temp_tree <- trees[[tree_num]]
      temp_name <- raw_tree_names[[tree_num]]
      
      for(unique_num in 1:unique_count){
        temp_unique <- unique_trees[[unique_num]]
        temp_unique_name <- unique_tree_names[[unique_num]]
        
        if(ape::all.equal.phylo(temp_tree,temp_unique,use.edge.length = FALSE)){
          unique_tree_list[[temp_unique_name]] <- c(unique_tree_list[[temp_unique_name]],temp_name)
        }
      }
    }
    
    # Get tree counts for each topolgogy
    topology_tallies <- purrr::map(.x=unique_tree_names,.f=function(x){length(unique_tree_list[[x]])}) %>% unlist()
    topology_groups <- purrr::map(.x=unique_tree_names,.f=function(x){Rboretum::semiSorter(unique_tree_list[[x]])}) %>% unlist()
    
    summary_df <- data.frame(Tree_Name = names(unique_trees),
                             Trees_with_Topology = topology_groups,
                             Tree_Count = as.integer(topology_tallies),
                             Tree_Percent = round((topology_tallies/as.numeric(raw_tree_count)*100),1),stringsAsFactors = FALSE)
    
    # Remove bootstrap values from pooled trees
    return_tree <- unique_trees
    names(return_tree) <- unique_tree_names
    
    pooled_tree_names <- summary_df %>% filter(Tree_Count > 1) %>% pull(Tree_Name)
    pooled_count <- length(pooled_tree_names)
    
    # Combine bootstrap value for pooled trees
    if(pooled_count>0){
      for(i in 1:pooled_count){
        pooled_tree_id <- pooled_tree_names[[i]]
        pooled_trees <- summary_df %>% filter(Tree_Name==pooled_tree_id) %>% pull(Trees_with_Topology) %>% semiVector()
        return_tree[[pooled_tree_id]] <- Rboretum::combineBootstraps(trees[pooled_trees])
      }
    }
    
    # Set names if requested
    if(tree_names){
      names(return_tree) <- pull(summary_df,Trees_with_Topology)
    }
    
    if(print_table){
      print(summary_df)
    }
    
    if(return_table){
      return(summary_df)
    } else{
      return(return_tree)
    }
  }
}