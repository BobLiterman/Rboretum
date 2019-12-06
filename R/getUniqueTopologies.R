#' Rboretum Unique Topology Fetcher
#' 
#' This function takes a multiPhylo where all trees share 3 or more taxa, and returns the unique topologies after pruning to a common set of taxa
#' @param trees Named, rooted multiPhylo object where all trees share at least three taxa
#' @param print_table OPTIONAL: If TRUE, return trees, and print summary table
#' @param return_table OPTIONAL: If TRUE, return summary table rather than multiPhylo
#' @return multiPhylo containing unique topologies
#' @export

getUniqueTopologies <- function(trees,print_table,return_table){
  
  if(!Rboretum::isMultiPhylo(trees,check_rooted = TRUE,check_named = TRUE,check_three_taxa = TRUE)){
    stop("'trees' must be a named, rooted multiPhylo object where all trees share at least three taxa.")
  } else if(Rboretum::checkSameTopology(trees)){ # One unique tree, return first
    return(trees[[1]])
  } else if(!Rboretum::checkSameTopology(trees,check_any = TRUE)){ # Trees are already unique
    return(trees)
  }

  if(missing(print_table)){
    print_table <- FALSE
  } else if(!is.logical(print_table)){
    print_table <- FALSE
  }
  
  if(missing(return_table)){
    return_table <- FALSE
  } else if(!is.logical(return_table)){
    return_table <- FALSE
  }
  
  tree_taxa <- Rboretum::getSharedTaxa(trees)
  
  if(!Rboretum::isMultiPhylo(trees,check_all_taxa = TRUE)){ # Trim to common taxa 
    trees <- Rboretum::treeTrimmer(trees,tree_taxa)
  }
  
  raw_tree_count <- length(trees)
  raw_tree_names <- names(trees)
  
  # Compare all tree topologies
  tree_a <- c()
  tree_b <- c()
  top_check <- c()

  for(i in 1:(raw_tree_count-1)){
    for(j in (i+1):raw_tree_count){
      tree_a <- c(tree_a,raw_tree_names[[i]])
      tree_b <- c(tree_b,raw_tree_names[[j]])
      top_check <- c(top_check,ape::all.equal.phylo(trees[[i]],trees[[j]],use.edge.length = FALSE))
    }
  }
  
  tree_compare <- data.frame(Tree_1=tree_a,Tree_2=tree_b,Same_Topology=top_check) %>% filter(Same_Topology)

  tree_groups <- list()
  grouped_trees <- c()
  unique_trees <- c()
  topology_names  <- c()
  
  for(i in 1:raw_tree_count){
    
    next_pos <- length(tree_groups) + 1

    focal_tree <- raw_tree_names[[i]]
    
    if(!focal_tree %in% grouped_trees){
      
      unique_trees <- c(unique_trees,trees[[i]])
      topology_names <- c(topology_names,paste(c('Topology_',next_pos),collapse =''))      
      
      tree_group <- tree_compare %>% filter(Tree_1 == focal_tree | Tree_2 == focal_tree)
      
      if(nrow(tree_group)==0){
        tree_groups[[next_pos]] <- focal_tree
        grouped_trees <- c(grouped_trees,focal_tree) %>% unique() %>% sort()
      }
      else{
        tree_groups[[next_pos]] <- as.vector(as.matrix(tree_compare[,c("Tree_1", "Tree_2")])) %>% unique() %>% sort()
        grouped_trees <- c(grouped_trees,as.vector(as.matrix(tree_compare[,c("Tree_1", "Tree_2")])) %>% unique() %>% sort()) %>% unique() %>% sort()
      }
    }
  }
  
  unique_count <- length(tree_groups)
  names(unique_trees) <- topology_names
  
  top_count <- c()
  top_trees <- c()
  
  for(i in 1:unique_count){
    top_count <- c(top_count,length(tree_groups[[i]]))
    top_trees <- c(top_trees,paste(tree_groups[[i]],collapse = ";"))
  }
  
  summary_df <- data.frame(Topology_ID = as.character(topology_names),
                           Trees_with_Topology = as.character(top_trees),
                           Tree_Count = as.integer(top_count),
                           Tree_Percent = round((top_count/as.numeric(raw_tree_count)*100),1))
  
  if(print_table){
    print(summary_df)
  }
  
  if(return_table){
    return(summary_df)
  } else{
    return(unique_trees)
  }
}