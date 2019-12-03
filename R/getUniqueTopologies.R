#' Rboretum Unique Topology Fetcher
#'
#' @param trees Named, rooted multiPhylo object
#' @param return_table OPTIONAL: If TRUE, return summary table rather than multiPhylo
#' @return multiPhylo containing unique topologies; or, summary table of the same data
#' @export

getUniqueTopologies <- function(trees,return_table){
  
  if(!Rboretum::isMultiPhylo(trees,check_rooted = TRUE,check_names = TRUE,check_three_taxa = TRUE)){
    stop("'trees' must be a named, rooted multiPhylo object where all trees share at least three taxa.")
  } else if(Rboretum::checkSameTopology(trees)){
    stop("All 'trees' have the same topology.")
  } else  if(Rboretum::isMultiPhylo(trees,check_unique = TRUE)){
    print("All 'trees' have a unique topology. Returning raw trees...")
    return(trees)
  } else if(!Rboretum::checkSameTaxa(trees)){
    print('Trimming trees to identical taxa sets...')
    trees <- Rboretum::treeTrimmer(trees,tree_taxa)
  }
  
  if(missing(return_table)){
    return_table <- FALSE
  } else if(!is.logical(return_table)){
    return_table <- FALSE
  }
  
  tree_compare <- compareTrees(trees) %>%
    filter(!Comparable) %>%
    select(Tree_1,Tree_2)
  
  tree_groups <- list()
  grouped_trees <- c()
  
  for(i in 1:nrow(tree_compare)){
    
    tree_1 <- as.character(tree_compare$Tree_1[i])
    tree_2 <- as.character(tree_compare$Tree_2[i])
    
    list_length <- length(tree_groups)
    topology_name <- paste(c('Topology_',list_length+1),collapse ='')
    
    if(list_length == 0){
      tree_groups[[1]] <- c(tree_1,tree_2)
      names(tree_groups) <- topology_name
      grouped_trees <- sort(c(tree_1,tree_2))
    } else{
      t1_check <- tree_1 %in% grouped_trees
      t2_check <- tree_2 %in% grouped_trees
      
      if(!(t1_check & t2_check)){
        if(t1_check){
          group <- names(tree_groups)[sapply(seq_along(tree_groups),function(x){tree_1 %in% tree_groups[[x]]})]
          tree_groups[[group]] <- sort(c(tree_groups[[group]],tree_2))
          grouped_trees <- sort(c(grouped_trees,tree_2))
        } else if(t2_check){
          group <- names(tree_groups)[sapply(seq_along(tree_groups),function(x){tree_2 %in% tree_groups[[x]]})]
          tree_groups[[group]] <- sort(c(tree_groups[[group]],tree_1))
          grouped_trees <- sort(c(grouped_trees,tree_1))
        } else{
          tree_groups[[topology_name]] <- c(tree_1,tree_2)
          grouped_trees <- sort(c(grouped_trees,tree_1,tree_2))
        }
      }
    }
  }
  
  tree_names <- names(trees)
  solo_trees <- tree_names[!tree_names %in% grouped_trees]
  
  if(length(solo_trees)>=1){
    for(i in 1:length(solo_trees)){
      list_length <- length(tree_groups)
      topology_name <- paste(c('Topology_',list_length+1),collapse ='')
      tree_groups[[topology_name]] <- solo_trees[i]
    }
  }
  
  unique_count <- length(tree_groups)
  
  unique_trees <- c()
  
  for(i in 1:unique_count){
    if(i == 1){
      unique_trees <- trees[[tree_groups[[1]][1]]]
    } else{
      unique_trees <- c(unique_trees,trees[[tree_groups[[i]][1]]])
    }
    tree_groups[[i]] <- c(length(tree_groups[[i]]),paste(sort(tree_groups[[i]]),collapse = ';'))
  }
  
  top_names <- names(tree_groups)
  top_count <- lapply(tree_groups,function(x){as.integer(x[[1]])}) %>% unlist() %>% as.integer()
  top_trees <- lapply(tree_groups,function(x){as.character(x[[2]])}) %>% unlist() %>% as.character()

  if(return_table){
    summary_df <- data.frame(Topology_ID = as.character(top_names),Trees_with_Topology = as.character(top_trees), Tree_Count = as.integer(top_count),Tree_Percent = round((top_count/as.numeric(length(trees))*100),1))
    return(summary_df)
  } else{
    names(unique_trees) <- top_names
    return(unique_trees)
  }
}