#' Rboretum Tree Dater
#'
#' This function takes a tree(s) with estimated branch lengths, and infers ages for each tree node using the 'chronos' function in 'ape'. 
#' @param tree Tree(s) to extract date information from. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees share 3+ taxa and support a single topology 
#' }
#' @param calibration_df Dataframe with 4 columns: (1) Taxon 1 (2) Taxon 2 (3) Min divergence time (4) Max divergence time
#' @param iterations How many times to estimate the age of each node prior to summarizing [Default: 1000]
#' @param return_cladogram OPTIONAL: If TRUE, return cladogram [Default: FALSE, return dataframe with node age information]
#' @return 
#' @export

treeDater <- function(tree,calibration_df,iterations,return_cladogram){
  
  if(missing(tree)){
    stop("treeDater requires a rooted phylo object, or a rooted set of multiPhylo trees that all support a common topology.")
  }
  
  # Ensure tree is valid
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_three_taxa = TRUE)){
    stop("treeDater requires a rooted phylo or multiPhylo object where all trees share 3+ taxa")  
  }
  
  # If multiPhylo, ensure a single topology
  if(Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_three_taxa = TRUE) & !Rboretum::isMultiPhylo(tree,check_all_equal = TRUE)){
    stop("Topologies vary among trees in supplied multiPhylo. getNodeAges only accepts a single tree, or a multiPhylo where all trees share a single topology.")
  }
  
  # Get tree taxa
  if(Rboretum::isPhylo(tree)){
    tree_taxa <- tree$tip.label
    tree <- c(tree,tree) # Create dummy multiPhylo for simple handling
    tree_count <- 1
  } else{
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree <- treeTrimmer(tree,tree_taxa)
      tree_count <- length(tree)
      return_tree <- tree
    } else{
      tree_count <- length(tree)
      return_tree <- tree
    }
  }
  
  # Ensure at least one calibration point
  if(missing(calibration_df)){
    stop("treeDater requires min/max divergence time estimates for at least one node in the tree to convert relative times to absoute times.")
  } else if(!is.data.frame(calibration_df)){
    stop("'calibration_df' should be a data frame.")
  } else if(!ncol(calibration_df)==4){
    stop("'calibration_df' should have 4 columns. (1) Taxon 1 (2) Taxon 2 (3) Min divergence time (4) Max divergence time")
  }
  colnames(calibration_df) <- c('Taxon_1','Taxon_2','Min','Max')
  calibration_df <- as.data.frame(calibration_df) %>%
    mutate_if(is.factor, as.character)
  
  # Check that calibration taxa exist in tree
  
  cal_taxa <- select(calibration_df,starts_with('Taxon')) %>% unlist() %>% unique()
  
  if(!all(cal_taxa %in% tree_taxa)){
    stop("Calibration data in 'calibration_df' contains information about taxa not present in 'tree'")
  } else{
    calibration_df$Two_Names <- as.character(paste(c(calibration_df$Taxon_1,calibration_df$Taxon_2),collapse = ";"))
    calibration_df <- calibration_df %>% rowwise() %>% mutate(Two_Names = semiSorter(Two_Names))
  }
  # Check that all min divergence times are <= max divergence times
  if(!all(calibration_df$Min <= calibration_df$Max)){
    stop("The minimum divergence time estimates for some calibration data in 'calibration_df' are greater than their associated maximum divergence time estimate")
  }
  
  # Set iterations
  if(missing(iterations)){
    iterations <- 1000
  } else if(!is.numeric(iterations)){
    iterations <- 1000    
  } else if(iterations < 1){
    iterations <- 1000
  }
  
  # Process node dating for each tree
  for(i in 1:tree_count){
    
    date_tree <- tree[[i]]
    
    # Get nodes for MRCA for each calibration point
    node <- purrr::map(.x=calibration_df$Two_Names,.f=function(x){ape::getMRCA(date_tree,tip=semiVector(x))}) %>% unlist()
    
    # Ensure all nodes are unique
    if(any(duplicated(node))){
      print(node[duplicated(node)])
      stop("The nodes printed above have two different sets of node calibrations")
    }
    
    age.min <- calibration_df$Min
    age.max <- calibration_df$Max
    tree_cal <- data.frame(node, age.min, age.max) %>% mutate(soft.bounds=FALSE)
    
    # Create dummy rows
    tree_edge_list <- compute.brlen(date_tree,1)$edge.length
    tree_branch_list <- branching.times(compute.brlen(date_tree,1))
    
    # Iterate chronos date estimation
    for(j in 1:iterations){
      tree_chronos_iter <- chronos(date_tree,calibration = tree_cal)
      tree_edge_list <- rbind(tree_edge_list, tree_chronos_iter$edge.length)
      tree_branch_list <- rbind(tree_branch_list, branching.times(tree_chronos_iter))
    }
    
    # Remove dummy rows
    tree_edge_list <- tree_edge_list[-1,]
    tree_branch_list <- tree_branch_list[-1,]
    
    tree_branches <- colnames(tree_branch_list)
    tree_medians <- c()
    
    for(j in 1:ncol(tree_branch_list)){
      tree_medians[colnames(tree_branch_list)[i]] <- median(tree_branch_list[,j])
    }
    
    tree_chronos_export <- compute.brlen(date_tree,1)
    
    tree_median_edge <- c()
    for(j in 1:ncol(tree_edge_list)){
      tree_median_edge<-c(tree_median_edge,median(tree_edge_list[,j]))
    }
    
    tree_chronos_export$edge.length <- tree_median_edge
    if(!is.ultrametric(tree_chronos_export)){
      tree_chronos_export <- force.ultrametric(tree_chronos_export,"extend")
    }
    
    if(tree_count == 1){
      if(missing(return_cladogram)){
        return(Rboretum::extractNodeAges(tree))
      } else if(return_cladogram){
        return(tree_chronos_export)
      }
    } else{
      return_tree[[i]] <- tree_chronos_export
    }
  }
  
  if(missing(return_cladogram)){
    return(Rboretum::extractNodeAges(return_tree,return_summary=TRUE))
  } else if(return_cladogram){
    return(return_tree)
  }
}