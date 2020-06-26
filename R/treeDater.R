#' Rboretum Tree Dater
#'
#' This function takes a tree(s) with estimated branch lengths, and infers ages for each tree node using the 'chronos' function in 'ape', or through an implementation of the 'RelTime' method. 
#' @param tree Tree(s) to extract date information from. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees support a single topology 
#' }
#' @param method Method for dating trees
#' \itemize{
#'   \item 'reltime' [Default]:  Date tree using the RelTime method, calibrating the root node at 1 (Relative dating; adapted from https://github.com/adamhockenberry/dca-weighting/tree/master/Code/supporting_functions.py)
#'   \item 'reltime_cal':  Use RelTime, and convert relative dates to absolute dates via one node calibration point (Absolute dating)
#'   \item 'chronos': Date tree using 'chronos' function from 'ape', calibrating the root node at 1 (Relative dating)
#'   \item 'chronos_cal': Use 'chronos', and convert relative dates to absolute dates via one or more node calibration points (Absolute dating)
#' }
#' @param calibration_df Dataframe with 4 columns: (1) Taxon 1 (2) Taxon 2 (3) Min divergence time (4) Max divergence time; NOTE: If using 'chronos_cal' dating, multiple calibrartion points are allowed, but only the first calibration point for 'reltime_cal', and the min and max estimates will be averaged together
#' @param iterations If using 'chronos', how many times to estimate the age of each node prior to summarizing [Default: 1000]
#' @return An ultrametric phylo object with branch lengths corresponding to time, or a multiPhylo of such trees.
#' @export

treeDater <- function(tree,method,calibration_df,iterations){
  
  # Ensure tree is present and valid
  if(missing(tree)){
    stop("treeDater requires a rooted phylo object, or a rooted set of multiPhylo trees that all support a common topology.")
  } else if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_all_taxa = TRUE,check_all_equal = TRUE)){
    stop("treeDater requires a rooted phylo or multiPhylo object where all trees have an identical topology")  
  }
  
  # Get tree taxa
  if(Rboretum::isPhylo(tree)){
    tree_taxa <- tree$tip.label
    tree <- c(tree,tree) # Create dummy multiPhylo for simple handling
    tree_count <- 1
  } else{
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    tree_count <- length(tree)
    return_tree <- tree
  }
  
  # Get root clades
  root_clades <- Rboretum::getTreeSplits(tree[[1]]) %>% filter(!is.na(Split_Node))
  root_A <- root_clades$Clade[[1]] %>% semiVector()
  root_B <- root_clades$Mirror_Clade[[1]] %>% semiVector()

  # Get tree dating method
  if(missing(method)){
    method <- 'reltime'
  } else if(!is.character(method)){
    warning("'method' should be 'reltime', 'reltime_cal', 'chronos', or 'chronos_cal'. Proceeding using uncalibrated RelTime")
    method  <- 'reltime'
  } else if(length(method)>1){
    warning("'method' should be 'reltime', 'reltime_cal', 'chronos', or 'chronos_cal'. Proceeding using uncalibrated RelTime")
    method  <- 'reltime'
  } else if(!method %in% c('reltime','reltime_cal','chronos')){
    warning("'method' should be 'reltime', 'reltime_cal', 'chronos', or 'chronos_cal'. Proceeding using uncalibrated RelTime")
    method  <- 'reltime'
  } else if(method %in% c('chronos','chronos_cal')){
    # Set iterations for chronos
    if(missing(iterations)){
      iterations <- 1000
    } else if(!is.numeric(iterations)){
      iterations <- 1000    
    } else if(iterations < 1){
      iterations <- 1000
    }
  } 
  
  # Ensure at least one calibration point if using 'reltime_cal' or 'chronos_cal'
  if(!method %in% c('reltime','chronos')){
    if(missing(calibration_df)){
      stop("If dating trees with 'reltime_cal' or 'chronos_cal', treeDater requires min/max divergence time estimates for at least one node in the tree to convert relative times to absoute times.")
    } else if(!is.data.frame(calibration_df)){
      stop("'calibration_df' should be a data frame.")
    } else if(!ncol(calibration_df)==4){
      stop("'calibration_df' should have 4 columns. (1) Taxon 1 (2) Taxon 2 (3) Min divergence time (4) Max divergence time")
    } else{
      
      # In case data.frame contains factors
      calibration_df <- as.data.frame(calibration_df) %>%
        mutate_if(is.factor, as.character)
      
      # Set colnames
      colnames(calibration_df) <- c('Taxon_1','Taxon_2','Min','Max')
      
      # Ensure all min time estimates are <= max time estimates
      if(!all(calibration_df$Min <= calibration_df$Max)){
        stop("The minimum divergence time estimates for some calibration data in 'calibration_df' are greater than their associated maximum divergence time estimate")
      }
    }
    
    # For RelTime, only use the top calibration point
    if(method == 'reltime_cal' & nrow(calibration_df)>1){
      warning("'reltime_cal' method can only handle a single node calibration. Using first row of 'calibration_df' as the calibrartion node...")
      calibration_df <- calibration_df %>% head(1)
    }
    
    # Check that calibration taxa exist in tree
    cal_taxa <- select(calibration_df,starts_with('Taxon')) %>% unlist() %>% unique()
    
    if(!all(cal_taxa %in% tree_taxa)){
      stop("Calibration data in 'calibration_df' contains information about taxa not present in 'tree'")
    } else{
      calibration_df$Two_Names <- as.character(paste(c(calibration_df$Taxon_1,calibration_df$Taxon_2),collapse = ";"))
      calibration_df <- calibration_df %>% rowwise() %>% mutate(Two_Names = semiSorter(Two_Names))
    }
  } else if(method %in% c('chronos')){
    
  }
  
  # Process chronos trees
  if(method == 'chronos'){
    
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
        return(tree_chronos_export)
      } else{
        return_tree[[i]] <- tree_chronos_export
      }
    }
      return(return_tree)
  } else{ # RelTime
    
    # Get nodes for MRCA for single calibration point
    cal_node <- ape::getMRCA(date_tree,tip=semiVector(calibration_df$Two_Names[[1]]))
    cal_age <- mean(c(calibration_df$Min[[1]],calibration_df$Max[[1]]))
    
    # Process node dating for each tree
    for(i in 1:tree_count){
    
      date_tree <- tree[[i]]
    
      # Set root node labels for Python
      date_tree$node.label <- paste0("Node_", LETTERS[1:length(subtrees(date_tree))])
      date_tree$node.label[[2]] <- 'Root_A' # Second entry in subtrees is the first half of root
      date_tree$node.label[[length(subtrees(date_tree))]] <- 'Root_B' # Last entry in subtrees is the second half of the root

      # Get RelTime tree
      date_reltime <- getRelTimeTree(write.tree(date_tree)) %>% read.tree(text = .) # Pass and retrieve trees through Reticulate via strings
      
      # FIGURE OUT SCALING USING NODE NUMBER OR TAXA
      primate_reltime$node.label <- paste0("Node_", LETTERS[1:length(subtrees(primate_reltime))])
      treePlotter(primate_reltime,branch_length = TRUE,xmax=0.09)
      
      focal_node <- 'Node_A'
      focal_node_age <- 43.35
      focal_node_rel_age <- branching.times(primate_reltime)[[focal_node]]
      absolute_node_estimates <- (branching.times(primate_reltime)/focal_node_rel_age)*focal_node_age
      primate_reltime$node.label <- round(absolute_node_estimates,2)
      reltime_plot <- treePlotter(primate_reltime,branch_length = TRUE,xmax=0.09,node_label_box = TRUE)

    
    }
  }
}