#' Rboretum Tree Dater
#'
#' This function takes a tree(s) with estimated branch lengths, and infers either relative or absolute node ages
#' @param tree Tree(s) from which to extract branch length information. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees share all taxa and support a single topology 
#' }
#' @param method Estimation of node ages via:
#' \itemize{
#'   \item 'reltime_uncal'[Default]: A RelTime-like implementation without node age calibration (Relative dating)
#'   \item 'reltime': A RelTime-like implementation, allowing for relative-to-absolute age conversion via a single node calibration supplied by 'calibration_df' (if more than one calibration is provided, only the first entry will be processed.)
#'   \item 'chronos': Estimate node age using 'chronos' from 'ape', allowing for relative-to-absolute age conversion via multiple node calibrations supplied by 'calibration_df'
#'   }
#' @param calibration_df Required if method is 'reltime' or 'chronos'; Dataframe with 4 columns: (1) Taxon 1 (2) Taxon 2 (3) Min divergence time (4) Max divergence time
#' @param iterations If using 'chronos', how many times to estimate the age of each node prior to summarizing [Default: 1000]
#' @return An ultrametric phylo object with branch lengths corresponding to time, or a multiPhylo of such trees
#' @export

treeDater <- function(tree,method,calibration_df,iterations){

  # Ensure tree is valid
  if(missing(tree)){
    stop("treeDater requires a rooted phylo object, or a rooted set of multiPhylo trees that all support a common topology, sharing all taxa.")
  } else if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_all_taxa = TRUE,check_all_equal = TRUE)){
    stop("treeDater requires a rooted phylo object, or a rooted set of multiPhylo trees that all support a common topology, sharing all taxa.")
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

  # Get method
  if(missing(method)){
    method <- 'reltime_uncal'
  } else if(!is.character(method)){
    warning("'method' should be 'reltime_uncal', 'reltime', or 'chronos'. Proceeding with uncalibrated RelTime...")
    method <- 'reltime_uncal'
  } else if(length(method)!=1){
    warning("'method' should be 'reltime_uncal', 'reltime', or 'chronos'. Proceeding with uncalibrated RelTime...")
    method <- 'reltime_uncal'
  } else if(!method %in% c('reltime_uncal','reltime_new','chronos','reltime_orig')){
    warning("'method' should be 'reltime_uncal', 'reltime', or 'chronos'. Proceeding with uncalibrated RelTime...")
    method <- 'reltime_uncal'
  } else if(!missing(calibration_df) & method=='reltime_uncal'){
    warning("Calibration information provided, but 'reltime_uncal' method was selected. Calibration information can only be used with 'reltime' or 'chronos' method. Ignoring calibration information and returning trees with relative dates...")
  }
  
  # Ensure at least one calibration point if using chronos or reltime
  if(method %in% c('reltime_new','reltime_orig','chronos')){
    
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
    
    # RelTime-like method can only take a single calibration point
    if(method == 'reltime' | method == 'reltime_new'){
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
    # Check that all min divergence times are <= max divergence times
    if(!all(calibration_df$Min <= calibration_df$Max)){
      stop("The minimum divergence time estimates for some calibration data in 'calibration_df' are greater than their associated maximum divergence time estimate")
    }
  }
  
  # Process if method is 'chronos'
  if(method == 'chronos'){
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
      tree_cal <- data.frame(node, age.min, age.max) %>% mutate(soft.bounds=FALSE) %>% `names<-`(c('node','age.min','age.max','soft.bounds'))
      
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
      
      # Taking the median value of a set of iterations can result in a slightly non-ultrametric tree. Sure up loose branches before returning.
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
  }

  # Process if method is 'reltime_uncal'
  if(method == 'reltime_uncal'){
    for(i in 1:tree_count){
      date_tree <- tree[[i]]
      
      # Replace 0 branch lengths with 1E-16 to avoid divide-by-zero errors
      date_tree$edge.length[date_tree$edge.length < 1E-16] <- 1E-16
      
      reltime_tree <- getRelTimeTree(write.tree(date_tree)) %>% read.tree(text=.)
      
      # Minor adjustments to force ultrametric tree
      if(!is.ultrametric(reltime_tree)){
        reltime_tree <- force.ultrametric(reltime_tree,"extend")
      }
      
      if(tree_count == 1){
        return(reltime_tree)
      } else{
        return_tree[[i]] <- reltime_tree
      }
    }
    return(return_tree)
  }
  
  # Rescale RelTime-like trees using calibration point if method is 'reltime'
  if(method == 'reltime_new'){
    for(i in 1:tree_count){
      date_tree <- tree[[i]]
      
      # Replace 0 branch lengths with 1E-16 to avoid divide-by-zero errors
      date_tree$edge.length[date_tree$edge.length < 1E-16] <- 1E-16      
      
      # Get RelTime-like tree
      reltime_tree <- getRelTimeTree_New(write.tree(date_tree)) %>% read.tree(text=.)
      
      # Minor adjustments to force ultrametric tree
      if(!is.ultrametric(reltime_tree)){
        reltime_tree <- force.ultrametric(reltime_tree,"extend")
      }
      
      # Get node ID for calibration node
      node <- ape::getMRCA(reltime_tree,tip=semiVector(calibration_df$Two_Names[[1]])) %>% unlist()      
      
      # Get tip labels
      tip_logical <- reltime_tree$edge[,2] <= length(reltime_tree$tip.label)
      tip_labels <- reltime_tree$edge[tip_logical,2]
      
      # Find shortest distance between calibration node and the nearest tip
      tip_distances <- purrr::map(.x = tip_labels,.f = function(x){dist.nodes(reltime_tree)[node, x]}) %>% unlist()
      shortest_path <- tip_distances[tip_distances==min(tip_distances)][[1]]
      node_cal_ratio <- shortest_path/node_cal_value
      
      # Scale RelTime-like tree to absolute node ages
      reltime_absolute <- reltime_tree
      reltime_absolute$node.label <- NULL
      reltime_absolute$edge.length <- reltime_absolute$edge.length/node_cal_ratio

      if(tree_count == 1){
        return(reltime_tree)
      } else{
        return_tree[[i]] <- reltime_tree
      }
    }
    return(return_tree)
  }
  
  if(method == 'reltime_orig'){
    for(i in 1:tree_count){
      date_tree <- tree[[i]]
      
      # Replace 0 branch lengths with 1E-16 to avoid divide-by-zero errors
      date_tree$edge.length[date_tree$edge.length < 1E-16] <- 1E-16      
      
      # Get RelTime-like tree
      reltime_tree <- getRelTimeTree_Orig(write.tree(date_tree)) %>% read.tree(text=.)
      
      # Minor adjustments to force ultrametric tree
      if(!is.ultrametric(reltime_tree)){
        reltime_tree <- force.ultrametric(reltime_tree,"extend")
      }
      
      # Get node ID for calibration node
      node <- ape::getMRCA(reltime_tree,tip=semiVector(calibration_df$Two_Names[[1]])) %>% unlist()      
      
      # Get tip labels
      tip_logical <- reltime_tree$edge[,2] <= length(reltime_tree$tip.label)
      tip_labels <- reltime_tree$edge[tip_logical,2]
      
      # Find shortest distance between calibration node and the nearest tip
      tip_distances <- purrr::map(.x = tip_labels,.f = function(x){dist.nodes(reltime_tree)[node, x]}) %>% unlist()
      shortest_path <- tip_distances[tip_distances==min(tip_distances)][[1]]
      node_cal_ratio <- shortest_path/node_cal_value
      
      # Scale RelTime-like tree to absolute node ages
      reltime_absolute <- reltime_tree
      reltime_absolute$node.label <- NULL
      reltime_absolute$edge.length <- reltime_absolute$edge.length/node_cal_ratio
      
      if(tree_count == 1){
        return(reltime_tree)
      } else{
        return_tree[[i]] <- reltime_tree
      }
    }
    return(return_tree)
  }
}