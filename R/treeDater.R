#' Rboretum Tree Dater
#'
#' This function takes a tree(s) with estimated branch lengths, and returns ultrametric chronograms generated using the 'chronos' function in 'ape'
#' @param tree Tree(s) to extract date information from. Options include:
#' \itemize{
#'   \item A single, rooted phylo object; or,
#'   \item A rooted multiPhylo object where all trees support a single topology 
#' }
#' @param calibration_df A 4-column dataframe/tibble with calibration information: (1) Taxon 1 (2) Taxon 2 [to get MRCA] (3) Min divergence time (4) Max divergence time; Multiple calibration points are allowed, so long as nodes are only calibrated once [Supercedes 'taxa' argument]
#' @param taxa A character vector (or semicolon-delimited set) of taxon IDs from which to find the MRCA and calibrate [One calibration point allowed; Superceded by 'calibration_df']
#' @param min_max If using 'node', a two-element numeric vector [e.g. c(50,75)] that provides the minimum and maximum age estimates for the focal calibration node [Superceded by 'calibration_df'; min <= max]
#' @param iterations How many times to estimate the age of each node prior to summarizing [Default: 1000]
#' @return An ultrametric phylo object with branch lengths corresponding to time, or a multiPhylo of such trees.
#' @export

treeDater <- function(tree,calibration_df,taxa,min_max,iterations){
  
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

  # Set iterations for chronos
  if(missing(iterations)){
    iterations <- 1000
  } else if(!is.numeric(iterations)){
    iterations <- 1000    
  } else if(iterations < 1){
    iterations <- 1000
  }

  # Get calibration information
  if(missing(calibration_df) & missing(taxa)){
    stop("treeDater requires min/max divergence time estimates for at least one node in the tree to convert relative times to absoute times. Use 'calibration_df' or 'taxa'/'min_max'...")
  }
  
  if(!missing(calibration_df)){
    
    if(!is.data.frame(calibration_df)){
      stop("'calibration_df' should be a data frame.")
    } else if(!ncol(calibration_df)==4){
      stop("'calibration_df' should have 4 columns. (1) Taxon 1 (2) Taxon 2 (3) Min divergence time (4) Max divergence time")
    } else{
      
      # In case data.frame contains factors
      calibration_df <- as.data.frame(calibration_df) %>%
        mutate_if(is.factor, as.character)
      
      # Set colnames
      colnames(calibration_df) <- c('Taxon_1','Taxon_2','Min','Max')
      
      # Check coltypes
      if(!is.character(calibration_df$Taxon_1) | !is.character(calibration_df$Taxon_2) | !is.numeric(calibration_df$Min) | !is.numeric(calibration_df$Max)){
        stop("'calibration_df' should have 4 columns. (1) Taxon 1 (2) Taxon 2 (3) Min divergence time (4) Max divergence time")
      }
      
      # Ensure all min time estimates are <= max time estimates
      if(!all(calibration_df$Min <= calibration_df$Max)){
        stop("The minimum divergence time estimates for some calibration data in 'calibration_df' are greater than their associated maximum divergence time estimate")
      }
    }
    
    # Check that calibration taxa exist in tree
    cal_taxa <- select(calibration_df,starts_with('Taxon')) %>% unlist() %>% unique()
    
    if(any(calibration_df$Taxon_1 == calibration_df$Taxon_2)){
      stop("Taxon_1 and Taxon_2 from 'calibration_df' cannot be identical...")
    } else if(!all(cal_taxa %in% tree_taxa)){
      stop("Calibration data in 'calibration_df' contains information about taxa not present in 'tree'")
    } else{
      calibration_df <- calibration_df %>% rowwise() %>% mutate(Taxa = semiSorter(c(Taxon_1,Taxon_2))) %>% ungroup() %>% select(Taxa,Min,Max)
    }
  } 
  
  # If no calibration_df is provided...
  if(missing(calibration_df)){
    
    if(missing(min_max)){ 
      stop("If using 'taxa' to specify calibration node, you must also supply min and max bounds for node date calibration via min_max")
    } else if(length(min_max)!=2){
      stop("'min_max' should be a two-element numeric vector")
    } else if(!is.numeric(min_max)){
      stop("'min_max' should be a two-element numeric vector")
    } else if(!(min_max[[1]]<=min_max[[2]])){
      stop("min of 'min_max' is > max...")
    } 
    
    if(!is.character(taxa)){
      stop("'taxa' should be a character vector of taxon IDs (or a semicolon-delimited set)...")
    } else if(!(length(taxa) >= 1)){
      stop("'taxa' cannot be 0-length...")
    } else if(length(taxa)==1){
      if(!semiChecker(taxa)){
        stop("Only one taxon ID provided, and it doesn't appear to be semicolon-delimited...")
      } else{
        cal_taxa <- unique(semiVector(taxa))
      }
    } else if(any(str_detect(taxa,";"))){
        stop("Multiple 'taxa' provided, but semicolons detected?")
      } else{
        cal_taxa <- unique(taxa)
      }

    if(!all(cal_taxa %in% tree_taxa)){
      stop("Calibration data provided by 'taxa' contains information about taxa not present in 'tree'")
    } else{
      calibration_df <- tibble(Taxa = semiSorter(cal_taxa),Min=min_max[[1]],Max=min_max[[2]])
    }
  }
  
  # Process chronos trees
  for(i in 1:tree_count){
    
    date_tree <- tree[[i]]
    
    # Get nodes for MRCA for each calibration point
    node <- purrr::map(.x=calibration_df$Taxa,.f=function(x){ape::getMRCA(date_tree,tip=semiVector(x))}) %>% unlist()
    
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
    
    tree_chronos_export$edge.length <- tree_median_edge
    if(!is.ultrametric(tree_chronos_export)){
      tree_chronos_export <- force.ultrametric(tree_chronos_export,"extend")
    }
    
    if(tree_count == 1){
      return_tree <- tree_chronos_export
    } else{
      return_tree[[i]] <- tree_chronos_export
    }
  }
    return(return_tree)
}