#' Rboretum Alignment Signal Support Mapper
#'
#' This function maps phylogenetic signal from a multiple-sequence alingment onto a rooted phylo or multiPhylo object
#' @param tree Tree(s) onto which alignment signal will be mapped. Options include:
#' \itemize{
#'   \item A rooted phylo object
#'   \item A rooted, named multiPhylo object where all trees share 3+ taxa
#' }
#' @param signal Output table from getAlignmentSignal() run with taxa against 'tree'
#' @param alignment_name OPTIONAL: Column name(s) for data being added [Default: Alignment name from signal dataframe + 'm_<MISSING>']
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column before it is not considered [Default: 0, no missing taxa allowed]
#' @param include_gap OPTIONAL: If TRUE, count sites with gap positions ('-') as informative signal; otherwise, count gaps as missing data [Default: FALSE: Gaps are treated as missing]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') [Default: FALSE]
#' @param include_singleton OPTIONAL: If TRUE, count sites with singletons as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: If TRUE, count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: If TRUE, count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: If TRUE, count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: If TRUE, count sites with pentallelic variation as part of total support [Default: TRUE]
#' @param existing_support OPTIONAL: Append these results to the output from getTreeSupport() run with the same 'tree' and different alignment options
#' @return The same split table from getTreeSplits(tree), but with (1) no root row, and (2) a support column for the specfied alignment/missing combination
#' @export

getTreeSupport <- function(tree,signal,alignment_name,max_missing,include_gap,only_gap,include_singleton,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic,existing_support){

  # Check if tree and signal are valid and compatible
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_named = TRUE,check_rooted = TRUE, check_three_taxa = TRUE)){
    stop("'tree' must be a rooted phylo object, or a named, rooted multiPhylo for getTreeSupport")
  } else if(!Rboretum::isAlignmentSignal(signal,tree)){
    stop("'signal' is either not the output from getAlignmentSignal(), or does not contain the same taxa as 'tree'")
  } else if(Rboretum::isPhylo(tree)){
    tree_count <- 1
    tree_taxa <- sort(tree$tip.label) # Get tree taxa
  } else{ # If tree is a multiPhylo
    
    tree_taxa <- Rboretum::getSharedTaxa(tree)
    
    # Reduce 'tree' to unique topologies
    if(Rboretum::isMultiPhylo(tree,check_all_equal=TRUE)){
      tree_count <- 1
      tree <- as.phylo(tree[[1]])
    } else if(Rboretum::isMultiPhylo(tree,check_all_unique = TRUE)){
      tree_count <- length(tree)
      tree_names <- names(tree)
    } else if(Rboretum::isMultiPhylo(tree,check_some_equal = TRUE)){
      tree <- Rboretum::getUniqueTopologies(tree)
      tree_count <- length(tree)
      tree_names <- names(tree)
    }
  }
  
  # Set maximum number of missing taxa allowed
  max_possible_missing <- length(tree_taxa) - 3
  
  if(missing(max_missing)){
    max_missing <- 0
  } else if(!is.integer(max_missing)){
    stop("'max_missing' must be an integer")
  } else{
    if(max_missing > max_possible_missing){
      max_missing <- max_possible_missing
    }
  }
  
  if(missing(include_gap)){
    include_gap <- FALSE
  } else if (!is.logical(include_gap)){
    include_gap <- FALSE
  }
  
  if(missing(only_gap)){
    only_gap <- FALSE
  } else if (!is.logical(only_gap)){
    only_gap <- FALSE
  }
  
  if(missing(include_singleton)){
    include_singleton <- TRUE
  } else if (!is.logical(include_singleton)){
    include_singleton <- TRUE
  }
  
  if(missing(include_biallelic)){
    include_biallelic <- TRUE
  } else if (!is.logical(include_biallelic)){
    include_biallelic <- TRUE
  }
  
  if(missing(include_triallelic)){
    include_triallelic <- TRUE
  } else if (!is.logical(include_triallelic)){
    include_triallelic <- TRUE
  }
  
  if(missing(include_quadallelic)){
    include_quadallelic <- TRUE
  } else if (!is.logical(include_quadallelic)){
    include_quadallelic <- TRUE
  }
  
  if(missing(include_pentallelic)){
    include_pentallelic <- TRUE
  } else if (!is.logical(include_pentallelic)){
    include_pentallelic <- TRUE
  }
  
  if(missing(existing_support)){
    add_support <- FALSE
  } else if(!Rboretum::isTreeSupport(existing_support,tree)){
    stop("'existing_support' is either not the output from getTreeSupport(), or does not contain split information from 'tree'")
  } else{
    add_support <- TRUE
  }

  informative_patterns <- c('biallelic','triallelic','quadallelic','pentallelic')
  
  signal <- signal %>% 
    filter(Non_Base_Count <= max_missing) %>%
    filter(Site_Pattern %in% informative_patterns)
  
  if(only_gap){
    signal <- signal %>%
      filter(Gap==TRUE)
  }
  
  if(!include_singleton){
    signal <- signal %>%
      filter(Singleton==FALSE)
  }
  
  if(!include_biallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'biallelic'))
  }
  
  if(!include_triallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'triallelic'))
  }
  
  if(!include_quadallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'quadallelic'))
  }
  
  if(!include_pentallelic){
    signal <- signal %>%
      filter(!str_detect(Site_Pattern,'pentallelic'))
  }
  
  # Figure out alignment count after filtering
  if(nrow(signal)==0){
    stop("No data fits the filtering criteria.")
  }
  
  raw_alignment_name <- unique(signal$Alignment_Name)
  alignment_count <- length(raw_alignment_name)

  if(alignment_count == 1){
    default_name <- paste(c(raw_alignment_name,"_m",max_missing),collapse = '')
  } else{
    default_name <- purrr::map(.x = raw_alignment_name,.f = function(x){paste(c(x,"_m",max_missing),collapse = '')}) %>% unlist() %>% as.character()
  }
  
  # Set alignment names to defaults if necessary
  if(missing(alignment_name)){
    alignment_name <- default_name
  } else if(!is.character(alignment_name)){
    alignment_name <- default_name
  } else if(length(alignment_name) != alignment_count){
    alignment_name <- default_name
  }

  if(tree_count == 1){
    support_df <- Rboretum::getTreeSupport_Worker(tree,signal,alignment_name)
    
    if(add_support){
      
      # Get old alignment names
      old_names <- names(exisiting_support)[4:ncol(existing_support)]
      
      # Can't have duplicate IDs
      if(any(alignment_name %in% old_names)){
        print(alignment_name[alignment_name %in% old_names])
        print("Columns above already in existing_support, returning unappendend results")
        return(support_df)
      } else{
        support_df <- left_join(existing_support,support_df,by=c('Clade','Mirror_Clade','Split_Node'))
        print(paste(c('Added results from:',paste(alignment_name,collapse = ";")),collapse = " "))
        return(support_df)
      }
    } else{
      return(support_df)
    }
  } else{ # If result is for multiple trees, return a named list
    support_list <- purrr::map(.x=tree,.f = function(x){Rboretum::getTreeSupport_Worker(x,signal,alignment_name)})
    names(support_list) <- tree_names
    
    if(add_support){
      # Get old alignment columns
      old_names <- names(existing_support[[1]])[4:ncol(existing_support[[1]])]
      
      # Can't have duplicate IDs
      if(any(alignment_name %in% old_names)){
        print(alignment_name[alignment_name %in% old_names])
        print("Columns above already in existing_support, returning unappendend results")
        return(support_list)
      } else{
        appended_list <- purrr::map(.x=support_list,.f=function(x){left_join(existing_support,x,by=c('Clade','Mirror_Clade','Split_Node'))})
        names(appended_list) <- names(support_list)
        print(paste(c('Added results from:',paste(alignment_name,collapse = ";")),collapse = " "))
        return(appended_list)
      }
    } else{
      return(support_list)
    }
  }
}