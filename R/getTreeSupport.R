#' Rboretum Alignment Signal Support Mapper
#'
#' This function maps phylogenetic signal from a multiple-sequence alingment onto a rooted phylo or multiPhylo object, or a set of specified clades
#' @param signal Output table from getAlignmentSignal()
#' @param tree OPTIONAL: Tree(s) onto which alignment signal will be mapped. Not considered if 'clade' argument is used. Tree options include:
#' \itemize{
#'   \item A rooted phylo object
#'   \item A rooted multiPhylo object where all trees share 3+ taxa
#' }
#' @param include_root OPTIONAL: If TRUE and using a 'tree',  return alignment support for root clades as well [Default: FALSE, don't include root clades]
#' @param clade OPTIONAL: Character vector of semicolon-separted taxa specifying specific clades of interest (i.e. get support for specific clade(s) rather than a whole tree). Supercedes 'tree' argument
#' @param dataset_name OPTIONAL: Character vector containing a new name for each alignment dataset in 'signal',  [Default: Alignment name from signal dataframe + 'm_<MISSING>']
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column before it is not considered [Default: Taxa Count - 3]
#' @param separate_signal OPTIONAL: If FALSE, return values as the sum of all datasets [Default: TRUE, return results separated by dataset]
#' @param include_singleton OPTIONAL: If FALSE, do not count sites with singletons as part of total support (i.e. parsimony-informative sites) [Default: TRUE]
#' @param include_gap OPTIONAL: If FALSE, count sites with gap positions ('-') as missing data; otherwise, count gaps as valid indel data [Default: TRUE: Gaps are treated as indel signal]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') [Default: FALSE]
#' @param only_biallelic OPTIONAL: If TRUE, only count sites with biiallelic variation as part of total support [Default: FALSE]
#' @param only_triallelic OPTIONAL: If TRUE, only count sites with triallelic variation as part of total support [Default: FALSE]
#' @param only_quadallelic OPTIONAL: If TRUE, only count sites with quadallelic variation as part of total support [Default: FALSE]
#' @param only_pentallelic OPTIONAL: If TRUE, only count sites with pentallelic variation as part of total support [Default: FALSE]
#' @param return_integer OPTIONAL: If TRUE, return the integer support summed across datasets [Default: FALSE, return results as a dataframe]
#' @param return_table OPTIONAL: If TRUE, return entire table of filtered signal counts by clade [FALSE: Return clade support counts]
#' @param existing_support OPTIONAL: Append these results to the output from getTreeSupport() run with the same 'tree' and different alignment options
#' @return A dataframe containing each monophyletic clade in 'tree', along with site support from all alignments in 'signal' as separate columns
#' @export

getTreeSupport <- function(signal,tree,include_root,clade,dataset_name,max_missing,separate_signal,include_singleton,include_gap,only_gap,only_biallelic,only_triallelic,only_quadallelic,only_pentallelic,return_integer,return_table,existing_support){
  
  # Process signal
  if(missing(signal)){
    stop("'getTreeSupport' requires 'signal' arguement.")
  } else if(!Rboretum::isAlignmentSignal(signal)){
    stop("'signal' is not the output of getAlignmentSignal.")
  } else{
    
    # Get taxa from signal
    signal_taxa <- Rboretum::isAlignmentSignal(signal,return_taxa = TRUE)
    
    # Set maximum number of missing taxa allowed
    max_possible_missing <- length(signal_taxa) - 3
    
    # Get max missing
    if(missing(max_missing)){
      max_missing <- max_possible_missing
    } else if(has_error(silent=TRUE,expr=as.integer(max_missing))){
      stop("'max_missing' should be an integer value")
    } else{
      max_missing <- as.integer(max_missing)
      
      # Ensure max_missing doesn't leave fewer than 3 taxa
      if(max_missing > max_possible_missing){
        max_missing <- max_possible_missing
      }
    }
  }
  
  # Check for tree or clade arguments
  if(missing(return_table)){
    if(missing(tree) & missing(clade)){
      stop("Must provide either a 'tree' or 'clade' argument, or ask that 'return_table' = TRUE")
    }
  }
  
  # Process toggle arguments
  if(missing(separate_signal)){
    separate_signal <- TRUE
  } else if(!is.logical(separate_signal)){
    separate_signal <- TRUE
  }
  
  if(missing(return_integer)){
    return_integer <- FALSE
  } else if(!is.logical(return_integer)){
    return_integer <- FALSE
  }
  
  if(missing(return_table)){
    return_table <- FALSE
  } else if(!is.logical(return_table)){
    return_table <- FALSE
  }
  
  if(missing(include_root)){
    include_root <- FALSE
  } else if(!is.logical(include_root)){
    include_root <- FALSE
  }
  
  if(missing(include_singleton)){
    include_singleton <- TRUE
  } else if (!is.logical(include_singleton)){
    include_singleton <- TRUE
  }
  
  if(missing(include_gap)){
    include_gap <- TRUE
  } else if (!is.logical(include_gap)){
    include_gap <- TRUE
  }
  
  if(missing(only_gap)){
    only_gap <- FALSE
  } else if (!is.logical(only_gap)){
    only_gap <- FALSE
  }
  
  if(missing(only_biallelic)){
    only_biallelic <- FALSE
  } else if (!is.logical(only_biallelic)){
    only_biallelic <- FALSE
  }
  
  if(missing(only_triallelic)){
    only_triallelic <- FALSE
  } else if (!is.logical(only_triallelic)){
    only_triallelic <- FALSE
  }
  
  if(missing(only_quadallelic)){
    only_quadallelic <- FALSE
  } else if (!is.logical(only_quadallelic)){
    only_quadallelic <- FALSE
  }
  
  if(missing(only_pentallelic)){
    only_pentallelic <- FALSE
  } else if (!is.logical(only_pentallelic)){
    only_pentallelic <- FALSE
  }
  
  # Get clade information (unless return_table = TRUE)
  if(!return_table){
    
    # Process clade instead of tree
    if(!missing(clade)){
      if(!Rboretum::semiChecker(clade)){
        stop("'clade' argument should be a charcter vector of taxa separated by semicolons...")
      } else{
        test_clade <- semiSorter(clade)
      }
    } else{ # If no 'clade' argument, process 'tree'
      
      if(Rboretum::isPhylo(tree,check_rooted = TRUE)){
        
        if(include_root){
          test_clade <- Rboretum::getTreeClades(tree,include_root = TRUE)
        } else{
          test_clade <- Rboretum::getTreeClades(tree)
        }
      } else if(Rboretum::isMultiPhylo(check_rooted = TRUE,check_three_taxa = TRUE)){
        
        # Ensure all trees have the same taxa
        if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
          tree <- treeTrimmer(tree)
        }
        
        # Autoname multiPhylo if necessary
        if(!Rboretum::isMultiPhylo(check_named = TRUE)){
          tree <- treeNamer(tree)
        }
        
        # Get clades from tree
        if(include_root){
          test_clade <- Rboretum::getTreeClades(tree,include_root = TRUE)
        } else{
          test_clade <- Rboretum::getTreeClades(tree)
        }
      } else{
        stop("'tree' should be a rooted phylo, or a rooted multiPhylo where all trees share 3+ taxa...")
      }
    }
    
    # Get test taxa
    test_taxa <- purrr::map(.x=test_clade,.f=function(x){semiVector(x)}) %>% unlist() %>% as.character() %>% unique() %>% naturalsort()
    
    # Validate signal and get signal taxa
    if(!all(test_taxa %in% signal_taxa)){
      stop("'signal' is either not the output from getAlignmentSignal() for the supplied clades/trees")
    }
  }
  
  # Extract informative sites
  informative_patterns <- c('biallelic','triallelic','quadallelic','pentallelic')
  signal <- signal %>%
    filter(Non_Base_Count <= max_missing) %>%
    filter(Site_Pattern %in% informative_patterns)
  
  # Ensure filtering has left some data
  if(nrow(signal)==0){
    stop("No data fits the filtering criteria.")
  }
  
  if(!include_singleton){
    signal <- signal %>%
      filter(is.na(Singleton_Taxa))
  }
  
  if(!include_gap){
    signal <- signal %>%
      filter(!str_detect("-",All_Site_Bases))
  }
  
  if(only_gap){
    signal <- signal %>%
      filter(str_detect("-",All_Site_Bases))
    include_gap <- TRUE
  }
  
  if(only_biallelic){
    signal <- signal %>%
      filter(Site_Pattern=='biallelic')
  }
  
  if(only_triallelic){
    signal <- signal %>%
      filter(Site_Pattern=='triallelic')
  }
  
  if(only_quadallelic){
    signal <- signal %>%
      filter(Site_Pattern=='quadallelic')
  }
  
  if(only_pentallelic){
    signal <- signal %>%
      filter(Site_Pattern=='pentallelic')
  }
  
  # Ensure filtering has left some data
  if(nrow(signal)==0){
    stop("No data fits the filtering criteria.")
  } else{
    final_signal_name <- unique(signal$Alignment_Name)
    final_signal_count <- length(final_signal_name)
  }
  
  if(!separate_signal){
    default_name <- paste(c('Total_m',max_missing),collapse = '')
  } else{
    # Generate default names based on signal object and number of missing taxa allowed (replaced by supplying dataset_name vector)
    default_name <- purrr::map(.x=final_signal_name,.f=function(x){paste(c(x,'_m',max_missing),collapse = '')}) %>% unlist() %>% as.character()
  }
  
  # Set alignment names to defaults if necessary
  if(missing(dataset_name)){
    dataset_name <- default_name
  } else if(!is.character(dataset_name)){
    dataset_name <- default_name
  } else if(!separate_signal & length(dataset_name)>1){
    print("'separate_signal' disabled, but more than one 'dataset_name' provided. Using default name...")
    dataset_name <- default_name
  } else if(separate_signal & length(dataset_name) != length(signal_name)){ # If number of alignments in 'signal' are different from the number of new names provided, use default names
    print("'signal' contains a different number of alignments than names provided by 'dataset_name'. Using default names...")
    dataset_name <- default_name
  } else{
    dataset_name <- dataset_name[final_signal_name %in% signal_name]
  }
  
  # Generate support tables
  if(!separate_signal | final_signal_count==1){ # If return results as a summation, or if only one alignment is present...
    support_table <- signal %>% select(starts_with('Split_')) %>% unlist() %>% table()
    
  } else{ # If splitting results up by dataset...
    
    support_table <- list()
    
    for(i in 1:final_signal_count){
      support_table[[i]] <- signal %>%
        filter(Alignment_Name == final_signal_name[i]) %>%
        select(starts_with('Split_')) %>%
        unlist() %>% table()
    }
    
    names(support_table) <- dataset_name
  }
  
  if(return_table){
    return(support_table)
  }
  
  clade_df <- data.frame(Clade = test_clade)
  
  # Add to existing support?
  if(missing(existing_support)){
    add_support <- FALSE
  } else if(!Rboretum::isTreeSupport(existing_support,test_clade,partial = TRUE)){
    print("'existing_support' is either not the output from getTreeSupport(), or does not contain clade information as requested from 'tree' or 'clade'. Returning unappended...")
    add_support <- FALSE
  } else{
    if(any(names(existing_support) %in% dataset_name)){
      print(names(existing_support)[names(existing_support) %in% dataset_name])
      print("'existing_support' already contains columns with the column names above. Cannot add a column with an identical names. Returning unappended...")
      add_support <- FALSE
    } else{
      add_support <- TRUE
    }
  }
  # Generate support counts
  if(!separate_signal | final_signal_count==1){ # If returning results as a summation, or if only one alignment is present...
    
    clade_support <- purrr::map(.x=test_clade,.f=function(x){Rboretum::tableCount(support_table,x)}) %>% unlist() %>% as.integer()
    
    if(return_integer){
      return(clade_support)
    }
    
    clade_df <- clade_df %>%
      mutate(Support = clade_support) %>%
      rename(!!dataset_name := Support)
    
  } else{ # If splitting results up by dataset...
    
    clade_support <- purrr::map(.x=dataset_name,.f=function(x){lapply(test_clade,function(y) Rboretum::tableCount(support_table[[x]],as.character(y))) %>% unlist() %>% as.integer()})
    names(clade_support) <- dataset_name
    
    if(return_integer){
      return(clade_support)
    }
    
    for(i in 1:final_signal_count){
      clade_df[,(i+1)] <- as.integer(clade_support[[i]])
    }
    
    names(clade_df) <- c('Clade',dataset_name)
  }
  
  if(add_support){
    existing_support <- existing_support %>%
      left_join(clade_df,by='Clade')
    print("Data added successfully...")
    return(existing_support)
  } else{
    return(clade_df)
  }
}