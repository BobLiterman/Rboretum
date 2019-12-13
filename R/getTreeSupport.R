#' Rboretum Alignment Signal Support Mapper
#'
#' This function maps phylogenetic signal from a multiple-sequence alingment onto a rooted phylo or multiPhylo object
#' @param signal Output table from getAlignmentSignal()
#' @param tree OPTIONAL: Tree(s) onto which alignment signal will be mapped. Not considered if 'clade' argument is used. Tree options include:
#' \itemize{
#'   \item A rooted phylo object
#'   \item A rooted, named multiPhylo object where all trees share 3+ taxa
#' }
#' @param clade OPTIONAL: Character vector including specific clades of interest as semicolon-separted taxon (i.e. get support for a specific clade set rather than a whole tree). Supercedes 'tree' argument
#' @param separate_signal OPTIONAL: If FALSE, return values as the sum of all datasets [Default: TRUE, return results separated by dataset]
#' @param return_integer OPTIONAL: If TRUE, and a specific set of clades are queried, return the integer support summed across datasets [Default: FALSE, return results as a dataframe]
#' @param include_root OPTIONAL: If TRUE and using a 'tree',  return alignment support for root clades as well [Default: FALSE, don't include root clades]
#' @param dataset_name OPTIONAL: Character vector containing a new name for each alignment dataset in 'signal',  [Default: Alignment name from signal dataframe + 'm_<MISSING>']
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column before it is not considered [Default: 0, no missing taxa allowed]
#' @param include_gap OPTIONAL: If TRUE, count sites with gap positions ('-') as informative signal; otherwise, count gaps as missing data [Default: FALSE: Gaps are treated as missing]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') [Default: FALSE]
#' @param include_singleton OPTIONAL: If TRUE, count sites with singletons as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: If TRUE, count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: If TRUE, count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: If TRUE, count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: If TRUE, count sites with pentallelic variation as part of total support [Default: TRUE]
#' @param return_table OPTIONAL: If TRUE, return entire table of filtered signal counts by clade [FALSE: Return clade support counts]
#' @param existing_support OPTIONAL: Append these results to the output from getTreeSupport() run with the same 'tree' and different alignment options
#' @return A dataframe containing each monophyletic clade in 'tree', along with site support from all alignments in 'signal' as separate columns
#' @export

getTreeSupport <- function(signal,tree,clade,separate_signal,return_integer,include_root,dataset_name,max_missing,include_gap,only_gap,include_singleton,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic,return_table,existing_support){
  
  # Check for tree or clade arguments
  if(missing(tree) & missing(clade) & missing(return_table)){
    stop("Must provide either a 'tree' or 'clade' argument, or ask that 'return_table' = TRUE")
  }
  
  # Validate signal and get signal taxa
  if(!Rboretum::isAlignmentSignal(signal)){
    stop("'signal' is either not the output from getAlignmentSignal()")
  } else{
    signal_taxa <- Rboretum::isAlignmentSignal(signal,return_taxa = TRUE)
    signal_name <- unique(signal$Alignment_Name)
  }
  
  # Set maximum number of missing taxa allowed
  max_possible_missing <- length(signal_taxa) - 3
  
  if(missing(max_missing)){
    max_missing <- 0
  } else if(has_error(silent=TRUE,expr=as.integer(max_missing))){
    stop("'max_missing' should be an integer value")
  } else{
    max_missing <- as.integer(max_missing)
  }
  
  if(max_missing > max_possible_missing){
    max_missing <- max_possible_missing
  }
  
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
  
  if(return_integer){
    separate_signal <- FALSE
  }
  
  if(missing(return_table)){
    return_table <- FALSE
  } else if(!is.logical(return_table)){
    return_table <- FALSE
  }
  
  if(!separate_signal){
    default_name <- paste(c('Total_m',max_missing),collapse = '')
  } else{
    # Generate default names based on signal object and number of missing taxa allowed (replaced by supplying dataset_name vector)
    default_name <- purrr::map(.x=signal_name,.f=function(x){paste(c(x,'_m',max_missing),collapse = '')}) %>% unlist() %>% as.character()
  }
  
  # Set alignment names to defaults if necessary
  if(missing(dataset_name)){
    dataset_name <- default_name
  } else if(!is.character(dataset_name)){
    dataset_name <- default_name
  } else if(!separate_signal & length(dataset_name)!=1){
    stop("Requested that data not be separated, yet provided multiple dataset names")
  } else if(separate_signal & length(dataset_name) != length(signal_name)){ # If number of alignments in 'signal' are different from the number of new names provided, use default names
    stop("'signal' contains a different number of alignments than names provided by 'dataset_name'")
  }
  
  if(!return_table){
    
    if(!missing(clade)){
      if(!is.character(clade)){
        stop("'clade' argument should be a charcter vector of semicolon-separated taxa")
      } else if(!purrr::map(.x=clade,.f=function(x){str_detect(x,";")}) %>% unlist() %>% all()){
        stop("'clade' argument should be a charcter vector of semicolon-separated taxa")
      } else{
        test_taxa <- purrr::map(.x=clade,.f=function(x){semiVector(x)}) %>% unlist() %>% as.character() %>% unique() %>% sort()
        test_clade <- purrr::map(.x=clade,.f=function(x){semiSorter(x)}) %>% unlist() %>% as.character() # Ensure clades are sorted internally
        if(!all(test_taxa %in% signal_taxa)){
          stop("'clade' contains taxon not present in signal")
        }
      }
      
    } else{ # If no 'clade' argument, process 'tree'
      
      if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_named = TRUE,check_rooted = TRUE, check_three_taxa = TRUE)){
        stop("'tree' must be a rooted phylo object, or a named, rooted multiPhylo for getTreeSupport")
      } else if(!Rboretum::isAlignmentSignal(signal,tree)){
        stop("'signal' is either not the output from getAlignmentSignal(), or does not contain the same taxa as 'tree'")
      }
      
      if(missing(include_root)){
        include_root <- FALSE
      } else if(!is.logical(include_root)){
        include_root <- FALSE
      }
      
      # Get clades from tree
      if(include_root){
        test_clade <- Rboretum::getTreeClades(tree,include_root = TRUE)
      } else{
        test_clade <- Rboretum::getTreeClades(tree)
      }
    }
    
    if(missing(existing_support)){
      add_support <- FALSE
    } else if(!Rboretum::isTreeSupport(existing_support,test_clade)){
      stop("'existing_support' is either not the output from getTreeSupport(), or does not contain identical clade information as requested from 'tree' or 'clade'")
    } else{
      
      if(any(names(existing_support) %in% dataset_name)){
        print(names(existing_support)[names(existing_support) %in% dataset_name])
        stop("'existing_support' already contains columns with the column names above. Cannot add a column with an identical names.")
      } else{
        add_support <- TRUE
      }
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
  
  informative_patterns <- c('biallelic','triallelic','quadallelic','pentallelic')
  
  signal <- signal %>% 
    filter(Non_Base_Count <= max_missing) %>%
    filter(Site_Pattern %in% informative_patterns)
  
  if(only_gap){
    signal <- signal %>%
      filter(Gap==TRUE)
    include_gap <- TRUE
  }
  
  if(!include_gap){
    signal <- signal %>%
      filter(Gap==FALSE)
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
  
  # Ensure filtering has left some data
  if(nrow(signal)==0){
    stop("No data fits the filtering criteria.")
  }
  
  surviving_alignments <- unique(signal$Alignment_Name)
  # Generate support tables
  if(!separate_signal | length(signal_name)==1){ # If return results as a summation, or if only one alignment is present...
    
    support_table <- signal %>% select(starts_with('Split_')) %>% unlist() %>% table()
    
  } else{ # If splitting results up by dataset...
    
    support_table <-purrr::map(.x = 1:length(signal_name),.f=function(x){ifelse(signal_name[x] %in% surviving_alignments,signal %>% filter(Alignment_Name == signal_name[x]) %>% select(starts_with('Split_')) %>% unlist() %>% table(),table("RBORTEUM_DUMMY"))})
    names(support_table) <- dataset_name
    
  }
  
  if(return_table){
    return(support_table)
  }
  
  clade_df <- data.frame(Clade = test_clade) %>% mutate(Clade = as.character(Clade))
  
  # Generate support counts
  if(!separate_signal | length(signal_name)==1){ # If returning results as a summation, or if only one alignment is present...
    
    clade_support <- purrr::map(.x=test_clade,.f=function(x){Rboretum::tableCount(support_table,x)}) %>% unlist() %>% as.integer()
    
    if(return_integer){
      return(clade_support)
    }
    
    clade_df <- clade_df %>%
      mutate(Support = clade_support) %>%
      rename(!!dataset_name := Support)
    
  } else{ # If splitting results up by dataset...
    
    clade_support <- purrr::map(.x=dataset_name,.f=function(x){lapply(test_clade,function(y) Rboretum::tableCount(support_table[[x]],as.chararacter(y))) %>% unlist() %>% as.integer()})
    names(clade_support) <- dataset_name
    
    if(return_integer){
      return(clade_support)
    }
    
    for(i in 1:(length(dataset_name))){
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