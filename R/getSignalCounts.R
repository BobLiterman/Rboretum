#' Rboretum Alignment Signal Support Tabler
#'
#' This function extracts phylogenetic signal from a multiple-sequence alignment and returns a table where the names are semicolon-separated supported clades
#' @param signal Output table from getAlignmentSignal()
#' @param clade OPTIONAL: Character vector containing all members (2+) of clade for which support is requested. If provided, return total intger support for this set of taxa, rather than the entire table. 
#' @param separate OPTIONAL: If true, separate results from different alignments (as coded by the Alignment_Name column) [Default: FALSE, return a single table/count]
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column before it is not considered [Default: 0, no missing taxa allowed]
#' @param include_gap OPTIONAL: If TRUE, count sites with gap positions ('-') as informative signal; otherwise, count gaps as missing data [Default: FALSE: Gaps are treated as missing]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') [Default: FALSE]
#' @param include_singleton OPTIONAL: If TRUE, count sites with singletons as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: If TRUE, count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: If TRUE, count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: If TRUE, count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: If TRUE, count sites with pentallelic variation as part of total support [Default: TRUE]
#' @return Table of counts for clade support among  all sites in multiple-sequence alignment; or support for a specific clade if requested.
#' @export

getSignalCounts <- function(signal,clade,separate,max_missing,include_gap,only_gap,include_singleton,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic){

  # Check if signal is valid
  if(!Rboretum::isAlignmentSignal(signal)){
    stop("'signal' is either not the output from getAlignmentSignal()")
  } else{
    signal_taxa <- Rboretum::isAlignmentSignal(signal,return_taxa = TRUE)
  }
  
  # Process clade information if present
  if(missing(clade)){
    return_clade <- FALSE
  } else if(!is.character(clade)){
    stop("'clade' should be a character vector of 2+ taxa")
  } else if(length(clade)==1){
    if(str_detect(clade,";")){
      clade <- sort(semiVector(clade))
      return_clade <- TRUE
    } else{
      stop("'clade' must contain 2+ taxa")
    }
  } else{
    clade <- sort(clade)
    return_clade <- TRUE
  }
  
  if(return_clade){
    if(!all(clade %in% signal_taxa)){
      print("Signal Taxa:")
      print(signal_taxa)
      print("Requested Clade:")
      print(clade)
      stop("Some/all species from 'clade' species not present in alignment signal")
    } else{
      clade <- vectorSemi(clade)
    }
  }
  
  # Separate results by alignment, or return total?
  if(missing(separate)){
    separate <- FALSE
  } else if(!is.logical(separate)){
    separate <- FALSE
  }
  
  # Set maximum number of missing taxa allowed
  max_possible_missing <- length(signal_taxa) - 3
  
  if(missing(max_missing)){
    max_missing <- 0
  } else if(max_missing == 0){
    max_missing <- as.integer(0)
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
  
  # Figure out alignment count after filtering
  if(nrow(signal)==0){
    stop("No data fits the filtering criteria.")
  }
  
  raw_alignment_name <- unique(signal$Alignment_Name)
  alignment_count <- length(raw_alignment_name)
  
  if(alignment_count == 1 | !separate){
    support_table <- signal %>% select(starts_with('Split_')) %>% unlist() %>% table()
    
    if(return_clade){
      return(Rboretum::tableCount(support_table,clade)) 
    } else{
      return(support_table)
    }
  } else if(alignment_count > 1 & separate){
    support_table <- purrr::map(.x = raw_alignment_name, .f = function(x){signal %>% filter(Alignment_Name == x) %>% select(starts_with('Split_')) %>% unlist() %>% table()})
    names(support_table) <- raw_alignment_name
    
    if(return_clade){
      support_counts <- purrr::map(.x = raw_alignment_name, .f = function(x){Rboretum::tableCount(support_table[[x]],clade)}) %>% unlist()
      support_df <- data.frame(Alignment=as.character(raw_alignment_name),Support=as.integer(support_counts))
      return(support_df)
    } else{
      return(support_table)
    }
  }
}