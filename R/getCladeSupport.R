#' Rboretum Clade Support Fetcher
#'
#' This function takes a character vector argument of taxon, and returns the total site support count for that clade from an alignment [output from getAlignmentSignal()]
#' @param signal Output table from getAlignmentSignal()
#' @param clade Character vector containing all taxa in clade of interest [Note: All requested taxa must appear in getAlignmentSignal data]
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column [Default: 0]
#' @param include_gap OPTIONAL: TRUE or FALSE; Count sites with gap positions ('-') as part of total support [Default: TRUE]
#' @param include_singleton OPTIONAL: TRUE or FALSE; Count sites with singletons as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: TRUE or FALSE; Count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: TRUE or FALSE; Count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: TRUE or FALSE; Count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: TRUE or FALSE; Count sites with pentallelic variation as part of total support [Default: TRUE]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') as part of total support [Default: FALSE]
#' @param as_root OPTIONAL: Process signal as if it were a root split (biallelic only + no missing taxa allowed) [Default: FALSE]
#' @return integer count of total sites of support for clade in available signal
#' @export

getCladeSupport <- function(signal,clade,max_missing,include_gap,include_singleton,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic,only_gap,as_root){
  
  if(!Rboretum::isAlignmentSignal(signal)){
    stop("'signal' does not appear to be the output from getAlignmentSignal()")
  }
  
  if(missing(max_missing)){
    max_missing <- 0
  }
  
  if(missing(include_gap)){
    include_gap <- TRUE
  } else if (!is.logical(include_gap)){
    include_gap <- TRUE
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
  
  if(missing(only_gap)){
    only_gap <- FALSE
  } else if (!is.logical(only_gap)){
    only_gap <- FALSE
  }
  
  if(missing(as_root)){
    as_root <- FALSE
  } else if(!is.logical(as_root)){
    as_root <- FALSE
  }
  
  if(length(clade)<2){
    stop("Clade must include 2+ taxon IDs")
  }
  
  signal_taxa <- signal %>%
    filter(!is.na(Split_1)) %>% 
    head(1) %>% 
    select(Singleton_Taxa,Non_Base_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5) %>% 
    select_if(~ !any(is.na(.))) %>% 
    unite(col = "Taxa",sep = ";") %>% 
    semiVector() %>% 
    sort()
  
  if(all(clade %in% signal_taxa)){

    informative_patterns <- c('biallelic','triallelic','quadallelic','pentallelic')
    
    if(as_root){
      max_missing <- 0
      informative_patterns <- c('biallelic')
    }
    
    signal <- signal %>% 
      filter(Non_Base_Count <= max_missing) %>%
      filter(Site_Pattern %in% informative_patterns)
    
    if(!include_gap){
      if(only_gap){
        stop("Cannot only use (only_gap) and exclude (include_gap) gap positions.")
      }
      signal <- signal %>%
        filter(Gap == FALSE)
    }
    
    if(!include_singleton){
      signal <- signal %>%
        filter(Singleton==FALSE)
    }
    
    if(!include_biallelic){
      if(as_root){
        stop("Can't run 'as_root' and exclude biallelic sites.")
      } else{
      signal <- signal %>%
        filter(!str_detect(Site_Pattern,'biallelic'))
      }
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
    
    if(only_gap){
      if(!include_gap){
        stop("Cannot only use (only_gap) and exclude (include_gap) gap positions.") 
      } else if(signal %>% filter(Gap==TRUE) %>% nrow() == 0){
        stop("Data contains no gap positions, but 'only_gap' was specified.")
      } else{
        signal <- signal %>%
          filter(Gap==TRUE)
      }
    }
    
    signal_table <- signal %>% 
      select(starts_with('Split_')) %>% 
      unlist() %>% 
      table()
    
    clade_semi <- sort(clade) %>% paste(collapse = ";")
    
    return(as.integer(Rboretum::tableCount(signal_table,clade_semi)))
    
    } else{ stop("Some taxa from 'clade' not present in getAlignmentSignal() data.") }
}