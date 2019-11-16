#' Rboretum Clade Support Fetcher
#'
#' This function takes a character vector argument of taxon, and returns the total site support count for that clade from an alignment [output from alignment.signal()]
#' @param signal Output table from alignment.signal()
#' @param clade Character vector containing all taxa in clade of interest [Note: All requested taxa must appear in alignment signal data]
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column [Default: 0]
#' @param include_gap OPTIONAL: TRUE or FALSE; Count sites with gap positions ('-') as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: TRUE or FALSE; Count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: TRUE or FALSE; Count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: TRUE or FALSE; Count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: TRUE or FALSE; Count sites with pentallelic variation as part of total support [Default: TRUE]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') as part of total support [Default: FALSE]
#' @return integer count of total sites of support for clade in available signal
#' @export
#' @examples
#' clade.support <- (signal,clade)
#' 
clade.support <- function(signal,clade,max_missing,include_gap,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic,only_gap){
  
  if(missing(max_missing)){
    max_missing <- 0
  }
  
  if(missing(include_gap)){
    include_gap <- TRUE
  } else if (!is.logical(include_gap)){
    include_gap <- TRUE
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
  
  if(length(clade)<2){
    stop("Clade must include 2+ taxon IDs")
  }
  
  if(all(names(signal) == c('Alignment_Name','Alignment_Position','Site_Pattern','Gap','Non_Base_Taxa','Non_Base_Count','Singleton_Taxa','Split_1','Split_2','Split_3','Split_4','Split_5'))){
    
    signal_taxa <- signal %>%
      filter(!is.na(Split_1)) %>%
      head(1) %>%
      select(Non_Base_Taxa,starts_with('Split_')) %>%
      select_if(~ !any(is.na(.))) %>%
      unite(col = "Taxa",sep = ";") %>%
      pull() %>% as.character() %>% str_split(pattern = ";") %>% unlist() %>% sort()
    
    if(all(clade %in% signal_taxa)){
      
      informative_patterns <- c('non_base_biallelic','non_base_gap_biallelic',
                                'non_base_triallelic','non_base_gap_triallelic',
                                'non_base_quadallelic','non_base_gap_quadallelic',
                                'non_base_pentallelic','non_base_gap_pentallelic',
                                'biallelic','gap_biallelic',
                                'triallelic','gap_triallelic',
                                'quadallelic','gap_quadallelic',
                                'pentallelic','gap_pentallelic')
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
      
      if(only_gap){
        if(!include_gap){
          stop("Cannot only use (only_gap) and exclude (include_gap) gap positions.")
        }
        signal <- signal %>%
          filter(Gap==TRUE)
      }
      
      signal_table <- signal %>% 
        select(starts_with('Split_')) %>% 
        unlist() %>% 
        table()
      
      clade_semi <- sort(clade) %>% paste(collapse = ";")
      
      return(as.integer(Rboretum::tableCount(signal_table,clade_semi)))
      
      } else{ stop("Some taxa from 'clade' not present in alignment signal data.") }
  } else{ stop("'signal' argument must be the direct output from alignment.signal()") }
}