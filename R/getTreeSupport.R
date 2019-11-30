#' Rboretum Alignment Signal Support Mapper
#'
#' This function computes the alignment support for a specified tree
#' @param tree phylo object onto which alignment signal will be mapped
#' @param signal Output table from getAlignmentSignal() run with taxa from 'tree'
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column [Default: 0]
#' @param alignment_name OPTIONAL: Column name for data being added [Default: Alignment name from signal dataframe + 'm_<MISSING>']
#' @param include_gap OPTIONAL: TRUE or FALSE; Count sites with gap positions ('-') as part of total support [Default: TRUE]
#' @param include_singleton OPTIONAL: TRUE or FALSE; Count sites with singletons as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: TRUE or FALSE; Count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: TRUE or FALSE; Count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: TRUE or FALSE; Count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: TRUE or FALSE; Count sites with pentallelic variation as part of total support [Default: TRUE]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') as part of total support [Default: FALSE]
#' @param existing_support OPTIONAL: Output from previous getTreeSupport() using the same tree, run with new alignment/missing combination
#' @return The same split table from getTreeSplits(tree), but with (1) no root row, and (2) a support column for the specfied alignment/missing combination
#' @export

getTreeSupport <- function(tree,signal,max_missing,alignment_name,include_gap,include_singleton,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic,only_gap,existing_support){

  if(has_error(silent=TRUE,expr=ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for tree.support")}
  
  if(!Rboretum::isAlignmentSignal(signal,tree)){
    stop("'signal' is either not the output from getAlignmentSignal(), or does not contain the same taxa as 'tree'")
  }
  
  if(missing(max_missing)){
    max_missing <- 0
  }
  
  if(missing(alignment_name)){
    alignment_name <- paste(c(as.character(signal$Alignment_Name[1]),'_m',as.character(max_missing)),collapse = '')
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
  
  if(!include_gap){
    if(only_gap){
      stop("Cannot only use (only_gap) and exclude (include_gap) gap positions.") }
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

  splits <- Rboretum::getTreeSplits(tree) %>%
    filter(!is.na(Split_Node))%>% 
    mutate(Clade = as.character(Clade),Mirror_Clade = as.character(Mirror_Clade))

  clades <- splits %>% pull(Clade) %>% as.character() %>% sort()
  mirror_clades <- splits %>% pull(Mirror_Clade) %>% as.character() %>% sort()
  
  all_signal_splits <- signal %>% select(starts_with('Split_')) %>% unlist() %>% table()

  clade_support <- c()
  for(clade in clades){
    clade_support <- c(clade_support,Rboretum::tableCount(all_signal_splits,clade))
  }
  
  support_df <- data.frame(Clade = clades,Support = clade_support) %>%
    mutate(Clade = as.character(Clade),Support = as.integer(Support)) %>% 
    left_join(splits,by='Clade') %>%
    select(Clade,Mirror_Clade,Split_Node,Split_Bootstrap,Support) %>%
    rename(!!alignment_name := Support)

  if(add_support){
    if(alignment_name %in% names(existing_support)){
        print("Column exists in 'existing_support' with that alignment/missing combination. Returning results for tree and alignment/missing combo specified.")
        return(support_df)
      } else{
        support_df <- left_join(existing_support,support_df,by=c("Clade", "Mirror_Clade", "Split_Node", "Split_Bootstrap"))
        print(paste(c("Added tree split data from",alignment_name),collapse = " "))
        return(support_df)
      }
  } else{ return(support_df) }
}