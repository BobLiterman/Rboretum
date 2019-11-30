#' Rboretum Alignment Signal Support Mapper
#'
#' This function computes the alignment support for a specified tree, processing signal from all alignment columns with <= 'max_missing' missing taxa [from output of getAlignmentSignal()]. Signal from multiple alignments/missing combinations can be added to the same table by passing the output of tree.support() to the 'existing_splits' argument
#' @param signal Output table from getAlignmentSignal()
#' @param tree phylo object where taxa EXACTLY match those from signal
#' @param max_missing OPTIONAL: Number of missing sites allowed in alignment column [Default: 0]
#' @param alignment_name OPTIONAL: Column name for data being added [Default: Alignment name from signal dataframe]
#' @param include_gap OPTIONAL: TRUE or FALSE; Count sites with gap positions ('-') as part of total support [Default: TRUE]
#' @param include_singleton OPTIONAL: TRUE or FALSE; Count sites with singletons as part of total support [Default: TRUE]
#' @param include_biallelic OPTIONAL: TRUE or FALSE; Count sites with biiallelic variation as part of total support [Default: TRUE]
#' @param include_triallelic OPTIONAL: TRUE or FALSE; Count sites with triallelic variation as part of total support [Default: TRUE]
#' @param include_quadallelic OPTIONAL: TRUE or FALSE; Count sites with quadallelic variation as part of total support [Default: TRUE]
#' @param include_pentallelic OPTIONAL: TRUE or FALSE; Count sites with pentallelic variation as part of total support [Default: TRUE]
#' @param only_gap OPTIONAL: TRUE or FALSE; Only count sites with gap positions ('-') as part of total support [Default: FALSE]
#' @param existing_splits OPTIONAL: Output from previous tree.support() using the same tree, run with new alignment/missing combination
#' @return The same split table from getTreeSplits(tree), but with a support column for the specfied alignment/missing combination
#' @export

getTreeSupport <- function(signal,tree,max_missing,alignment_name,include_gap,include_singleton,include_biallelic,include_triallelic,include_quadallelic,include_pentallelic,only_gap,existing_splits){

  if(!all(names(signal) == c('Alignment_Name','Alignment_Position','Site_Pattern','Gap','Singleton','Singleton_Taxa','Non_Base_Taxa','Non_Base_Count','Split_1','Split_2','Split_3','Split_4','Split_5'))){
    stop("'signal' argument must be output from getAlignmentSignal()")
  }
  
  if(has_error(silent=TRUE,expr=ape::is.rooted(tree))){
    stop("Error in ape::is.rooted. Is 'tree' a phylo object?")
  } else if(!ape::is.rooted(tree)){
    stop("Tree must be rooted for tree.support")}
  
  if(missing(max_missing)){
    max_missing <- 0
  }
  
  if(missing(alignment_name)){
    alignment_name <-  paste(c(as.character(signal$Alignment_Name[1]),'_m',as.character(max_missing)),collapse = '')
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
  
  if(missing(existing_splits)){
    existing_splits <- FALSE
  } else if(is.data.frame(existing_splits) & all(names(existing_splits)[1:4] == c('Clade','Mirror_Clade','Split_Node','Split_Bootstrap'))){
    
    old_splits <- existing_splits
    old_clades <- old_splits %>% pull(Clade) %>% as.character() %>% sort()
    old_mirror_clades <- old_splits %>% pull(Mirror_Clade) %>% as.character() %>% sort()
    existing_splits <- TRUE
  
    } else{
    print("Argument passed to 'existing_splits' should be output from tree.support(). Returning results for tree and alignment/missing combo specified.")
    existing_splits <- FALSE
  }
  
  signal_taxa <- signal %>%
    filter(!is.na(Split_1)) %>% 
    head(1) %>% 
    select(Singleton_Taxa,Non_Base_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5) %>% 
    select_if(~ !any(is.na(.))) %>% 
    unite(col = "Taxa",sep = ";") %>% 
    semiVector() %>% 
    sort()

  tree_taxa <- sort(tree$tip.label)

  if(!all(tree_taxa == signal_taxa)){
    print("Tree Taxa:")
    print(tree_taxa)
    print("Signal Taxa:")
    print(signal_taxa)
    stop("Taxa from signal analysis doesn't match that from tree.")
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
    clade_support <- c(clade_support,tableCount(all_signal_splits,clade))
  }
  
  support_df <- data.frame(Clade = clades,Support = clade_support) %>%
    mutate(Clade = as.character(Clade),Support = as.integer(Support)) %>% 
    left_join(splits,by='Clade') %>%
    select(Clade,Mirror_Clade,Split_Node,Split_Bootstrap,Support) %>%
    rename(!!alignment_name := Support)

  if(existing_splits){
    if(all(clades == old_clades) & all(mirror_clades == old_mirror_clades)){
      if(alignment_name %in% names(old_splits)){
        print("Column exists in 'existing_splits' with that alignment/missing combination. Returning results for tree and alignment/missing combo specified.")
        return(support_df)
      } else{
        support_df <- left_join(old_splits,support_df,by=c("Clade", "Mirror_Clade", "Split_Node", "Split_Bootstrap"))
        print(paste(c("Added tree split data from",alignment_name),collapse = " "))
        return(support_df)
      }
    } else{
      print("Clades don't match between pre-existing splits and tree provided. Returning results for tree and alignment provided.")
      return(support_df)}
    } else{ return(support_df)}
}