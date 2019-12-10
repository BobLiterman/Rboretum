#' Rboretum Alignment Signal Support Mapper
#'
#' This function maps phylogenetic signal from a multiple-sequence alingment onto a rooted phylo or multiPhylo object
#' @param tree A rooted phylo object onto which alignment signal will be mapped.
#' @param signal Output table from getAlignmentSignal() run with taxa against 'tree'
#' @param alignment_name Column name(s) for data being added [Default: Alignment name from signal dataframe + 'm_<MISSING>']
#' @return The same split table from getTreeSplits(tree), but with (1) no root row, and (2) a support column for the specfied alignment/missing combination
#' @export

getTreeSupport_Worker <- function(tree,signal,alignment_name){

  # Check if tree and signal are valid and compatible
  if(!Rboretum::isPhylo(tree,check_rooted = TRUE)){
    stop("'tree' must be a rooted phylo object for getTreeSupport_Worker")
  } else if(!Rboretum::isAlignmentSignal(signal,tree)){
    stop("'signal' is either not the output from getAlignmentSignal(), or does not contain the same taxa as 'tree'")
  } else {
    raw_alignment_name <- unique(signal$Alignment_Name)
    alignment_count <- length(raw_alignement_name)
    if(alignment_count != length(alignment_name)){
      stop("'alignment_name' and number of detected alignments do not match")
    }
  } 
  
  support_df <- Rboretum::getTreeSplits(tree) %>%
    filter(!is.na(Split_Node)) %>% 
    mutate(Clade = as.character(Clade),Mirror_Clade = as.character(Mirror_Clade))
    
  clades <- support_df %>% pull(Clade) %>% as.character() %>% sort()

  for(i in 1:alignment_count){
    
    filter_name <- raw_alignment_name[i]
    alignment_col <- alignment_name[i]
    
    all_signal_splits <- signal %>% filter(Alignment_Name == filter_name) %>% select(starts_with('Split_')) %>% unlist() %>% table()
    
    clade_support <- c()
    for(clade in clades){
      clade_support <- c(clade_support,tableCount(all_signal_splits,clade))
    }
    
    temp_df <- data.frame(Clade = clades,Support = clade_support) %>%
      mutate(Clade = as.character(Clade),Support = as.integer(Support)) %>%
      rename(!!alignment_col := Support)
    
    support_df <- support_df %>%
      left_join(temp_df,by='Clade')
  }
  return(support_df)
}
