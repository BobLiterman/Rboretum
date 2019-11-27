#' Rboretum Concordance Counter
#'
#' This function takes the signal from an alignment file [via getAlignmentSignal()] and counts the number of splits that provide concordant signal, relative to a given tree
#' @param signal Output table from getAlignmentSignal()
#' @param tree Rooted phylo object
#' @return Dataframe containing the subset of 'signal' that has split information, and the count of concordant splits for each site
#' @export
#'

concordance.counter <- function(signal,tree){
  
  if(!all(names(signal)==c('Alignment_Name','Alignment_Position','Site_Pattern','Gap','Singleton','Singleton_Taxa','Non_Base_Taxa','Non_Base_Count','Split_1','Split_2','Split_3','Split_4','Split_5'))){
    stop("'signal' must be output from getAlignmentSignal()")
  } else{
    signal_taxa <- signal %>%
      filter(!is.na(Split_1)) %>% 
      head(1) %>% 
      select(Singleton_Taxa,Non_Base_Taxa,Split_1,Split_2,Split_3,Split_4,Split_5) %>% 
      select_if(~ !any(is.na(.))) %>% 
      unite(col = "Taxa",sep = ";") %>% 
      semiVector() %>% 
      sort()
  }
  
  if(!Rboretum::isPhylo(tree)){
    stop("'tree' does not appear to be a valid phylo object")
  } else if(!ape::is.rooted(tree)){
    stop("'tree' must be rooted for concordance.counter()")
  } else{
    tree_taxa <- sort(tree$tip.label)
  }
  
  if(!all(tree_taxa == signal_taxa)){
    stop("Taxa from signal analysis doesn't match that from tree.")
  }
  
  signal$Signal_Taxa <- signal %>%
    select(starts_with('Split')) %>%
    unite(col='Taxa',sep=";",na.rm=TRUE) %>% 
    pull()
  
  signal <- signal %>%
    rowwise() %>%
    mutate(Signal_Taxa = semiSorter(Signal_Taxa)) %>%
    ungroup()
  
  unique_taxa_sets <- unique(signal$Signal_Taxa)
  
  if(length(unique_taxa_sets) == 1){
    tree <- Rboretum::treeTrimmer(tree,semiVector(unique_taxa_sets[1]))
    oneTree <- TRUE
  } else{
    trees <- c(rtree(10),rtree(10))
    for(i in 1:length(unique_taxa_sets)){
      trees[[i]] <-Rboretum::treeTrimmer(tree,semiVector(unique_taxa_sets[i]))
    }
    names(trees) <- unique_taxa_sets
    oneTree <- FALSE
  }
  
  signal$Conc_Sites <- 0
  biallelic_signal <- signal %>% filter(Site_Pattern == 'biallelic')
  triallelic_signal <- signal %>% filter(Site_Pattern == 'triallelic')
  quadallelic_signal <- signal %>% filter(Site_Pattern == 'quadallelic')
  pentallelic_signal <- signal %>% filter(Site_Pattern == 'pentallelic')
  
  if(oneTree){
    print("Calculating concordance counts for a single set taxa...")
    if(nrow(biallelic_signal) > 0){
      biallelic_signal <- biallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=tree,tips=semiVector(Split_1)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_2))))) %>%
        ungroup()
      print("Finished biallelic sites...")
    }
    
    if(nrow(triallelic_signal) > 0){
      triallelic_signal <- triallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=tree,tips=semiVector(Split_1)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_2)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_3))))) %>%
        ungroup()
      print("Finished triallelic sites...")
    }
    
    if(nrow(quadallelic_signal) > 0){
      quadallelic_signal <- quadallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=tree,tips=semiVector(Split_1)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_2)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_3)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_4))))) %>%
        ungroup()
      print("Finished quadallelic sites...")
    }
    
    if(nrow(pentallelic_signal) > 0){
      pentallelic_signal <- pentallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=tree,tips=semiVector(Split_1)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_2)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_3)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_4)),ape::is.monophyletic(phy=tree,tips=semiVector(Split_5))))) %>%
        ungroup()
      print("Finished pentallelic sites...")
    }
    
    all_signal <- rbind(biallelic_signal,triallelic_signal,quadallelic_signal,pentallelic_signal)
  } else{
    print("Calculating concordance counts for different sets of available taxa...")
    
    if(nrow(biallelic_signal) > 0){
      biallelic_signal <- biallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_1)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_2))))) %>%
        ungroup()
      print("Finished biallelic sites...")
    }
    
    if(nrow(triallelic_signal) > 0){
      triallelic_signal <- triallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_1)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_2)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_3))))) %>%
        ungroup()
      print("Finished triallelic sites...")
    }
    
    if(nrow(quadallelic_signal) > 0){
      quadallelic_signal <- quadallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_1)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_2)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_3)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_4))))) %>%
        ungroup()
      print("Finished quadallelic sites...")
    }
    
    if(nrow(pentallelic_signal) > 0){
      pentallelic_signal <- pentallelic_signal %>%
        rowwise() %>%
        mutate(Conc_Sites = sum(c(ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_1)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_2)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_3)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_4)),ape::is.monophyletic(phy=trees[[Signal_Taxa]],tips=semiVector(Split_5))))) %>%
        ungroup()
      print("Finished pentallelic sites...")
    }
    
    all_signal <- rbind(biallelic_signal,triallelic_signal,quadallelic_signal,pentallelic_signal)
  }
  
  return(all_signal)
}