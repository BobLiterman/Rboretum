#' Compare Splits for Two Unique Topologies
#'
#' This function returns a dataframe with all splits in tree_1 and tree_2, with common and unique splits noted in the 'Split_Group' column [Shared,Unique_Tree_1,Unique_Tree_2]
#' @param tree_1 Phylo object
#' @param tree_2 Phylo object
#' @return Dataframe with all splits in tree_1 and tree_2, with common and unique splits noted in the 'Split_Group' column [Shared,Unique_Tree_1,Unique_Tree_2]
#' @export
#' @examples
#' compareSplits(tree_1,tree_2)
#'

compareSplits <- function(tree_1,tree_2){

  # Ensure trees are comparable
  if(!Rboretum::checkComparable(tree_1,tree_2)){
    stop("Trees not comparable (Trees either share <3 species, have identical topology after pruning, or cannot be processed.)")
  }

  # Prune trees to shared species
  shared_species <- Rboretum::getSharedSpecies(tree_1,tree_2)
  tree_1 <- Rboretum::getTrimmedTree(tree_1,shared_species)
  tree_2 <- Rboretum::getTrimmedTree(tree_2,shared_species)

  tree_1_splits <- Rboretum::getAllSplits(getTrimmedTree(tree_1,shared_species))
  tree_1_clade <- tree_1_splits %>% pull(Clade) %>% as.character()
  tree_1_mirror <- tree_1_splits %>% pull(Mirror_Clade) %>% as.character()
  tree_1_allSplits <- sort(c(tree_1_clade,tree_1_mirror))

  tree_2_splits <- Rboretum::getAllSplits(getTrimmedTree(tree_2,shared_species))
  tree_2_clade <- tree_2_splits %>% pull(Clade) %>% as.character()
  tree_2_mirror <- tree_2_splits %>% pull(Mirror_Clade) %>% as.character()
  tree_2_allSplits <- sort(c(tree_2_clade,tree_2_mirror))


  shared_splits <- tree_1_splits %>% filter(Clade %in% tree_2_allSplits) %>% select(-Split_Node) %>% `row.names<-`(NULL) %>% mutate(Split_Group = 'Shared')
  tree_1_splits <- tree_1_splits %>% filter(!(Clade %in% tree_2_allSplits)) %>% select(-Split_Node) %>% `row.names<-`(NULL) %>% mutate(Split_Group = 'Unique_Tree_1')
  tree_2_splits <- tree_2_splits %>% filter(!(Clade %in% tree_1_allSplits)) %>% select(-Split_Node) %>% `row.names<-`(NULL) %>% mutate(Split_Group = 'Unique_Tree_2')

  return(rbind(shared_splits,tree_1_splits,tree_2_splits))
}
