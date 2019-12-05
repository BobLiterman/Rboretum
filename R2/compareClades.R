#' Rboretum Tree Clade Comparison
#'
#' This function takes a named multiPhylo object, and queries support for all uniquely identified monophyletic clades
#' @param trees Named multiPhylo object [Set tree names by using names(trees) <- c('Tree1_Name','Tree2_Name',etc.)]
#' @return Dataframe with clade support information 
#' @export

compareClades <- function(trees){
  
  if(!Rboretum::isMultiPhylo(trees)){
    stop("'trees' does not appear to be a valid multiPhylo object with 2+ trees")
  } else if(!Rboretum::checkSharedTaxa(trees)){
    stop("Trees do not share at least three common species.")
  } else if(Rboretum::checkSameTopology(trees)){
    stop("Tree topologies are identical, and thus cannot be compared.")
  } else if(is.null(names(trees))){
    stop("Trees in multiPhylo must be named for compareClades. Name trees via names(trees) <- c('Tree1_Name','Tree2_Name',etc.)")
  }
  
  tree_count <- length(trees)
  tree_names <- names(trees)

  if(!Rboretum::checkSameTaxa(trees)){
    shared_species <- Rboretum::getSharedTaxa(trees)
    trees <- Rboretum::treeTrimmer(trees,shared_species)
    names(trees) <- tree_names
  }

  # Tally splits
  all_clades <- c()

  for(i in 1:tree_count){
    all_clades <- c(all_clades,Rboretum::getTreeClades(trees[[i]]))
  }

  tallied_splits <- as.data.frame(table(all_clades))

  split_df <- tallied_splits %>%
    rename(Clade = 'all_clades',Tree_Count = 'Freq') %>%
    mutate(Clade_Size = (str_count(Clade,';')+1)) %>%
    mutate(Clade = as.character(Clade),Tree_Count = as.integer(Tree_Count),Clade_Size = as.integer(Clade_Size)) %>%
    mutate(Tree_Percent = (as.numeric(Tree_Count)/as.numeric(tree_count))*100)
  
  tree_column <- c()

  for(clade in split_df$Clade){
    tree_list <- c()
    for(i in 1:tree_count){
      if(ape::is.monophyletic(trees[[i]],Rboretum::semiVector(as.character(clade)))){
        tree_list <- c(tree_list,tree_names[i])
      }
    }
    tree_column <- c(tree_column,paste(tree_list,collapse = ";"))
  }
  split_df$Trees_with_Clade <- tree_column
  return(split_df)
}
