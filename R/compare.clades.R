#' Rboretum Tree Clade Comparison
#'
#' This function takes a multiPhylo object, and queries support for all uniquely identified monophyletic clades
#' @param trees multiPhylo object (If named, names will be used in analysis)
#' @param return_shared_only OPTIONAL: If TRUE, returns only information about clades supported by all trees; (Default: FALSE)
#' @return Dataframe with clade support information 
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#'
#' # Return information on all possible splits
#' compare.clades(trees)
#'
#' # Return information on only those splits shared in all trees
#' compare.clades(trees,return_shared_only=TRUE)
#'

compare.clades <- function(trees,return_shared_only){
  
  if(!Rboretum::is.multiPhylo(trees)){
    stop("'trees' does not appear to be a valid multiPhylo object with 2+ trees")
  } else if(!Rboretum::check.shared(trees)){
    stop("Trees do not share at least three common species.")
  } else if(Rboretum::same.topology(trees)){
    stop("Tree topologies are identical, and thus cannot be compared.")
  } 
  
  tree_count <- length(trees)
  
  if(is.null(names(trees))){
    tree_names <- unlist(lapply(X = 1:tree_count,function(x) paste(c("Tree",x),collapse = "_")))
  } else{
    tree_names <- names(trees)
  }

  if(missing(return_shared_only)){
    return_shared_only <- FALSE
  } else if(!is.logical(return_shared_only)){
    return_shared_only <- FALSE
    }
  
  if(!Rboretum::same.taxa(trees)){
    shared_speces <- Rboretum::get.shared(trees)
    trees <- Rboretum::trim.tree(trees,shared_speces)
  }

  # Tally splits
  all_clades <- c()

  for(i in 1:tree_count){
    all_clades <- c(all_clades,Rboretum::get.clades(trees[[i]]))
  }

  tallied_splits <- as.data.frame(table(all_clades))

  split_df <- tallied_splits %>%
    rename(Clade = 'all_clades',Tree_Count = 'Freq') %>%
    mutate(Clade_Size = (str_count(Clade,';')+1)) %>%
    mutate(Clade = as.character(Clade),Tree_Count = as.integer(Tree_Count),Clade_Size = as.integer(Clade_Size)) %>%
    filter(Clade_Size > 1)

  # If return_shared_only = TRUE, just return clades shared by all trees. Otherwise, return all clades
  if(return_shared_only){
    split_df <- split_df %>% 
      filter(Tree_Count == tree_count) %>%
      select(Clade,Clade_Size)
    return(split_df)
  } else{
    
    split_df <- split_df %>%
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
}
