#' Get Shared Species List from Multiphylo Object
#'
#' This function takes a multiphylo object and returns a sorted vector of taxa shared among all trees
#' @param trees Multiphylo object
#' @return Sorted vector of taxa shared among all trees in multiphylo
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#'
#' getSharedSpeciesMulti(trees)
#'

getSharedSpeciesMulti <- function(trees){

  # Ensure multiphylo has >= 2 trees
  tree_count <- length(trees)

  if(!tree_count>=2){
    stop("At least two trees are required for comparison.")
  }

  # Check that all trees contain all species
  species_check <- c()
  for(i in 1:tree_count){
    temp_tree <- trees[[i]]
    species_check <- c(species_check,temp_tree$tip.label)
  }

  species_check_df <- as.data.frame(table(species_check))
  shared_species <- species_check_df %>% filter(Freq == tree_count) %>% pull(species_check) %>% unlist() %>% as.character()

  if(length(shared_species)>=3){
    return(sort(shared_species))
  } else{
    stop("Trees share fewer than 3 species.")
  }
}
