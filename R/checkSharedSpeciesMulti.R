#' Check for 3+ Shared Species Among Trees in Multiphylo
#'
#' This function takes a multiphylo object and returns TRUE if 3 or more species are shared among all trees; otherwise, FALSE
#' @param trees Multiphylo object
#' @return TRUE (all trees share >= 3 species) or FALSE
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#' checkSharedSpeciesMulti(trees)
#'
checkSharedSpeciesMulti <- function(trees){

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

  if(tree_count %in% species_check_df$Freq){
    shared_species <- species_check_df %>% filter(Freq == tree_count) %>% pull(species_check) %>% unlist() %>% as.character()
  }
  else{
    return(FALSE)
  }

  if(length(shared_species)>=3){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
