#' Check for Identical Species Among Trees in Multiphylo
#'
#' This function takes a multiphylo object and returns TRUE if all species are shared among all trees; otherwise, FALSE
#' @param trees Multiphylo object
#' @return TRUE if all trees have all species in common; else, FALSE
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3,...)
#' checkIdenticalSpeciesMulti(trees)
#'
checkIdenticalSpeciesMulti <- function(trees){

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

  if(any(species_check_df$Freq < tree_count)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}
