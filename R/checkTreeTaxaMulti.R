#' Check Tree Taxa in Multiphylo
#'
#' This function returns TRUE if all 'taxa' are present in all trees in a multiphylo, FALSE otherwise
#' @param trees Multiphylo object
#' @param taxa Vector containing taxa to check
#' @return TRUE if 'taxa' in all trees, else, FALSE
#' @export
#' @examples
#' trees <- c(tree_1,tree_2,tree_3)
#' checkTreeTaxaMulti(trees,taxa)
#'

checkTreeTaxaMulti <- function(trees,taxa){
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

  if(all(taxa %in% shared_species)){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}
