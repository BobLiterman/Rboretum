#' Check Monophyly of Taxa
#'
#' This function returns TRUE if  'taxa' are monophyletic in 'tree', FALSE otherwise
#' @param tree Phylo object
#' @param taxa Vector containing taxa to check
#' @return TRUE if 'taxa' monophyletic in 'tree', else, FALSE
#' @export
#' @examples
#' checkMonophyletic(tree,taxa_to_check)
#'

# Returns TRUE if 'taxa' are monophyletic in 'tree'
checkMonophyletic <- function(tree,taxa){
  if(!(all(taxa %in% tree$tip.label))){
    stop("Taxa missing from tree.")
  }
  if(ape::is.monophyletic(tree,taxa)){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
