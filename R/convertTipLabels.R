#' Convert Tip Labels
#'
#' This function takes a tree and changes tip labels based on a supplied dataframe
#' @param tree Phylo object
#' @param name_df Dataframe with name equivalencies
#' @param from Column name with current tree tip IDs
#' @param to Column name with desired IDs
#' @return Tree objecct with converted names
#' @export
#' @examples
#' convertTipLabels(tree,name_df,from,to)
#'

convertTipLabels <- function(tree,name_df,from,to){
  tree_taxa <- tree$tip.label

  to <- as.character(to)
  from <- as.character(from)

  current_ids <- as.character(pull(name_df,from))
  new_ids <- as.character(pull(name_df,to))

  if(!all(tree_taxa %in% current_ids)){
    stop("Specified 'from' column does not contain all IDs from tree")
  }

  if(any(duplicated(current_ids))){
    stop("Specified 'from' column contains duplicate IDs")
  }

  new_id_list <- c()
  for(old_id in tree_taxa){
    new_id_list <- c(new_id_list,new_ids[match(old_id,current_ids)])
  }

  tree$tip.label <- new_id_list
  return(tree)
}
