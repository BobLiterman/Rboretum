#' Rboretum Tip Label Swapper
#'
#' This function takes a tree and changes tip labels based on a supplied dataframe
#' @param tree phylo object
#' @param name_df Dataframe with name equivalencies
#' @param from Column name with current tree tip IDs (if missing, first column of name_df used as default)
#' @param to Column name with desired IDs (if missing, second column of name_df used as default)
#' @return phylo object with converted names
#' @export
#' @examples
#' convert.tips(tree,name_df,from,to)
#'

convert.tips <- function(tree,name_df,from,to){
  
  if(missing(tree)){
    stop("No 'tree' provided")
  } else if(!Rboretum::is.phylo(tree)){
    stop("'tree' does not appear to be a valid phylo object.")
  } else{ tree_taxa <- tree$tip.label }
  
  if(missing(name_df)){
    stop("No name table provided to name_df")
  } else if(ncol(name_df)<2){
    stop("'name_df' must at least have two columns [current IDs, new IDs]")
  }
  
  if(missing(from)){
    print("No 'from' column name provided. Defaulting to column 1 of 'name_df'...")
    from <- names(name_df)[1]
  } else{ from <- as.character(from) }
  
  if(missing(to)){
    print("No 'to' column name provided. Defaulting to column 2 of 'name_df'...")
    to <- names(name_df)[2]
  } else{ to <- as.character(to) } 
  
  if(!all(c(to,from) %in% names(name_df))){
    stop("Either 'to' or 'from' column not found in 'name_df'")
  } else{
    current_ids <- as.character(pull(name_df,from))
    new_ids <- as.character(pull(name_df,to)) 
  }

  if(!all(tree_taxa %in% current_ids)){
    stop("Specified 'from' column does not contain all IDs from tree")
  } else if(any(duplicated(current_ids))){
    stop("Specified 'from' column contains duplicate IDs")
  }

  new_id_list <- c()
  for(old_id in tree_taxa){
    new_id_list <- c(new_id_list,new_ids[match(old_id,current_ids)])
  }

  tree$tip.label <- new_id_list
  return(tree)
}
