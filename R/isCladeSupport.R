#' Rboretum Clade Support Checker
#'
#' This function returns TRUE if all the clades in 'tree' are present in the clade comparsion output from getTreeClades(return_counts=TRUE); otherwise, FALSE
#' @param test_object R object to check
#' @param tree Check that passed support object contains information about a specific set of trees. Options include:
#' \itemize{
#'   \item A rooted phylo object; or,
#'   \item A named, rooted multiPhylo object where all trees share 3+ taxa
#' }
#' @return TRUE if output of getTreeClades(tree,return_count=TRUE); otherwise, FALSE
#' @export

isCladeSupport <- function(test_object,tree){
  
  # Check clade support data frame
  if(!is.data.frame(test_object)){
    return(FALSE) # Clade support is a data frame
  } else if(ncol(test_object)!=3){
    return(FALSE) # Clade support has 3 columns
  } else if(!all(names(test_object)==c('Clade','Count','Trees'))){
    return(FALSE) # Check column names
  } else{
    support_clades <- test_object %>% pull(Clade) %>% as.character()
  }
  
  # Get tree clades
  tree_clades <- getTreeClades(tree)
 
  # If all tree clades are in the clade comparison, return TRUE; otherwise FALSE
  if(has_error(silent=TRUE,expr=!all(tree_clades %in% support_clades))){
    return(FALSE)
  } else if(!all(tree_clades %in% support_clades)){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
  