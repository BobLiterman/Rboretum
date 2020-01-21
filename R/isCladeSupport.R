#' Rboretum Clade Support Checker
#'
#' This function returns TRUE if all the clades in 'tree' are present in the clade comparsion output from getTreeClades(return_counts=TRUE); otherwise, FALSE
#' @param test_object R object to check
#' @param tree Check that passed support object contains information about a specific set of trees. Options include:
#' \itemize{
#'   \item A rooted phylo object; or,
#'   \item A named, rooted multiPhylo object where all trees share 3+ taxa
#' }
#' @param partial OPTIONAL: If TRUE, as long as all clades from 'tree' are present, return TRUE [Default: FALSE, require that 'tree' also includes all clades in support object]
#' @return TRUE if output of getTreeClades(tree,return_count=TRUE); otherwise, FALSE
#' @export

isCladeSupport <- function(test_object,tree,partial){
  
  # Check clade support data frame
  if(!is.data.frame(test_object)){
    return(FALSE) # Clade support is a data frame
  } else if(ncol(test_object)<3){
    return(FALSE) # Clade support has 3 columns
  } else if(!all(names(test_object)[1:3]==c('Clade','Count','Trees'))){
    return(FALSE) # Check column names
  } else{
    support_clades <- test_object %>% pull(Clade) %>% as.character()
    support_clades <- purrr::map(.x=support_clades,.f = function(x){Rboretum::semiSorter(x)}) %>% unlist() %>% as.character() %>% sort() # Sort clades
  }
  
  # Get tree clades
  tree_clades <- getTreeClades(tree) %>% as.character()
  tree_clades <- purrr::map(.x=tree_clades,.f = function(x){Rboretum::semiSorter(x)}) %>% unlist() %>% as.character() %>% sort() # Sort clades
  
 
  # Allow partial matches?
  if(missing(partial)){
    partial <- FALSE
  } else if(!is.logical(partial)){
    partial <- FALSE
  }
  
  if(partial){ # Partial matches are true if all clades from tree are present in support
    
    if(all(tree_clades %in% support_clades)){
      return(TRUE) 
    } else{
      return(FALSE) 
    }
    
  } else{ 
    if(length(tree_clades) != length(support_clades)){
      return(FALSE) 
    } else{
      if(all(tree_clades == support_clades)){
        return(TRUE)
      } else{
        return(FALSE)
      }
    }
  }
}
  