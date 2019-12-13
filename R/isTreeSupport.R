#' Rboretum Tree Support Checker
#'
#' This function returns TRUE if the passed object is the output of getTreeSupport; otherwise, FALSE
#' @param test_object R object to check
#' @param tree Check that passed support object contains information about a specific set of trees. Options include:
#' \itemize{
#'   \item A rooted phylo object; or,
#'   \item A named, rooted multiPhylo object where all trees share 3+ taxa
#' }
#' @param partial OPTIONAL: Allow for parital matches (as long as all clades from 'tree' are present, return TRUE) [Default: FALSE, test all clades]
#' @return TRUE if output of getTreeSupport(tree,...); otherwise, FALSE
#' @export

isTreeSupport <- function(test_object,tree,partial){
  
  # Ensure dataframe columns match expected and that data exists
  if(!is.data.frame(test_object)){
    return(FALSE) # Tree support is a df
  } else if(names(test_object)[1] != "Clade"){
    return(FALSE) # Tree support starts with a clade column
  } else if(!Rboretum::isPhylo(tree,check_rooted = TRUE) & !Rboretum::isMultiPhylo(tree,check_rooted = TRUE,check_named = TRUE,check_three_taxa = TRUE)){
    stop("'tree' is not a valid rooted phylo or named, rooted multiPhylo object")
  }
  
  # Allow partial matches?
  if(missing(partial)){
    partial <- FALSE
  } else if(!is.logical(partial)){
    partial <- FALSE
  }
  
  # Get tree clades and clades from getTreeSupport output
  support_clades <- sort(as.character(test_object$Clade))
  tree_clades <- sort(as.character(Rboretum::getTreeClades(tree)))
  
  if(partial){ # Partial matches are true if all clades from tree are present in support
    
    if(all(tree_clades %in% support_clades)){
      return(TRUE) 
    } else{
      return(FALSE) 
    }
    
  } else{ # Complete matches are true if 'tree' and support contain identical clade information
    
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