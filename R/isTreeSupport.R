#' Rboretum Tree Support Checker
#'
#' This function returns TRUE if the passed object is the output of getTreeSupport() and contains information about clades in 'test_clade'; otherwise, FALSE
#' @param test_object R object to check
#' @param test_clade Source for clades to check. Function checks that passed support object contains information about these clades. Options include:
#' \itemize{
#'   \item Rooted phylo or multiPhylo object
#'   \item vector of semicolon-separated clades. 
#' }
#' @param partial OPTIONAL: If TRUE, as long as all clades from 'test_clade' are present, return TRUE [Default: FALSE, require that test_clades also includes all clades in support object]
#' @return TRUE if output of getTreeSupport() with information on 'test_clade'; otherwise, FALSE
#' @export

isTreeSupport <- function(test_object,test_clade,partial){
  
  if(missing(test_object)){
    stop("'test_object' not provided")
  }
  
  if(missing(test_clade)){
    stop("'test_clade' required")
  }
  
  # Allow partial matches?
  if(missing(partial)){
    partial <- FALSE
  } else if(!is.logical(partial)){
    partial <- FALSE
  }
  
  # Ensure dataframe columns match expected and that data exists
  if(!is.data.frame(test_object)){
    return(FALSE) # Tree support is a df
  } else if(names(test_object)[1] != "Clade"){
    return(FALSE) # Tree support starts with a clade column
  }
    
  # Get sorted clades from getTreeSupport
  support_clades <- naturalsort(as.character(test_object$Clade))
  
  if(!purrr::map(.x=support_clades,.f=function(x){str_detect(x,";")}) %>% unlist() %>% all()){
    stop("Clades in 'test_object' should be semicolon-separated character strings")
  } 
  
  support_clades <- purrr::map(.x=support_clades,.f = function(x){Rboretum::semiSorter(x)}) %>% unlist() %>% as.character() %>% naturalsort() # Sort clades
  
  # Get test clades 
  if(Rboretum::isPhylo(test_clade,check_rooted = TRUE)){
    test_clade <- Rboretum::getTreeClades(test_clade,include_root = TRUE)
  } else if(Rboretum::isMultiPhylo(test_clade,check_rooted = TRUE,check_three_taxa = TRUE)){
    
    # Reduce to common taxa if necessary
    if(!Rboretum::isMultiPhylo(tree,check_all_taxa = TRUE)){
      tree <- treeTrimmer(tree)
    }
    test_clade <- Rboretum::getTreeClades(tree,include_root = TRUE)
  } else if(!is.character(test_clade)){ # Ensure 'test_clade' is a character vector of semicolon-separated clades
    stop("'test_clade' must be a phylo, multiPhylo, or character vector")
  } else if(!purrr::map(.x=test_clade,.f=function(x){str_detect(x,";")}) %>% unlist() %>% all()){
    stop("Clades in 'test_clade' should be semicolon-separated character strings")
  } else{
    test_clade <- purrr::map(.x=test_clade,.f = function(x){Rboretum::semiSorter(x)}) %>% unlist() %>% as.character() %>% naturalsort() # Sort clades
  }
  
  if(partial){ # Partial matches are true if all clades from tree are present in support
    
    if(all(test_clade %in% support_clades)){
      return(TRUE) 
    } else{
      return(FALSE) 
    }
    
  } else{ # Complete matches are true if 'test_clade' and support contain identical clade information
    
    if(length(test_clade) != length(support_clades)){
      return(FALSE) 
    } else{
      if(all(test_clade == support_clades)){
        return(TRUE)
      } else{
        return(FALSE)
      }
    }
  }
}