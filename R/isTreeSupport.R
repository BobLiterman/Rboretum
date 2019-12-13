#' Rboretum Tree Support Checker
#'
#' This function returns TRUE if the passed object is the output of getTreeSupport() and contains information about clades in 'test_clade'; otherwise, FALSE
#' @param test_object R object to check
#' @param test_clade Character vector of semicolon-separated clades. Function checks that passed support object contains information about these clades.
#' @param partial OPTIONAL: If TRUE, as long as all clades from 'test_clade' are present, return TRUE [Default: FALSE, require that test_clades also includes all clades in support object]
#' @return TRUE if output of getTreeSupport() with information on 'test_clade'; otherwise, FALSE
#' @export

isTreeSupport <- function(test_object,test_clade,partial){
  
  # Ensure dataframe columns match expected and that data exists
  if(!is.data.frame(test_object)){
    return(FALSE) # Tree support is a df
  } else if(names(test_object)[1] != "Clade"){
    return(FALSE) # Tree support starts with a clade column
  } else{
    
    # Get sorted clades from getTreeSupport
    support_clades <- sort(as.character(test_object$Clade))
    
    if(!purrr::map(.x=support_clades,.f=function(x){str_detect(x,";")}) %>% unlist() %>% all()){
      stop("Clades in 'test_object' should be semicolon-separated character strings")
    } 
    
    support_clades <- purrr::map(.x=support_clades,.f = function(x){Rboretum::semiSorter(x)}) %>% unlist() %>% as.character() %>% sort() # Sort clades
    
  }
  
  # Ensure 'test_clade' is a character vector of semicolon-separated clades
  
  if(!is.character(test_clade)){
    stop("'test_clade' must be a character vector")
  } else if(!purrr::map(.x=test_clade,.f=function(x){str_detect(x,";")}) %>% unlist() %>% all()){
    stop("Clades in 'test_clade' should be semicolon-separated character strings")
  } else{
    test_clade <- purrr::map(.x=test_clade,.f = function(x){Rboretum::semiSorter(x)}) %>% unlist() %>% as.character() %>% sort() # Sort clades
  }
  
  # Allow partial matches?
  if(missing(partial)){
    partial <- FALSE
  } else if(!is.logical(partial)){
    partial <- FALSE
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