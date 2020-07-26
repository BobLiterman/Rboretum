#' Rboretum Semi Checker
#'
#' This function returns TRUE if the passed object is a character vector containing semicolon-separated taxa. If return_count = TRUE, it will return the number of elements in each element after splitting
#' @param test_object Object to check
#' @param return_count OPTIONAL: If TRUE and test_object is a semicolon-separated character, return the number of elements after splitting [Default: FALSE]
#' @return TRUE if test_object is a semicolon-separated character (or element count if return_count = TRUE), or FALSE
#' @export

semiChecker <- function(test_object,return_count){
  
  # Ensure test_object is present and a character
  if(missing(test_object)){
    stop("'test_object' not provided...")
  } else if(!is.character(test_object)){
    return(FALSE)
  }
  
  # Return logical or count?
  if(missing(return_count)){
    return_count <- FALSE
  } else if(!is.logical(return_count)){
    return_count <- FALSE
  } else if(length(return_count)!=1){
    return_count <- FALSE    
  }
  
  # Process single values
  if(length(test_object)==1){
    
    if(str_detect(test_object,";")){ # If semicolons are detected...
      
      # Return count of items or TRUE
      if(return_count){
        return(length(semiVector(test_object)))
        } else{
          return(TRUE)
        }
      } else{
        return(FALSE)
      }
  } else{ # If character vector contains multiple elements, ensure all are semicolon-separated
      semi_check <- purrr::map(.x=test_object,.f=function(x){str_detect(x,";")}) %>% unlist() %>% all()
      
      if(semi_check){
        clade_counts <- purrr::map(.x=test_object,.f=function(x){length(semiVector(x))}) %>% unlist()
        
        # Return count of items or TRUE
        if(return_count){
          return(clade_counts)
        } else{
          return(TRUE)
        }
      } else{
        return(FALSE)
      }
    }
  }