#' Rboretum Semicolon Sorter
#'
#' This function returns a sorted, semicolon-delimited character vector from a plain character vector, or a character vector of semicolon-delimited elements
#' @param string_to_split Character vector
#' @return Sorted semicolon delimited character vector
#' @examples 
#' myShuffledString <- 'c;d;a;b'
#' semiSorter(myShuffledString)
#' > 'a;b;c;d'
#' @export

semiSorter <- function(string_to_sort){
  
  if(missing(string_to_sort)){
    stop("semiSorter expects a character argument")
  }
  
  # If NA is passed, return NA
  if(length(string_to_sort)==1){
    if(is.na(string_to_sort)){
      return(NA)
    }
  }
  
  if(!is.character(string_to_sort)){
    stop("semiSorter expects a character argument")
  }
  
  # If string_to_sort has a single element...
  if(length(string_to_sort)==1){
    
    # If string_to_sort is already semicolon-separated, release, sort, and rejoin
    if(semiChecker(string_to_sort)){
      sorted_string <- semiVector(string_to_sort) %>% naturalsort() %>% vectorSemi()
      return(sorted_string)
    } else{
      # 'string_to_sort' is a single-element character with no ";", return...
      return(string_to_sort)
    }
  }
  
  # If string_to_sort has > 1 element....
  if(length(string_to_sort)>1){
    
    # If string_to_sort is a vector of semicolon-separated elements, sort all and return
    if(semiChecker(string_to_sort)){
      sorted_string <- purrr::map(.x=string_to_sort,.f=function(x){semiVector(x) %>% naturalsort() %>% vectorSemi()}) %>% unlist()
      return(sorted_string)
    } else{
      
      # If string_to_sort is a vector of characters, sort, join, and return
      sorted_string <- naturalsort(string_to_sort) %>% vectorSemi()
    }
  }
}