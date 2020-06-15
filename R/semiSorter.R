#' Rboretum Semicolon Sorter
#'
#' This function returns a sorted, semicolon-delimited string from a character vector or an unsorted semicolon-delimited string
#' @param string_to_split Character vector or unsorted semicolon delimited string
#' @return Sorted semicolon delimited string
#' @examples 
#' myShuffledString <- 'c;d;a;b'
#' semiSorter(myShuffledString)
#' > 'a;b;c;d'
#' @export

semiSorter <- function(string_to_sort){
  
  if(!is.character(string_to_sort)){
    stop("semiSorter expects a character argument")
  }
  if(length(string_to_sort)==1){
    if(suppressWarnings(str_detect(string_to_sort,";"))){ # If already semicolon-separated, return sorted character vector
      return(paste(naturalsort(semiVector(string_to_sort)),collapse = ";"))
    } else  # Return single elements
      return(string_to_sort)
  } else{ # Return sorted + joined character
    return(paste(naturalsort(string_to_sort),collapse = ";"))
  }
}
