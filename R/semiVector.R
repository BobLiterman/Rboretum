#' Rboretum Semicolon Vectorizer
#'
#' This function takes a string separated by semicolons (;) and returns a character vector
#' @param string_to_split Semicolon delimited character vector of length 1
#' @return Character vector
#' @examples
#' myString <- 'a;b;c;d'
#' semiVector(myString)
#' > ['a','b','c','d']
#' @export

semiVector <- function(string_to_split){
  
  # Ensure character vector of length 1
  if(!is.character(string_to_split)){
    stop("semiVector requires a character argument")
  } else if(length(string_to_split)!=1){
    stop("semiVector requires a character argument with a single element")
  }
  
  # If no ';', return string
  if(!stringr::str_detect(string_to_split[1],";")){
    return(string_to_split)
  } else{ # Return vector split by ;
    return(stringr::str_split(string_to_split[1],pattern = ';') %>% unlist())
  }
}