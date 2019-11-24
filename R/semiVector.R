#' Rboretum Semicolon Vectorizer
#'
#' This function returns a vector of items from a string separated by semicolons (;)
#' @usage semiVector(string_to_split)
#' @param string_to_split Semicolon delimited string
#' @return Character vector
#' @export
#' @examples
#' myString <- 'a;b;c;d'
#' semiVector(myString)
#' > ['a','b','c','d']
#'

# Given semicolon delimeted string, return single vector
semiVector <- function(string_to_split){
  if(!str_detect(string_to_split,";")){
    return(string_to_split) # Return strings with no semicolon
  } else{
  return(stringr::str_split(string_to_split,pattern = ';') %>% unlist())
  }
}
