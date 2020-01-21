#' Rboretum Semicolon Vectorizer
#'
#' This function takes a string separated by semicolons (;) and returns a character vector
#' @param string_to_split Semicolon delimited string
#' @return Character vector
#' @examples
#' myString <- 'a;b;c;d'
#' semiVector(myString)
#' > ['a','b','c','d']
#' @export

semiVector <- function(string_to_split){
  return(stringr::str_split(string_to_split,pattern = ';') %>% unlist())
}
