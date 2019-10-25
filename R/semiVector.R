#' Create Vector from Semicolon Delimited String
#'
#' This function returns a vector of items from a string separated by semicolons (;)
#' @param string_to_split Semicolon delimited string
#' @return Vector of items
#' @export
#' @examples
#' semiVector(string_to_split)
#'

# Given semicolon delimeted string, return single vector
semiVector <- function(string_to_split){
  return(stringr::str_split(string_to_split,pattern = ';') %>% unlist())
}
