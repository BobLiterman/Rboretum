#' Rboretum Vector Semicolonizer
#'
#' This function takes a character vector and retuns a string separated by semicolons (;)
#' @param char_vec Character vector
#' @return Semicolon-separated character vector of length 1
#' @examples
#' myVec <- c('a','b','c','d')
#' vectorSemi(myVec)
#' > 'a;b;c;d'
#' @export

vectorSemi <- function(char_vec){
  return(paste(char_vec,collapse = ";"))
}
