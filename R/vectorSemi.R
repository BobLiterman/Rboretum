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
  
  if(missing(char_vec)){
    stop("vectorSemi expects a character argument")
  } else if(is.na(char_vec)){
    return(NA)
  } else if(!is.character(char_vec)){
    stop("vectorSemi expects a character argument")
  }
  
  # If vector has only 1 element, return element
  if(length(char_vec)==1){
    return(char_vec)
  } else{
    return(paste(char_vec,collapse = ";"))
  }
}
