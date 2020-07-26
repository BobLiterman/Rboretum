#' Rboretum Vector Semicolonizer
#'
#' This function takes a character vector, or a list of vectors, and retuns a string separated by semicolons (;)
#' @param char_vec Character vector or list of character vectors
#' @return Semicolon-separated strings for each passed vector
#' @examples
#' myVec <- c('a','b','c','d')
#' vectorSemi(myVec)
#' > 'a;b;c;d'
#' @export

vectorSemi <- function(to_semi){
  
  # Check argument
  if(missing(to_semi)){
    stop("vectorSemi expects a character or list argument")
  } 
  
  # If NA is passed, return NA
  if(length(string_to_split)==1){
    if(is.na(string_to_split)){
      return(NA)
    }
  }
  
  if(!is.list(to_semi) & !(is.character(to_semi))){
    stop("vectorSemi expects a character or list argument")
  }
  
  # Process character argument
  if(is.character(to_semi)){
    return(paste(to_semi,collapse = ";"))
  }
  
  # Process list argument
  if(is.list(to_semi)){
    semi_list <- purrr::map(.x=to_semi,.f=function(x){paste(x,collapse = ";")})
    return(semi_list)
  }
}
