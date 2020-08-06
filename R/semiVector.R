#' Rboretum Semicolon Vectorizer
#'
#' This function takes a string separated by semicolons (;) and returns a character vector
#' @param string_to_split Semicolon delimited character
#' @return If string_to_split is a single item, semiVector returns a character vector split by ';'. If string_to_split contains multiple items, function returns a list of separated vectors
#' @examples
#' myString <- 'a;b;c;d'
#' semiVector(myString)
#' myStrings <- c('a;b;c;d','e;f;g;h)
#' semiVector(myStrings)
#' @export

semiVector <- function(string_to_split){
  
  if(missing(string_to_split)){
    stop("semiVector expects a character argument")
  }
  
  # If NA is passed, return NA
  if(length(string_to_split)==1){
    if(is.na(string_to_split)){
      return(NA)
    }
  }
  
  # 'string_to_split' should be a character vector
  if(!is.character(string_to_split)){
    stop("semiVector expects a character argument")
  }
  
  # If 'string_to_split' is one element, return bits as a character vector
  if(length(string_to_split)==1){
    return(stringr::str_split(string_to_split[1],pattern = ';') %>% unlist())
  }
  
  # If 'string_to_split' is > one element, return character vectors as a list
  if(length(string_to_split)>1){
    sep_list = purrr::map(.x=string_to_split,.f=function(x){stringr::str_split(x,pattern = ';')[[1]]})
    return(sep_list)
  }
}