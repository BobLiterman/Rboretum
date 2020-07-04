#' Rboretum Semicolon Vectorizer
#'
#' This function takes a string separated by semicolons (;) and returns a character vector
#' @param string_to_split Semicolon delimited character
#' @return If string_to_split is a single item, semiVector returns a character vector split by ';'. If string_to_split contains multiple items, function returns a list of separated vectors
#' @examples
#' myString <- 'a;b;c;d'
#' semiVector(myString)
#' > ["a","b","c","d"]
#' myStrings <- c('a;b;c;d','e;f;g;h)
#' semiVector(myStrings)
#' > [[1]]
#' > [1] ["a","b","c","d"]
#' > [[2]]
#' > [1] ["e","f","g","h"]
#' @export

semiVector <- function(string_to_split){
  
  if(missing(string_to_split)){
    stop("semiVector expects a character argument")
  } else if(!is.character(string_to_split)){
    if(length(string_to_split)>1){
      stop("semiVector expects a character argument")
    } else{
      if(is.na(string_to_split)){
        return(NA)
      } else{
        stop("semiVector expects a character argument")
      } 
    }
  }
  
  if(length(string_to_split)==1){
      return(stringr::str_split(string_to_split[1],pattern = ';') %>% unlist())
  } else{
    sep_list = purrr::map(.x=string_to_split,.f=function(x){stringr::str_split(x,pattern = ';')[[1]]})
    return(sep_list)
  }
}