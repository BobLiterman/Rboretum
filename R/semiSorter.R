#' Rboretum Semicolon Sorter
#'
#' This function returns a sorted, semicolon-delimited string from an unsorted semicolon-delimited string
#' @param string_to_split Unsorted semicolon delimited string
#' @return Sorted semicolon delimited string
#' @examples 
#' myShuffledString <- 'c;d;a;b'
#' semiSorter(myShuffledString)
#' > 'a;b;c;d'
#' @export

semiSorter <- function(string_to_sort){
    return(paste(sort(semiVector(string_to_sort)),collapse = ";"))
}
