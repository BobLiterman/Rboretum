#' Rboretum Table Search
#'
#' This function returns the count of items in an R table, or 0 if item does not occur in table
#' @param search_table R Table (Result of 'table()' call)
#' @param name Name of item you want the count for
#' @return Number of times that value occurs in the table
#' @export
tableCount <- function(search_table,name){
  if(!is.integer(search_table) | is.null(names(search_table))){
    stop("'search_table' is not a named integer list.")
  } else if(name %in% names(search_table)){
    return(as.integer(search_table[name]))
  } else{
    return(as.integer(0))
  }
}
