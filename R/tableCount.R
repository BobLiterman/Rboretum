#' Search Table for Value
#'
#' This function returns the count of items in a table, or 0 if item does not occur in table
#' @param table_to_search Table (Result of 'table' call)
#' @param name_to_search Title of item you want the count for
#' @return Number of times that value occurs in the table
#' @export
#' @examples
#' tableCount(table_to_search,name_to_search)
#'

tableCount <- function(table_to_search,name_to_search){
  if(name_to_search %in% names(table_to_search)){
    return(as.integer(table_to_search[name_to_search]))
  }
  else{
    return(as.integer(0))
  }
}
