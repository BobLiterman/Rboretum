#' Rboretum Plot Checker
#'
#' This function returns TRUE if the passed object is of class ggtree or ggplot; Otherwise, FALSE
#' @param test_object R object to check
#' @return TRUE if plot, otherwise FALSE
#' @export
#'
isPlot <- function(test_object){
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){
    return(FALSE)
  } else{
    test_object_class <- unlist(attributes(test_object)$class)

    if(!any(c('ggtree','ggplot') %in% test_object_class)){
      return(FALSE)
    } else{
      return(TRUE)
    }
  }
}
