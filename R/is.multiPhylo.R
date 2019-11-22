#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo; Otherwise, FALSE
#' @param test_object R object to check
#' @return TRUE if multiPhylo, otherwise FALSE
#' @export
#'
is.multiPhylo <- function(test_object){
  if(has_error(unlist(attributes(test_object)))){
    return(FALSE)
  } else{
    test_object_class <- unlist(attributes(test_object)$class)

    if('multiPhylo' %in% test_object_class){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}
