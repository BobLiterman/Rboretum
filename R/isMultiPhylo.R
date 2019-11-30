#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo and has 2+ trees; Otherwise, FALSE
#' @param test_object R object to check
#' @return TRUE if multiPhylo with 2+ trees, otherwise FALSE
#' @export
#'
isMultiPhylo <- function(test_object){
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){
    return(FALSE)
  } else{
    test_object_class <- unlist(attributes(test_object)$class)

    if('multiPhylo' %in% test_object_class & length(test_object)>=2){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}
