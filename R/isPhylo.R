#' Rboretum Phylo Checker
#'
#' This function returns TRUE if the passed object is of class phylo; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_rooted OPTIONAL: If TRUE, also check if tree is rooted. Default: FALSE, don't check rootedness
#' @return TRUE if phylo, otherwise FALSE
#' @export

isPhylo <- function(test_object,check_rooted){
  
  if(missing(check_rooted)){
    check_rooted <- FALSE
  } else if(!is.logical(check_rooted)){
    check_rooted <- FALSE
  }
  
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){
    return(FALSE)
  } else{
    test_object_class <- unlist(attributes(test_object)$class)

    if('phylo' %in% test_object_class){
      if(!check_rooted){
        return(TRUE)
      } else{
        if(ape::is.rooted(test_object)){
          return(TRUE)
        } else{
          return(FALSE)
        }
      }
    } else{
      return(FALSE)
    }
  }
}
