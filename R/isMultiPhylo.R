#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo and has 2+ trees; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_rooted OPTIONAL: If TRUE, also check if tree is rooted. Default: FALSE, don't check rootedness
#' @return TRUE if multiPhylo with 2+ trees, otherwise FALSE
#' @export

isMultiPhylo <- function(test_object,check_rooted){
  
  if(missing(check_rooted)){
    check_rooted <- FALSE
  } else if(!is.logical(check_rooted)){
    check_rooted <- FALSE
  }
  
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){
    return(FALSE)
  } else{
    test_object_class <- unlist(attributes(test_object)$class)

    if('multiPhylo' %in% test_object_class & length(test_object)>=2){
      if(!check_rooted){
        return(TRUE)
      } else{
        root_check  <- c()
        for(i in 1:length(test_object)){
          root_check <- c(root_check,ape::is.rooted(test_object[[i]]))
        }
        if(all(root_check)){
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
