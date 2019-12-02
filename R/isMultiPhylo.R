#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo and has 2+ trees; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_rooted OPTIONAL: If TRUE, also check if all trees are rooted. Default: FALSE, don't check rootedness
#' @param check_named OPTIONAL: If TRUE, also check if trees are named. Default: FALSE, don't check names
#' @return TRUE if all flagged conditions match, otherwise FALSE
#' @export

isMultiPhylo <- function(test_object,check_rooted,check_names){
  
  if(missing(check_rooted)){
    check_rooted <- FALSE
  } else if(!is.logical(check_rooted)){
    check_rooted <- FALSE
  }
  
  if(missing(check_names)){
    check_names <- FALSE
  } else if(!is.logical(check_names)){
    check_names <- FALSE
  }
  
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){
    return(FALSE)
  } else{
    test_object_class <- unlist(attributes(test_object)$class)

    if('multiPhylo' %in% test_object_class & length(test_object)>=2){
      
      root_check  <- c()
      for(i in 1:length(test_object)){
        root_check <- c(root_check,ape::is.rooted(test_object[[i]]))
      }
      name_check <- is.null(names(test_object))
      
      if(check_names & check_rooted){
        if(all(root_check) & !name_check){
          return(TRUE)
        } else{
          return(FALSE)
        }
      } else if(check_names & !check_rooted){
        if(!name_check){
          return(TRUE)
        } else{
          return(FALSE)
        }
      } else if(!check_names & check_rooted){
        if(all(root_check)){
          return(TRUE)
        } else{
          return(FALSE)
        }
      } else if(!check_names & !check_rooted){
        return(TRUE)
      }
    }
    
  }
}