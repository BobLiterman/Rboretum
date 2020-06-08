#' Rboretum Phylo Checker
#'
#' This function returns TRUE if the passed object is of class phylo; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_rooted OPTIONAL: If TRUE, also check if tree is rooted. [Default: FALSE, don't check rootedness]
#' @return TRUE if phylo, otherwise FALSE
#' @examples
#' isPhylo(tree) # Check if 'tree' is a phylo object
#' isPhylo(tree,check_rooted = TRUE) # Check if 'tree' is a rooted phylo object
#' @export

isPhylo <- function(test_object,check_rooted){
  
  if(missing(check_rooted)){
    check_rooted <- FALSE
  } else if(!is.logical(check_rooted)){
    check_rooted <- FALSE
  }
  
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){ # Can attributes be unlisted?
    return(FALSE) # Object attributes can't be unlisted --> FALSE
  } else{

    if('phylo' %in% unlist(attributes(test_object))['class']){ # Is 'phylo' in $class?
      if(!check_rooted){
        return(TRUE) # Is a tree --> TRUE
      } else{
        if(ape::is.rooted(test_object)){ # Is tree rooted?
          return(TRUE) # Is a tree, and rooted --> TRUE
        } else{
          return(FALSE) # Is a tree, but not rooted --> FALSE
        }
      }
    } else{
      return(FALSE) # 'phylo' not in $class attributes --> FALSE
    }
  }
}
