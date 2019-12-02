#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo and has 2+ trees; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_rooted OPTIONAL: If TRUE, also check if all trees are rooted. Default: FALSE, don't check rootedness
#' @param check_named OPTIONAL: If TRUE, also check if trees are named. Default: FALSE, don't check names
#' @param check_taxa OPTIONAL: If TRUE, also check if trees share identical taxa. Default: FALSE, don't check taxa
#' @param check_unique OPTIONAL: If TRUE, also check if trees all have a unique topology. Default: FALSE, don't check topology
#' @return TRUE if all flagged conditions match, otherwise FALSE
#' @export

isMultiPhylo <- function(test_object,check_rooted,check_names,check_taxa,check_unique){
  
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
  
  if(missing(check_taxa)){
    check_taxa <- FALSE
  } else if(!is.logical(check_taxa)){
    check_taxa <- FALSE
  }
  
  if(missing(check_unique)){
    check_unique <- FALSE
  } else if(!is.logical(check_unique)){
    check_unique <- FALSE
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
      
      root_check <- all(root_check)
      name_check <- !is.null(names(test_object))
      taxa_check <- Rboretum::checkSameTaxa(test_object)
      
      if(has_error(silent = TRUE,Rboretum::compareTrees(test_object))){
        unique_check <- FALSE
      } else{
      
        compare_vector <- Rboretum::compareTrees(test_object) %>% pull(Comparable)
      
        if(all(compare_vector)){
          unique_check <- TRUE
        } else{
          unique_check <- FALSE
        }
      }
      
      if(!check_rooted){
        root_check <- TRUE
      }
      if(!check_names){
        name_check <- TRUE
      }
      if(!check_taxa){
        taxa_check <- TRUE
      }
      if(!check_unique){
        unique_check  <- TRUE
      }
      if(all(c(root_check,name_check,taxa_check,unique_check))){
        return(TRUE)
      } else{
        return(FALSE)
      }
    } else{
      return(FALSE)
    }
  }
}