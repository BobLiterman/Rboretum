#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo and has 2+ trees; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_named OPTIONAL: If TRUE, also check if trees are named. [Default: FALSE, don't check names]
#' @param check_rooted OPTIONAL: If TRUE, also check if all trees are rooted. [Default: FALSE, don't check rootedness]
#' @param check_three_taxa OPTIONAL: If TRUE, also check if all trees share at least three taxa. [Default: FALSE, don't check taxa]
#' @param check_all_taxa OPTIONAL: If TRUE, also check if all trees share all taxa. [Default: FALSE, don't check taxa]
#' @param verbose OPTIONAL: If TRUE: Print messages along the way; else execute silently. [Default: FALSE]
#' @return TRUE if 'test_object' is a multiPhylo with 2+ trees; otherwise, FALSE
#' @export

isMultiPhylo <- function(test_object,check_named,check_rooted,check_three_taxa,check_all_taxa,verbose){
  
  if(missing(verbose)){
    verbose <- FALSE
  } else if(!is.logical(verbose)){
    verbose <- FALSE
  }
  
  if(missing(check_named)){
    check_named <- FALSE
  } else if(!is.logical(check_named)){
    check_named <- FALSE
  }
  
  if(missing(check_rooted)){
    check_rooted <- FALSE
  } else if(!is.logical(check_rooted)){
    check_rooted <- FALSE
  }
  
  if(missing(check_three_taxa)){
    check_three_taxa <- FALSE
  } else if(!is.logical(check_three_taxa)){
    check_three_taxa <- FALSE
  }
  
  if(missing(check_all_taxa)){
    check_all_taxa <- FALSE
  } else if(!is.logical(check_all_taxa)){
    check_all_taxa <- FALSE
  }
  
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){
    if(verbose){ print("Can't execute unlist(attributes(test_object))")}
    return(FALSE)
  } else{
    
    test_object_class <- unlist(attributes(test_object)$class)
    
    if(!'multiPhylo' %in% test_object_class){
      if(verbose){ print("'multiPhylo' not an attribute of test_object")}
      return(FALSE)
    } else if(length(test_object) < 2){
      if(verbose){ print("'multiPhylo' contains fewer than two trees.")}
      return(FALSE)
    } else if(check_named & is.null(names(test_object))){
      if(verbose){ print("'test_object' has no names assigned.")}
      return(FALSE)
    } else if(check_rooted & !all(unlist(purrr::map(.x = test_object,.f = ape::is.rooted)))){
      if(verbose){ print("One or more trees from 'test_object' are not rooted.")}
      return(FALSE)
    } else if(check_three_taxa & !purrr::map(.x = test_object,.f = function(x){return(x$tip.label)}) %>% unlist() %>% table() %>% (function(x){x==length(test_object)}) %>% sum() >= 3){
      if(verbose){ print("Trees share fewer than three common taxa.")}
      return(FALSE)
    } else if(check_all_taxa & !purrr::map(.x = test_object,.f = function(x){return(x$tip.label)}) %>% unlist() %>% table() %>% (function(x){x==length(test_object)}) %>% sum() == purrr::map(.x = test_object,.f = function(x){return(x$tip.label)}) %>% unlist() %>% unique() %>% length()){
      if(verbose){ print("Trees do not share all taxa.")}
      return(FALSE)
    } else{
      return(TRUE)
    }
  }
}
