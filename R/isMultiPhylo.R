#' Rboretum MultiPhylo Checker
#'
#' This function returns TRUE if the passed object is of class multiPhylo and has 2+ trees; Otherwise, FALSE
#' @param test_object R object to check
#' @param check_named OPTIONAL: If TRUE, also check if trees are named. [Default: FALSE, don't check names]
#' @param check_rooted OPTIONAL: If TRUE, also check if all trees are rooted. [Default: FALSE, don't check rootedness]
#' @param check_three_taxa OPTIONAL: If TRUE, also check if all trees share at least three taxa. [Default: FALSE, don't check taxa]
#' @param check_all_taxa OPTIONAL: If TRUE, also check if all trees share all taxa. [Default: FALSE, don't check taxa]
#' @return TRUE if 'test_object' is a multiPhylo with 2+ trees; otherwise, FALSE
#' @examples 
#' isMultiPhylo(myTrees) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees
#' isMultiPhylo(myTrees,check_named=TRUE) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees have names
#' isMultiPhylo(myTrees,check_rooted=TRUE) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees are rooted
#' isMultiPhylo(myTrees,check_three_taxa) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees share at least three taxa
#' isMultiPhylo(myTrees,check_all_taxa) # Check if 'myTrees' is a valid multiPhylo object with 2+ trees, and that all trees share all taxa
#' @export

isMultiPhylo <- function(test_object,check_named,check_rooted,check_three_taxa,check_all_taxa){

  # Check if trees are named?
  if(missing(check_named)){
    check_named <- FALSE
  } else if(!is.logical(check_named)){
    check_named <- FALSE
  }
  
  # Check if all trees are rooted?
  if(missing(check_rooted)){
    check_rooted <- FALSE
  } else if(!is.logical(check_rooted)){
    check_rooted <- FALSE
  }
  
  # Check if all trees share at least three taxa?
  if(missing(check_three_taxa)){
    check_three_taxa <- FALSE
  } else if(!is.logical(check_three_taxa)){
    check_three_taxa <- FALSE
  }
  
  # Check if all trees share all taxa?
  if(missing(check_all_taxa)){
    check_all_taxa <- FALSE
  } else if(!is.logical(check_all_taxa)){
    check_all_taxa <- FALSE
  }
  
  if(has_error(silent=TRUE,expr=unlist(attributes(test_object)))){ # Can attributes be unlisted?
    return(FALSE) # Object attributes can't be unlisted --> FALSE
  } else if(!'multiPhylo' %in% unlist(attributes(test_object)$class)){
      return(FALSE) # 'multiPhylo' not in object class --> FALSE
    } else if(length(test_object) < 2){
      return(FALSE) # 'test_object' has length < 2 --> FALSE
    } else if(!purrr::map(.x = test_object,.f = function(x){Rboretum::isPhylo(x)}) %>% unlist() %>% all()){
      return(FALSE) # 'test_object' not made up of phylo objects --> FALSE
    } else{
      
      # Summarize taxa in multiPhylo 'test_object'
      tree_taxa <- purrr::map(.x = test_object,.f = function(x){x$tip.label}) %>% unlist() %>% unique() %>% sort()
      taxa_count <- length(tree_taxa)
      tree_count <- length(test_object)
      
      # Get tips from each tree, table results, and count tip labels that occur in all trees
      shared_count <- purrr::map(.x = test_object,.f = function(x){x$tip.label}) %>% unlist() %>% table() %>% (function(x){x==tree_count}) %>% sum()
      
      has_names <- !is.null(names(test_object))
      if(has_names){
        tree_names <- names(test_object)
        name_length <- length(names(test_object))
        name_error <- any(is.na(tree_names)) | any(is.null(tree_names)) | name_length != tree_count
      }
      
      if(check_named & (!has_names | name_error)){
        return(FALSE) # 'test_object' is a valid multiPhylo but trees are not named --> FALSE
      } else if(check_rooted & !all(unlist(purrr::map(.x = test_object,.f = ape::is.rooted)))){
        return(FALSE) # 'test_object' is a valid multiPhylo but trees are not rooted --> FALSE
    } else if(check_three_taxa & shared_count < 3){
        return(FALSE) # 'test_object' is a valid multiPhylo but trees share fewer than three taxa --> FALSE
      } else if(check_all_taxa & shared_count < taxa_count){
        return(FALSE) # 'test_object' is a valid multiPhylo but trees do not have identical taxa --> FALSE
      } else{ 
        return(TRUE) # All checks passed --> TRUE
      }
    }
}
