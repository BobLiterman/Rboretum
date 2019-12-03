#' Rboretum Clade Support Checker
#'
#' This function returns TRUE if the passed object is the output of getCladeSupport; otherwise, FALSE
#' @param test_object R object to check
#' @param tree OPTIONAL: phylo object [Check that passed clade support contains information about a specific tree.]
#' @return TRUE if output of getCladeSupport(); otherwise, FALSE
#' @export


isCladeSupport <- function(test_object,tree){
    
  if(!is.data.frame(test_object)){
    return(FALSE)
  } else if(nrow(test_object) == 0){
    return(FALSE)
  } else if(has_error(silent=TRUE,expr=all(names(test_object)==c('Clade','Tree_Count','Clade_Size','Tree_Percent','Trees_with_Clade')))){
    return(FALSE)
  } else if(!all(names(test_object)==c('Clade','Tree_Count','Clade_Size','Tree_Percent','Trees_with_Clade'))){    
    return(FALSE)
  } 
  
  if(missing(tree)){
    return(TRUE)
  } else if(!Rboretum::isPhylo(tree)){
    stop("'tree' argument does not appear to be a valid phylo object")
  } else if(!ape::is.rooted(tree)){
    stop("'tree' must be rooted for isCladeSupport()")
  } else{
    
    clades <- test_object %>% pull(Clade) %>% as.character()
    tree_clades <- Rboretum::getTreeClades(tree) %>% as.character()

    if(all(tree_clades %in% clades)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}
