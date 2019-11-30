#' Rboretum Tree Support Checker
#'
#' This function returns TRUE if the passed object is the output of getTreeSupport; otherwise, FALSE
#' @param test_object R object to check
#' @param tree OPTIONAL: phylo object [Check that passed tree support contains information about a specific tree.]
#' @return TRUE if output of getTreeSupport(); otherwise, FALSE
#' @export
#'
isTreeSupport <- function(test_object,tree){
  
  if(!is.data.frame(test_object)){
    return(FALSE)
  } else if(nrow(test_object) == 0){
    return(FALSE)
  } else if(has_error(all(names(existing_splits)[1:4] == c('Clade','Mirror_Clade','Split_Node','Split_Bootstrap')))){
    return(FALSE)
  } else if(!all(names(existing_splits)[1:4] == c('Clade','Mirror_Clade','Split_Node','Split_Bootstrap'))){    
    return(FALSE)
  } 
  
  if(missing(tree)){
    return(TRUE)
  } else if(!Rboretum::isPhylo(tree)){
    stop("'tree' argument does not appear to be a valid phylo object")
  } else if(!ape::is.rooted(tree)){
    stop("'tree' must be rooted for isTreeSupport()")
  } else{
    
    clades <- test_object %>% pull(Clade) %>% as.character() %>% sort()
    mirror_clades <- test_object %>% pull(Mirror_Clade) %>% as.character() %>% sort()
    
    tree_splits <- Rboretum::getTreeSplits(tree) %>% filter(!is.na(Split_Node))
    tree_clades <- tree_splits %>% pull(Clade) %>% as.character() %>% sort()
    tree_mirror_clades <- tree_splits %>% pull(Mirror_Clade) %>% as.character() %>% sort()
    
    if(all(clades == tree_clades) & all(mirror_clades == tree_mirror_clades)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}
