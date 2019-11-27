#' Rboretum Unique Clade Fetcher
#'
#' This function takes a named multiPhylo object and a focal tree, and returns clades that occur in the focal clade but not all other trees in the multiPhylo
#' @param trees Named multiPhylo object
#' @param focal_tree Name of focal tree
#' @return Character vector of clades that occur in the focal tree, and are missing from at least one other tree in the multiPhylo
#' @export
#'

get.uniqueClades <- function(trees,focal_tree){
  
  if(!Rboretum::is.multiPhylo(trees)){
    stop("'trees' does not appear to be a valid multiPhlyo object with 2+ trees")
  } else if(!Rboretum::check.shared(trees)){
    stop("'trees' don't appear to share at least 3 taxa in common.")
  } else if(Rboretum::sameTopology(trees)){
    stop("'trees' have identical topology. No comparison possible.")
  } else if(is.null(names(trees))){
    stop("'trees' must be named for get.uniqueClades. Name multiPhylo by assigning through names(trees) <- c('name1','name2',etc)")
  }
  
  tree_count <- length(trees)
  
  if(missing(focal_tree)){
    stop("'focal_tree' name not assigned.")
  } else{
    focal_tree <- as.character(focal_tree)
  }
  if(has_error(tree <- trees[[focal_tree]])){
    stop("'focal_tree' must be the name of a tree in the multiPhylo. Check names(trees).")
  }
  
  clade_compare <- Rboretum::compare.clades(trees) %>%
    filter(Tree_Percent < 100)
  
  clades <- c()
  for(i in 1:nrow(clade_compare)){
    
    if(str_detect(clade_compare$Trees_with_Clade[i],";")){
      check <- semiVector(clade_compare$Trees_with_Clade[i])
    } else{ check <- clade_compare$Trees_with_Clade[i] }
    if(focal_tree %in% check){
      clades <- c(clades,clade_compare$Clade[i])
    }
  }
  
  return(clades)
}
